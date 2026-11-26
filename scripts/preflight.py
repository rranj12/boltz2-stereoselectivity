"""Check everything that can fail twenty minutes into a GPU session.

Ordered by how expensive the failure is to discover late. Each check prints
PASS, WARN or FAIL and the script exits non-zero if anything FAILs, so it can
gate a run script.

The one that matters most is the CCD collision check. Boltz resolves a `ccd:`
code against its canonical dictionary FIRST and only then falls back to the
cache directory (`parse/schema.py:628-643`). If one of our five-character
codes happens to name a real CCD entry, Boltz silently substitutes a
completely different molecule, every arm of that pair is wrong, and nothing in
the output says so. That is the single worst failure mode in the pipeline
because it produces plausible numbers.

Run with the Boltz environment where possible -- the deeper checks need it:

    PYTHONPATH=src path/to/boltz-venv/bin/python scripts/preflight.py
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import pandas as pd

# Weights, CCD dictionary and the molecule cache. Boltz downloads these on
# first run, and a part-written download is a confusing failure.
APPROXIMATE_DOWNLOAD_GB = 6.0
HEADROOM_GB = 10.0

status = {"fail": 0, "warn": 0}


def report(level: str, check: str, detail: str) -> None:
    if level in status:
        status[level] += 1
    print(f"  {level.upper():4}  {check}: {detail}")


def check_disk(cache: Path) -> None:
    anchor = cache if cache.exists() else cache.parent
    free_gb = shutil.disk_usage(anchor).free / 1e9
    needed = APPROXIMATE_DOWNLOAD_GB + HEADROOM_GB
    if free_gb < APPROXIMATE_DOWNLOAD_GB:
        report("fail", "disk", f"{free_gb:.1f} GB free, weights alone need ~{APPROXIMATE_DOWNLOAD_GB:.0f} GB")
    elif free_gb < needed:
        report("warn", "disk", f"{free_gb:.1f} GB free; weights plus outputs want ~{needed:.0f} GB")
    else:
        report("pass", "disk", f"{free_gb:.1f} GB free")


def check_boltz() -> object | None:
    try:
        import boltz
    except ImportError:
        report("warn", "boltz import", "not importable in this interpreter; deeper checks skipped")
        return None
    report("pass", "boltz import", f"{getattr(boltz, '__version__', 'unknown version')}")
    return boltz


def check_sequences(path: Path) -> dict:
    if not path.exists():
        report("fail", "sequences", f"{path} missing -- run scripts/fetch_sequences.py")
        return {}
    sequences = json.loads(path.read_text())
    report("pass", "sequences", f"{len(sequences)} targets resolved")
    return sequences


def check_manifest(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        report("warn", "manifest", f"{path} missing -- run scripts/build_run.py")
        return None
    manifest = pd.read_csv(path).fillna({"msa": "", "ligand_code": ""})
    report("pass", "manifest", f"{len(manifest)} inputs")
    return manifest


def check_inputs_exist(manifest: pd.DataFrame) -> None:
    missing = [row.yaml for row in manifest.itertuples() if not Path(row.yaml).exists()]
    if missing:
        report("fail", "input files", f"{len(missing)} referenced YAMLs do not exist")
    else:
        report("pass", "input files", "all present")


def check_msas(manifest: pd.DataFrame) -> None:
    declared = [m for m in manifest.msa.astype(str) if m]
    if not declared:
        report(
            "warn",
            "msa",
            "no cached MSAs; arms of a pair may get different alignments "
            "(scripts/fetch_msas.py)",
        )
        return
    missing = {m for m in declared if not Path(m).exists()}
    if missing:
        report("fail", "msa", f"{len(missing)} referenced MSA files do not exist")
    elif len(declared) != len(manifest):
        report("warn", "msa", f"only {len(declared)} of {len(manifest)} inputs carry one")
    else:
        report("pass", "msa", f"{len(set(declared))} cached alignments, shared across arms")


def check_mols_installed(manifest: pd.DataFrame, cache: Path, run_dir: Path) -> None:
    codes = {c for c in manifest.ligand_code.astype(str) if c}
    if not codes:
        report("pass", "substituted mols", "none required (native mode)")
        return
    installed = {p.stem for p in (cache / "mols").glob("*.pkl")} if (cache / "mols").exists() else set()
    absent = codes - installed
    if absent:
        report(
            "fail",
            "substituted mols",
            f"{len(absent)} of {len(codes)} not in {cache / 'mols'} "
            f"-- cp {run_dir / 'mols'}/*.pkl {cache / 'mols'}/",
        )
    else:
        report("pass", "substituted mols", f"{len(codes)} installed")


def check_ccd_collision(manifest: pd.DataFrame, cache: Path) -> None:
    """A code that names a real CCD entry silently substitutes another molecule."""
    codes = {c for c in manifest.ligand_code.astype(str) if c}
    if not codes:
        report("pass", "ccd collision", "no substituted codes to check")
        return
    try:
        from boltz.data.mol import load_canonicals

        canonical = set(load_canonicals(str(cache / "mols")))
    except Exception as exc:  # noqa: BLE001 - any failure here is a skipped check
        report("warn", "ccd collision", f"could not load the CCD dictionary ({exc})")
        return
    collisions = codes & canonical
    if collisions:
        report(
            "fail",
            "ccd collision",
            f"{sorted(collisions)} are real CCD codes -- Boltz would load the "
            "WRONG molecule and the run would look fine",
        )
    else:
        report("pass", "ccd collision", f"{len(codes)} codes, none shadow a CCD entry")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", default="runs/pilot", type=Path)
    ap.add_argument("--sequences", default="data/pairs/target_sequences.json", type=Path)
    ap.add_argument("--cache", default=Path.home() / ".boltz", type=Path)
    args = ap.parse_args()

    print(f"preflight for {args.run}\n")

    check_disk(args.cache)
    check_boltz()
    check_sequences(args.sequences)

    manifest = check_manifest(args.run / "manifest.csv")
    if manifest is not None and len(manifest):
        check_inputs_exist(manifest)
        check_msas(manifest)
        check_mols_installed(manifest, args.cache, args.run)
        check_ccd_collision(manifest, args.cache)

    print(f"\n{status['fail']} failures, {status['warn']} warnings")
    raise SystemExit(1 if status["fail"] else 0)


if __name__ == "__main__":
    main()
