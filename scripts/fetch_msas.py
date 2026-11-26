"""Generate one MSA per target, once, so every arm of a pair shares it.

Why this exists
---------------

`build_run.py` emits four arms x N replicates per pair, all with an identical
protein sequence. Run with `--use_msa_server`, Boltz generates an MSA
separately for every one of those input files. The server is not promised to
be deterministic, so the arms of a pair can end up conditioned on different
alignments.

That is a problem specific to this experiment rather than a general
inefficiency. The design is a paired comparison in which the protein is held
fixed and only the ligand's handedness moves. An MSA that varies between arms
puts variance into exactly the difference being measured, and it is invisible
in the output -- the prediction just comes back slightly different, and nothing
says why.

Generating the MSA once per target and pointing every arm at the same file
removes that degree of freedom entirely.

How
---

By calling Boltz's own `run_mmseqs2` and writing the CSV in Boltz's own
format, rather than reimplementing either. The point is to be byte-identical
to what `--use_msa_server` would have produced for a single-chain input, so
this script is a caching layer and not a second implementation that can drift.
Single-chain only: with one protein entity there is no pairing step, which is
why `main.py` sets `paired_msas` to empty in that case and why this does too.

Needs the Boltz environment, since it imports Boltz:

    PYTHONPATH=src path/to/boltz-venv/bin/python scripts/fetch_msas.py

Network-only. No GPU.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd


def msa_filename(target_chembl_id: str) -> str:
    return f"{target_chembl_id}.csv"


def write_msa_csv(rows: list[str], path: Path) -> None:
    """Boltz's single-chain MSA format: `key,sequence`, key -1 for unpaired."""
    lines = ["key,sequence"] + [f"-1,{seq}" for seq in rows]
    path.write_text("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pairs", default="data/pairs/pilot.csv", type=Path)
    ap.add_argument("--sequences", default="data/pairs/target_sequences.json", type=Path)
    ap.add_argument("--out", default="data/msa", type=Path)
    ap.add_argument(
        "--max-targets",
        type=int,
        default=8,
        help="the MSA server is a free academic service; take only what the "
        "run needs rather than the whole pair set",
    )
    ap.add_argument("--host-url", default="https://api.colabfold.com")
    args = ap.parse_args()

    try:
        from boltz.data import const
        from boltz.data.msa.mmseqs2 import run_mmseqs2
    except ImportError:
        sys.exit(
            "boltz is not importable. Run this with the Boltz environment:\n"
            "  PYTHONPATH=src path/to/boltz-venv/bin/python scripts/fetch_msas.py"
        )

    if not args.sequences.exists():
        sys.exit(f"{args.sequences} does not exist -- run scripts/fetch_sequences.py first")

    sequences = json.loads(args.sequences.read_text())
    pairs = pd.read_csv(args.pairs)

    # Order by the pair set so the targets we fetch are the ones the pilot run
    # will actually use, not an arbitrary slice of 150.
    ordered, seen = [], set()
    for target in pairs.target_chembl_id:
        if target in sequences and target not in seen:
            seen.add(target)
            ordered.append(target)
    wanted = ordered[: args.max_targets]

    args.out.mkdir(parents=True, exist_ok=True)
    print(f"{len(ordered)} targets have a sequence; fetching {len(wanted)}")

    written, skipped = 0, 0
    for i, target in enumerate(wanted, 1):
        destination = args.out / msa_filename(target)
        if destination.exists():
            print(f"  [{i}/{len(wanted)}] {target} already cached")
            skipped += 1
            continue

        sequence = sequences[target]
        print(f"  [{i}/{len(wanted)}] {target} ({len(sequence)} aa) ...", flush=True)

        # Single chain: no pairing step, matching main.py's len(data) == 1 path.
        unpaired = run_mmseqs2(
            [sequence],
            str(args.out / f"{target}_tmp"),
            use_env=True,
            use_pairing=False,
            host_url=args.host_url,
        )

        # a3m comes back as alternating header/sequence lines.
        rows = unpaired[0].strip().splitlines()[1::2][: const.max_msa_seqs]
        if not rows:
            print(f"      no alignment returned, skipping {target}")
            continue
        write_msa_csv(rows, destination)
        written += 1
        print(f"      {len(rows)} sequences -> {destination}")

    print(f"\nwrote {written}, already cached {skipped}, into {args.out}")
    print(
        "Pass --msa-dir to build_run.py so every arm of a pair points at the "
        "same file."
    )


if __name__ == "__main__":
    main()
