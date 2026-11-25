"""Materialise a Boltz-2 run: YAMLs, substituted mols, and a manifest.

Builds the crossover described in the README's *Pair-set design*. For each pair
it emits four arms (two graphs x two reference-conformer handednesses) and one
input per replicate seed.

Two modes:

  --mode native       plain SMILES YAMLs. Boltz embeds its own conformer, so
                      only `a_native` and `b_native` exist. This is the
                      baseline arm and it runs on a stock install.
  --mode crossover    all four arms, with the reference conformer supplied by
                      substituting a pickled mol. Needs the substitution route
                      verified first -- see scripts/verify_substitution.py.

The substitution route works by writing `<CODE>.pkl` into Boltz's molecule
cache directory and referencing `<CODE>` as a `ccd` entry in the YAML. Boltz
resolves a ccd code against its canonical dictionary first and only then falls
back to the cache directory (`parse/schema.py:628-643`), so codes here are five
characters and prefixed to make a collision with a real CCD entry unlikely. A
collision would silently substitute the WRONG molecule, so
`verify_substitution.py` checks that what comes back is what went in.

Nothing here runs Boltz or needs a GPU.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from stereobench.boltz_io import build_yaml
from stereobench.conformer import (
    build_crossover,
    build_reference_conformer,
    conformer_labels,
    graph_labels,
    write_mol_pickle,
)

LIGAND_CODE_PREFIX = "SB"


def ligand_code(index: int) -> str:
    """Five-character pseudo-CCD code, e.g. SB042."""
    return f"{LIGAND_CODE_PREFIX}{index:03d}"


def build_ccd_yaml(sequence: str, code: str, *, msa_path: str | None = None) -> str:
    """A YAML whose ligand is a cache-resolved code rather than a SMILES.

    Mirrors `build_yaml`, but the ligand entity carries `ccd:` so Boltz loads
    our mol -- and therefore our reference conformer -- instead of embedding
    one from a SMILES string.
    """
    lines = [
        "version: 1",
        "sequences:",
        "  - protein:",
        "      id: A",
        f"      sequence: {sequence}",
    ]
    if msa_path is not None:
        lines.append(f"      msa: {msa_path}")
    lines += [
        "  - ligand:",
        "      id: B",
        f"      ccd: {code}",
        "properties:",
        "  - affinity:",
        "      binder: B",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pairs", default="data/pairs/pilot.csv", type=Path)
    ap.add_argument("--sequences", default="data/pairs/target_sequences.json", type=Path)
    ap.add_argument("--out", default="runs/pilot", type=Path)
    ap.add_argument("--mode", choices=("native", "crossover"), default="crossover")
    ap.add_argument("--replicates", type=int, default=3)
    ap.add_argument("--max-pairs", type=int, default=12)
    ap.add_argument("--embed-seed", type=int, default=0xC0FFEE)
    args = ap.parse_args()

    pairs = pd.read_csv(args.pairs)
    sequences = json.loads(args.sequences.read_text()) if args.sequences.exists() else {}

    usable = pairs[pairs.target_chembl_id.isin(sequences)].head(args.max_pairs)
    skipped = len(pairs) - len(pairs[pairs.target_chembl_id.isin(sequences)])
    print(f"pairs with a sequence: {len(usable)} (skipped {skipped} without one)")

    yaml_dir = args.out / "inputs"
    mol_dir = args.out / "mols"
    yaml_dir.mkdir(parents=True, exist_ok=True)
    if args.mode == "crossover":
        mol_dir.mkdir(parents=True, exist_ok=True)

    manifest, code_index = [], 0

    for pair in usable.itertuples():
        sequence = sequences[pair.target_chembl_id]
        try:
            arms = build_crossover(
                pair.eutomer_smiles, pair.distomer_smiles, seed=args.embed_seed
            )
        except (ValueError, RuntimeError) as exc:
            print(f"  skipping {pair.pair_id}: {exc}")
            continue

        if args.mode == "native":
            arms = [a for a in arms if not a.is_crossed]

        for arm in arms:
            code = ligand_code(code_index)
            code_index += 1

            if args.mode == "crossover":
                write_mol_pickle(arm.mol, mol_dir / f"{code}.pkl")
                content = build_ccd_yaml(sequence, code)
            else:
                content = build_yaml(sequence, arm.smiles)

            for replicate in range(args.replicates):
                # Boltz's --seed covers both the diffusion noise and the
                # featuriser's roto-translation, so replicates vary it rather
                # than rerunning unseeded. Same quantity, but attributable.
                name = f"{pair.pair_id[:40]}__{arm.arm_id}__r{replicate}".replace("/", "_")
                path = yaml_dir / f"{name}.yaml"
                path.write_text(content)
                manifest.append(
                    {
                        "run_name": name,
                        "pair_id": pair.pair_id,
                        "arm_id": arm.arm_id,
                        "replicate": replicate,
                        "boltz_seed": 1000 + replicate,
                        "ligand_code": code if args.mode == "crossover" else "",
                        "graph_source": arm.graph_source,
                        "conformer_source": arm.conformer_source,
                        "is_crossed": arm.is_crossed,
                        "graph_labels": str(graph_labels(arm.mol)),
                        "conformer_labels": str(conformer_labels(arm.mol)),
                        "target_chembl_id": pair.target_chembl_id,
                        "assay_chembl_id": pair.assay_chembl_id,
                        "eutomer_pact": pair.eutomer_pact,
                        "distomer_pact": pair.distomer_pact,
                        "delta_pact": pair.delta_pact,
                        "group": pair.group,
                        "yaml": str(path),
                    }
                )

    manifest_df = pd.DataFrame(manifest)
    manifest_path = args.out / "manifest.csv"
    manifest_df.to_csv(manifest_path, index=False)

    print(f"\ninputs written : {len(manifest_df)}")
    print(f"distinct pairs : {manifest_df.pair_id.nunique() if len(manifest_df) else 0}")
    print(f"arms per pair  : {manifest_df.arm_id.nunique() if len(manifest_df) else 0}")
    print(f"manifest       : {manifest_path}")

    if args.mode == "crossover":
        print(
            "\nBefore running, copy the substituted mols into Boltz's cache:\n"
            f"  cp {mol_dir}/*.pkl ~/.boltz/mols/\n"
            "and confirm the route first with scripts/verify_substitution.py."
        )
    print(
        "\nRun with --use_msa_server, the SAME msa for both members of a pair, "
        "and note that steering (--use_potentials) is OFF by default."
    )


if __name__ == "__main__":
    main()
