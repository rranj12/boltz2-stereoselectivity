"""Build the replicate inputs for the noise-floor measurement.

This is the go/no-go gate and it should run before any other GPU work. The
project's claims are all differences between paired predictions. A difference
is only readable against the scatter of the measurement itself: if repeat
predictions of the SAME input vary by 0.5 pAct, then a near-zero arm differing
by 0.006 cannot distinguish "the model is correctly flat" from "the model is
noise", and no amount of curation rescues that.

Two arms, because there are two different noise floors and the project needs
to know both.

  substituted   the reference conformer is pinned by a cached pickled mol, so
                only Boltz's own stochasticity varies: diffusion noise and the
                featuriser's roto-translation, both covered by --seed. This is
                the floor the crossover experiment is measured against.

  native        plain SMILES, so Boltz embeds its own conformer. ETKDG is
                unseeded inside Boltz (`compute_3d_conformer` never sets
                `options.randomSeed`), so this arm additionally carries
                reference-conformer variance that --seed does not control.

The gap between the two is not overhead -- it is a measurement of how much the
unseeded conformer contributes on its own, which is the quantity the whole
mechanistic arm is about. Running only the substituted arm would leave that
unmeasured.

Nothing here runs Boltz or needs a GPU.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from stereobench.boltz_io import build_yaml
from stereobench.conformer import build_reference_conformer, write_mol_pickle

from build_run import build_ccd_yaml, ligand_code, write_run_script


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pairs", default="data/pairs/pilot.csv", type=Path)
    ap.add_argument("--sequences", default="data/pairs/target_sequences.json", type=Path)
    ap.add_argument("--msa-dir", type=Path, default=None)
    ap.add_argument("--out", default="runs/noise_floor", type=Path)
    ap.add_argument("--ligands", type=int, default=3)
    ap.add_argument("--replicates", type=int, default=8)
    ap.add_argument("--base-seed", type=int, default=2000)
    ap.add_argument("--embed-seed", type=int, default=0xC0FFEE)
    args = ap.parse_args()

    if not args.sequences.exists():
        raise SystemExit(
            f"{args.sequences} does not exist -- run scripts/fetch_sequences.py first"
        )

    sequences = json.loads(args.sequences.read_text())
    pairs = pd.read_csv(args.pairs)
    pairs = pairs[pairs.target_chembl_id.isin(sequences)]

    # Distinct targets: a noise floor measured on three ligands against one
    # protein would describe that protein, not the model.
    chosen = pairs.drop_duplicates(subset="target_chembl_id").head(args.ligands)
    if len(chosen) < args.ligands:
        print(f"only {len(chosen)} distinct targets available")

    yaml_dir, mol_dir = args.out / "inputs", args.out / "mols"
    yaml_dir.mkdir(parents=True, exist_ok=True)
    mol_dir.mkdir(parents=True, exist_ok=True)

    def msa_for(target: str) -> str | None:
        if args.msa_dir is None:
            return None
        candidate = args.msa_dir / f"{target}.csv"
        return str(candidate.resolve()) if candidate.exists() else None

    manifest = []
    for index, pair in enumerate(chosen.itertuples()):
        sequence = sequences[pair.target_chembl_id]
        msa_path = msa_for(pair.target_chembl_id)
        smiles = pair.eutomer_smiles

        mol = build_reference_conformer(smiles, seed=args.embed_seed)
        code = ligand_code(index)
        write_mol_pickle(mol, mol_dir / f"{code}.pkl")

        arms = {
            "substituted": build_ccd_yaml(sequence, code, msa_path=msa_path),
            "native": build_yaml(sequence, smiles, msa_path=msa_path),
        }

        for arm, content in arms.items():
            for replicate in range(args.replicates):
                name = f"nf{index}__{arm}__r{replicate}"
                path = yaml_dir / f"{name}.yaml"
                path.write_text(content)
                manifest.append(
                    {
                        "run_name": name,
                        "ligand_index": index,
                        "arm": arm,
                        "replicate": replicate,
                        # Seeds run in the same range across arms so the two
                        # are compared over the same stochasticity, not over
                        # two different draws.
                        "boltz_seed": args.base_seed + replicate,
                        "ligand_code": code if arm == "substituted" else "",
                        "target_chembl_id": pair.target_chembl_id,
                        "smiles": smiles,
                        "yaml": str(path),
                        "msa": msa_path or "",
                    }
                )

    manifest_df = pd.DataFrame(manifest)
    manifest_path = args.out / "manifest.csv"
    manifest_df.to_csv(manifest_path, index=False)

    write_run_script(
        manifest_df,
        args.out / "run.sh",
        out_root=(args.out / "outputs").resolve(),
        uses_cached_msa=bool(len(manifest_df)) and bool((manifest_df.msa != "").all()),
    )

    print(f"ligands        : {len(chosen)}")
    print(f"inputs written : {len(manifest_df)} "
          f"({len(chosen)} x 2 arms x {args.replicates} replicates)")
    print(f"manifest       : {manifest_path}")
    print(f"run script     : {args.out / 'run.sh'}")
    print(
        f"\nCopy the pinned mols into Boltz's cache first:\n"
        f"  cp {mol_dir}/*.pkl ~/.boltz/mols/\n"
        "\nsigma is the per-arm standard deviation of predicted affinity across\n"
        "replicates. Read it before committing GPU time to anything else: if it\n"
        "is comparable to the effects being chased, the experiment cannot see\n"
        "them and that is worth knowing on day one rather than day nine."
    )


if __name__ == "__main__":
    main()
