"""Turn StereoAtlas candidates into verified pairs via ChEMBL.

Three filters, cheapest first, because each one is a network round trip we
would rather not spend:

1. SMILES (batched) -> `classify_pair` must say `enantiomers`. Kills meso
   compounds and diastereomers that the descriptor-string test waves through.
2. Atom count -> Boltz-2's affinity head takes at most 128 atoms.
3. Activities (two calls per surviving pair) -> both members must appear in
   ONE assay, uncensored.

Writes a table matching data/pairs/SCHEMA.md.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from stereobench.boltz_io import (
    AFFINITY_MAX_ATOMS,
    AFFINITY_RECOMMENDED_MAX_ATOMS,
    affinity_atom_count,
)
from stereobench.chembl import fetch_activities, fetch_smiles, find_shared_assay
from stereobench.curate import looks_like_a_cutoff, pact_from_nm
from stereobench.stereo import classify_pair


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--candidates", default="data/pairs/candidates.csv", type=Path)
    ap.add_argument("--out", default="data/pairs/pilot.csv", type=Path)
    ap.add_argument("--rejects", default="data/pairs/rejected.csv", type=Path)
    ap.add_argument("--per-group", type=int, default=150, help="candidates to attempt per arm")
    ap.add_argument("--max-stereocenters", type=int, default=2)
    args = ap.parse_args()

    candidates = pd.read_csv(args.candidates)
    pool = candidates[
        (candidates.n_stereocenters <= args.max_stereocenters) & candidates.is_human
    ]
    shortlist = pd.concat(
        [group.head(args.per_group) for _, group in pool.groupby("group")]
    ).reset_index(drop=True)
    print(f"attempting {len(shortlist)} candidates", flush=True)

    ids = sorted(
        set(shortlist.eutomer_chembl_id) | set(shortlist.distomer_chembl_id)
    )
    print(f"fetching SMILES for {len(ids)} compounds...", flush=True)
    smiles = fetch_smiles(ids)
    print(f"  got {len(smiles)}", flush=True)

    rows, rejects = [], []

    def reject(row, reason: str) -> None:
        # No cutoff flag here: a rejected row has no agreed measurement to
        # flag. An earlier version read `eu_value`/`di_value` from the
        # enclosing scope, which are bound further down the loop -- so the
        # first rejection raised NameError and every later one silently
        # recorded the PREVIOUS pair's values.
        rejects.append(
            {
                "eutomer_chembl_id": row.eutomer_chembl_id,
                "distomer_chembl_id": row.distomer_chembl_id,
                "target_chembl_id": row.target_chembl_id,
                "group": row.group,
                "reason": reason,
            }
        )

    # Stage 1 and 2: structure filters, no network beyond the batch above.
    structural = []
    for row in shortlist.itertuples():
        eu, di = smiles.get(row.eutomer_chembl_id), smiles.get(row.distomer_chembl_id)
        if not eu or not di:
            reject(row, "no_smiles")
            continue
        try:
            relationship = classify_pair(eu, di)
        except ValueError:
            reject(row, "unparseable_smiles")
            continue
        if relationship != "enantiomers":
            reject(row, f"not_enantiomers:{relationship}")
            continue
        try:
            n_atoms = max(affinity_atom_count(eu), affinity_atom_count(di))
        except ValueError:
            reject(row, "unparseable_smiles")
            continue
        if n_atoms > AFFINITY_MAX_ATOMS:
            reject(row, f"too_many_atoms:{n_atoms}")
            continue
        structural.append((row, eu, di, n_atoms))

    print(
        f"passed structure filters: {len(structural)} "
        f"(of {len(shortlist)})",
        flush=True,
    )

    # Stage 3: the same-assay requirement. Two calls per pair.
    for i, (row, eu_smiles, di_smiles, n_atoms) in enumerate(structural, 1):
        if i % 25 == 0:
            print(f"  assay lookup {i}/{len(structural)}", flush=True)
        try:
            eu_acts = fetch_activities(
                row.eutomer_chembl_id, row.target_chembl_id, row.activity_type
            )
            di_acts = fetch_activities(
                row.distomer_chembl_id, row.target_chembl_id, row.activity_type
            )
        except RuntimeError as exc:
            print(f"  lookup failed: {exc}", file=sys.stderr)
            reject(row, "lookup_failed")
            continue

        # Distinguish "ChEMBL returned nothing" from "nothing matched". The
        # first is usually transient and recoverable; conflating them hides a
        # shrinking pair set behind a plausible-looking reason.
        if not eu_acts or not di_acts:
            reject(row, "no_activities_returned")
            continue

        shared = find_shared_assay(eu_acts, di_acts)
        if shared is None:
            reject(row, "no_shared_uncensored_assay")
            continue
        eu_m, di_m = shared.eutomer, shared.distomer

        # Potency is decided by the shared assay, not by the archive's mean,
        # so the eutomer may swap sides here. Lower nM is stronger.
        eu_pact, di_pact = (
            pact_from_nm(eu_m.standard_value),
            pact_from_nm(di_m.standard_value),
        )
        if di_pact > eu_pact:
            eu_smiles, di_smiles = di_smiles, eu_smiles
            eu_pact, di_pact = di_pact, eu_pact
            eutomer_id, distomer_id = row.distomer_chembl_id, row.eutomer_chembl_id
            eu_value, di_value = di_m.standard_value, eu_m.standard_value
        else:
            eutomer_id, distomer_id = row.eutomer_chembl_id, row.distomer_chembl_id
            eu_value, di_value = eu_m.standard_value, di_m.standard_value

        rows.append(
            {
                "pair_id": f"{eutomer_id}__{distomer_id}__{eu_m.assay_chembl_id}",
                "relationship": "enantiomers",
                "target_chembl_id": row.target_chembl_id,
                "target_name": row.target_name,
                "assay_chembl_id": eu_m.assay_chembl_id,
                "activity_type": row.activity_type,
                "eutomer_smiles": eu_smiles,
                "distomer_smiles": di_smiles,
                "eutomer_chembl_id": eutomer_id,
                "distomer_chembl_id": distomer_id,
                "eutomer_value_nm": eu_value,
                "distomer_value_nm": di_value,
                "eutomer_pact": eu_pact,
                "distomer_pact": di_pact,
                "delta_pact": eu_pact - di_pact,
                "n_stereocenters": row.n_stereocenters,
                "n_atoms": n_atoms,
                "within_recommended_atoms": n_atoms <= AFFINITY_RECOMMENDED_MAX_ATOMS,
                "group": row.group,
                # Flagged, not dropped: an "inactive at top concentration"
                # reading often reaches ChEMBL with relation "=", so the
                # qualifier filter cannot catch it.
                "suspicious_cutoff": looks_like_a_cutoff(eu_value)
                or looks_like_a_cutoff(di_value),
                "document_chembl_id": eu_m.document_chembl_id,
                # How much choosing was involved. 1 shared assay means the
                # pair was not selected at all; several means it was, and the
                # reader should be able to tell those apart.
                "n_shared_assays": shared.n_shared_assays,
                "n_measurements": shared.n_measurements,
            }
        )

    pairs = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(args.out, index=False)
    pd.DataFrame(rejects).to_csv(args.rejects, index=False)

    print(f"\nverified pairs: {len(pairs)}")
    if not pairs.empty:
        # Regroup on the measured delta: the archive's tier came from a mean
        # across assays, and a same-assay delta can land elsewhere entirely.
        pairs["measured_group"] = pd.cut(
            pairs.delta_pact,
            [-0.01, 0.3, 1.0, 99],
            labels=["near_zero", "middle", "tail"],
        )
        print(pairs.measured_group.value_counts().to_string())
        print(f"\ndelta_pact: median {pairs.delta_pact.median():.2f}, "
              f"max {pairs.delta_pact.max():.2f}")
        print(f"targets: {pairs.target_chembl_id.nunique()}")
        print(f"within recommended 56 atoms: {pairs.within_recommended_atoms.sum()}")
        print(f"flagged on a cutoff value  : {pairs.suspicious_cutoff.sum()}")
        near_two = (pairs.delta_pact.sub(0.30103).abs() < 0.01).sum()
        print(f"within 0.01 of exactly 2x  : {near_two}"
              "   <- a large count means the arms were sampled at a tier edge")
        chose = (pairs.n_shared_assays > 1).sum()
        print(f"had >1 candidate assay     : {chose} "
              f"({chose / len(pairs):.0%}) <- the only rows the assay-choice "
              "rule could have biased")
    print(f"\nrejections: {len(rejects)}")
    if rejects:
        print(pd.DataFrame(rejects).reason.str.split(":").str[0].value_counts().to_string())
    print(f"\nwrote {args.out} and {args.rejects}")


if __name__ == "__main__":
    main()
