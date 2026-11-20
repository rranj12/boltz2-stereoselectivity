# Pair-set column contract

One row per pair. `pAct` columns are −log10(molar), so **higher = stronger**, matching
`boltz_io.Affinity.pic50`.

| Column | Notes |
|---|---|
| `pair_id` | stable identifier |
| `relationship` | from `stereo.classify_pair`: `enantiomers` \| `diastereomers` |
| `target_chembl_id`, `target_uniprot` | |
| `assay_chembl_id` | **identical for both members** — the same-assay requirement |
| `activity_type` | `Ki` \| `IC50` \| `Kd`; not mixed within a pair |
| `eutomer_smiles`, `distomer_smiles` | isomeric SMILES, stereo explicit |
| `eutomer_chembl_id`, `distomer_chembl_id` | |
| `eutomer_pact`, `distomer_pact` | eutomer is the more potent member by construction |
| `delta_pact` | `eutomer_pact − distomer_pact`, ≥ 0 |
| `n_stereocenters` | |
| `group` | `tail` (large Δ) \| `near_zero` (model should be flat) \| `control` |

Excluded at curation time: censored/qualified values (`>`, `<`), pairs spanning assays,
ligands over Boltz-2's affinity atom limit, and anything `classify_pair` does not confirm.
