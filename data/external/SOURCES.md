# Upstream sources

This file is the provenance contract: anything landing in `data/external/` gets a row, with
version, licence and checksum, before it is used.

## Vendored

Downloaded 2026-11-24 from Zenodo record 23158759, MD5s verified against the record manifest.
Lives in `data/external/stereoatlas/` (gitignored).

| File | Bytes | MD5 | Role |
|---|---|---|---|
| `FAMILY_MINMAX_FOLD.csv` | 10,688,200 | `212237060c57bf25ad5a2657587a8fc4` | fold + tier per family×target×metric |
| `TRUE_ISOMER_LITERATURE_CLEAN.csv` | 32,392,970 | `7e999e3b0e5cb0025cf72c32bb1389e1` | **per-member table — the pair members** |
| `STEREOATLAS_MASTER_LITE.csv` | 31,248,228 | `4343ec175714f97a4fda4d58594b77b2` | superset of the fold columns + TSI/DCS |
| `TARGET_PRIORITY_RANKING.csv` | 1,127,593 | `65848a8a43b33ebe78423522a78e624f` | target selection |
| `by_metric.zip` | 8,252,463 | `99233f9709afa2546a73c414a6594599` | per-metric splits; adds pairwise-heatmap views |
| `MIXED_SINGLE_FAMILIES.csv` | 94,762,881 | `5b77ef340eb957dc0f8bf54ecf150eab` | **not needed — see below** |
| `ZENODO_README.md` | 11,250 | `19986b72c53359a3323fc854ab2efa9f` | upstream pipeline doc |

Record metadata: published 2026-10-05, licence `cc-by-4.0`, no version string set, 18 files.
Title as deposited: *StereoAtlas: A Large-Scale Dataset of Chiral Compounds with Distinct
Bioactivities Curated from ChEMBL37*.

## Candidate pair sources

### StereoAtlas — primary candidate
- Thangaraj & Murugan, *Scientific Data* 2026. Preprint: `10.26434/chemrxiv.15010093/v1`
- Data: Zenodo `10.5281/zenodo.23158759` (v1, 2026-10-05, 364.9 MB, 18 files) — **CC BY 4.0**
- Code: https://github.com/neurocare-iiitd/StereoAtlas — MIT
- 155,526 isomer families from ChEMBL37; 420,721 compounds; 4,519 targets
**Schema verified 2026-11-24 by direct inspection. The original “only two files are needed”
note was wrong, and so was a first correction to it — both are superseded by this.**

- `FAMILY_MINMAX_FOLD.csv` (16 cols, 65,780 rows, 34,940 distinct families) and
  `STEREOATLAS_MASTER_LITE.csv` (49 cols, 72,812 rows; a column superset) are **family-level**.
  Each row names its two members only as stereo-descriptor strings in `min_variant` /
  `max_variant` (e.g. `R,R,R,R,R,R,S,S,S,S`). `compound_chembl_id` is a single family
  representative, and `scaffold_smiles` is stereo-stripped. Neither file identifies the two
  members individually.
- **`TRUE_ISOMER_LITERATURE_CLEAN.csv` is the per-member table** (18 cols, 155,526 rows — one
  per stereoform, matching the headline “155,526 true isomer pairs”; 35,597 families). Join
  `FAMILY_MINMAX_FOLD` → this on `stereochemical_family_id` + `stereochemical_variant` to turn
  each `min_variant` / `max_variant` into a `compound_chembl_id`. Measured: all 65,780 fold rows
  have their family present here, and **60,462 resolve both members to distinct ChEMBL ids**.
- **`MIXED_SINGLE_FAMILIES.csv` is the wrong file and should not be used.** Despite carrying a
  `canonical_smiles` column, the upstream README defines it as “114,573 families (n≥2) where
  only one mixed R/S stereoform was tested” — i.e. families with no pair. Its family-id space
  is **disjoint** from the fold tables (153,222 vs 34,940 families, intersection 0). Same for
  `ONLY_R_FAMILIES.csv` / `ONLY_S_FAMILIES.csv`.
- **No StereoAtlas file carries a per-member isomeric SMILES for paired families.** SMILES must
  be fetched from ChEMBL by `compound_chembl_id`. This is the same API pass already required for
  `assay_id` and the activity qualifier, so it costs nothing extra — but it means
  `stereo.classify_pair` cannot run until that pass is done.
- **No `assay_chembl_id` and no activity qualifier in any StereoAtlas file.** The same-assay
  requirement in `data/pairs/SCHEMA.md` is unmet by this archive alone; the ChEMBL re-join is
  mandatory, not an optimisation.
- `by_metric.zip` extracts to 39 CSVs across 13 metrics. The `*_PAIRWISE_HEATMAP.csv` variants
  carry `variant_a`/`value_a`/`variant_b`/`value_b`/`direction`, which is a more convenient
  per-comparison shape than MinMax for the near-zero arm. Still no SMILES.

### Usable-pair census (measured, `FAMILY_MINMAX_FOLD.csv`)

Plausibility filter = both values within 1e-4 – 1e7 nM. Only 161 of 65,780 rows fail it, so the
upstream censoring problem is real but far smaller than the 195,000,000× headline suggests.

| Tier | All | Plausible | + both members resolvable |
|---|---|---|---|
| None (<2×) | 31,575 | 31,538 | 26,229 |
| Mild (2–5×) | 13,882 | 13,851 | 13,851 |
| Moderate (5–10×) | 6,252 | 6,231 | 6,231 |
| Strong (10–50×) | 8,559 | 8,537 | 8,537 |
| Extreme (>50×) | 5,512 | 5,462 | 5,462 |
| **Total** | **65,780** | **65,619** | **60,310** |

Of the 60,310, **32,481 have every stereocentre inverted** between the two variant strings —
enantiomer *candidates* only. The variant string is not sufficient to call enantiomerism
(meso compounds), so `stereo.classify_pair` must still confirm each one on real SMILES.

Arms are not supply-limited: ~14,000 tail candidates (Strong + Extreme) and ~26,000 near-zero
candidates (None), before same-assay and ≤128-atom filtering.

### ChEMBL REST API — for requalification
Used to attach `assay_id`, the activity qualifier and the raw value to the shortlist only.
Avoids a multi-GB ChEMBL dump: a few hundred requests covers the candidates that survive
tail-enrichment.

## Related work — cite, do not depend on

- **Chiral Cliffs** — Schneider, Lewis, Fechner & Ertl, ChemMedChem 2018;13(13):1315–1324,
  `10.1002/cmdc.201700798`. ~3,800 chiral pairs from ChEMBL23, same biochemical binding assay.
  The abstract states dataset, descriptors and code ship with the paper, but **availability of
  the SI was not confirmed** (publisher blocks automated access). Older ChEMBL, binary
  cliff/non-cliff orientation.
- **Ashraf 2026**, *Target-family differences in stereochemical discrimination: a same-assay
  analysis of enantiomer pairs in ChEMBL*, `10.64336/001c.170283`. Reports 18,346 same-assay
  enantiomer pairs from ChEMBL37, median Δ ~3×, 27.7% >10×, 7.5% >100×.
  **Journal of High School Science 10(3):608–640, single author. No data availability statement,
  no data link, CC BY-NC-SA.** Useful for the same-assay framing; the benchmark should not
  depend on its numbers. Re-derive the Δ distribution from data we hold.
- **RIKEN Boltz-2 ChEMBL benchmark** — Zenodo `10.5281/zenodo.17960302`, CC BY 4.0. Boltz-2
  affinity predictions over a ChEMBL34-derived set. Not stereo-focused and ships no structures,
  but `Data_S1_target_aaseq.csv` maps UniProt → PDB → sequence, and the numeric predictions are
  a sanity baseline for "is our harness producing values in the right range".

### Note on combining sources
StereoAtlas and Ashraf 2026 are **both ChEMBL37 derivatives**. Merging them double-counts the
same underlying measurements under two incompatible pair definitions and inflates apparent
sample size without adding independent evidence. Use one archive for candidate rows and
reimplement the same-assay discipline directly.
