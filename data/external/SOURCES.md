# Upstream sources

Nothing here is vendored yet. This file is the provenance contract: anything landing in
`data/external/` gets a row, with version, licence and checksum, before it is used.

## Candidate pair sources

### StereoAtlas — primary candidate
- Thangaraj & Murugan, *Scientific Data* 2026. Preprint: `10.26434/chemrxiv.15010093/v1`
- Data: Zenodo `10.5281/zenodo.23158759` (v1, 2026-10-05, 364.9 MB, 18 files) — **CC BY 4.0**
- Code: https://github.com/neurocare-iiitd/StereoAtlas — MIT
- 155,526 isomer families from ChEMBL37; 420,721 compounds; 4,519 targets
- Only two files are needed: `FAMILY_MINMAX_FOLD.csv` (10.7 MB, primary analysis file,
  carries `target_chembl_id` and `stereoselectivity_class`) and `STEREOATLAS_MASTER_LITE.csv`
  (31.2 MB, pre-filtered to defined stereo class and fold ≥2×). ~42 MB, not 365 MB.

**Caveats, all of which require work downstream:**
- Pairs are grouped by same target + organism + **metric** — *not* same assay. Re-joining to
  ChEMBL for `assay_id` is required.
- **No enantiomer/diastereomer flag.** "Isomer families" explicitly includes partial inversion,
  so `stereo.classify_pair` has to run over every candidate regardless.
- Column schema is not published; `curate` must be written against the real headers.
- Reported maximum fold difference is 195,000,000×, which is almost certainly a censored
  "inactive at max concentration" value treated as a measurement. Qualifier filtering is mandatory.

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
