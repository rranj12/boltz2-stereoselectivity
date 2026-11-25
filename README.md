# boltz2-stereoselectivity

**Does Boltz-2's affinity head respond to stereochemistry?**

Enantiomers are the cleanest possible probe of a structure-based affinity model.
Identical atoms, identical bonds, identical formula, identical 2D graph, but mirror-image
geometry. Experimentally, one enantiomer frequently binds far more tightly than the other.
If Boltz-2 assigns them the same affinity, the head is reading chemistry and not shape.

There is a second, already-documented failure to attach this to: Boltz-2 often generates the
wrong enantiomer from a stereo-specified SMILES. That allows the paired question — *when the
geometry comes out mirrored, does the predicted affinity follow the geometry or the input
string?*

Existing probes of Boltz-2's coordinate sensitivity have all worked the protein side
(mutate the pocket, shuffle targets) or substituted whole poses. 

## Status

The verification and evaluation core is implemented and tested. The upstream candidate archive
is downloaded and its schema verified. The curated pair set is not built yet, and no inference
has been run.

**Done**

- `stereobench.stereo` / `boltz_io` / `metrics` — implemented, unit-tested.
- StereoAtlas acquired (Zenodo 23158759, CC BY 4.0), MD5s checked against the record manifest,
  provenance recorded in `data/external/SOURCES.md`. Gitignored; re-fetch from the DOI.
- Schema verified against the real headers rather than the paper. Two corrections matter:
  the fold tables are **family-level and do not identify the two members**;
  `TRUE_ISOMER_LITERATURE_CLEAN.csv` is the per-member table, joined on
  `stereochemical_family_id` + `stereochemical_variant`. 60,462 of 65,780 fold rows resolve
  both members to distinct ChEMBL ids.
- Supply measured: ~14,000 tail candidates (Strong + Extreme) and ~26,000 near-zero candidates,
  of which 32,481 are enantiomer candidates pending `classify_pair` confirmation. Upstream
  censoring affects only 161 of 65,780 rows.
- **Boltz source audited** at `v2.2.1` (`b1ebfc4`). Both open questions are answered and a
  third channel was found that the original framing missed — see *Source findings*.

**Blocking, in order**

1. **Affinity-stage coordinate capture.** Confirmed necessary: the affinity head scores
   coordinates it diffuses itself, and nothing writes them to disk. Without a patched Boltz
   there is no checkpoint 3 and the key conditional cannot be evaluated.
2. **ChEMBL requalification pass.** No StereoAtlas file carries a per-member SMILES, an
   `assay_chembl_id` or an activity qualifier. One API pass over the shortlist supplies all
   three. `classify_pair` cannot run until it is done, so there is no `curate` module yet.
3. **Noise floor σ.** Measure before curating the full set; it sets the replicate count for
   every later run.

## What's here

| Module | Does |
|---|---|
| `stereobench.stereo` | CIP labelling from SMILES and from 3D; enantiomer / diastereomer classification; requested-vs-realized comparison |
| `stereobench.boltz_io` | Boltz-2 YAML with an affinity request; parse `affinity_<name>.json`; enforce the affinity ligand-size limits |
| `stereobench.metrics` | Paired metrics, the noise floor, and the key conditional |

Three design points worth knowing, because they are where this kind of benchmark usually
goes wrong:

**Enantiomers are defined exactly, not by counting.** A pair is enantiomeric iff the mirror
image of A is canonically equal to B. The tempting rule — "all stereocenters inverted" —
miscalls meso compounds, which are their own mirror image. `classify_pair` has a test for this.

**CIP labels are keyed by graph-canonical rank**, not atom index. Atom indices are not
comparable between two mols parsed from different SMILES, and the whole benchmark rests on
comparing stereo descriptors across molecules.

**`ambiguous` is a verdict, not a failure.** If a requested stereocenter can't be read off
the predicted geometry, the pair is excluded rather than scored. An unreadable center is
not evidence in either direction.

## Metrics

Per pair, in pAct space (higher = stronger, so predicted and experimental deltas share a
sign convention — `Affinity.pic50` converts):

- **sign accuracy** — does the predicted direction match experiment?
- **Δ correlation** — predicted vs experimental magnitude
- **the key conditional** — when the *realized* geometry is identical between the two
  members, does |Δpredicted| collapse to ~0? A geometry-driven head must be flat here.
  A non-zero Δ means the number came from somewhere other than the coordinates.

### The noise floor is not optional

The headline result is a null, and a null is unreadable without a noise floor. Boltz-2 **does**
expose `--seed` (`main.py:919`), but it **defaults to None — no seeding** (`main.py:1102`), so
out of the box repeated runs of one input differ by some σ. If σ is 0.4 log units and the pairs
differ by 0.3, there is no finding in either direction. `noise_sigma` measures it, and every
threshold is derived from it rather than chosen by hand. Measure σ before curating the full pair
set — it sets the replicate count for every later run.

Two sources of run-to-run variance, and the seed covers both: the diffusion noise, and the
**random choice of input conformer** (`featurizerv2.py:1193-1198`, unseeded `random.choice`)
plus the **random roto-translation applied to it** (`featurizerv2.py:1494-1499`). Measure σ by
varying `--seed` over a fixed input rather than by rerunning unseeded — same quantity, but
reproducible and attributable.

## Pair-set design

Two groups, not one. Testing only that Boltz-2 *responds* to stereochemistry would be passed
by a model that emits a different number for any perturbed input:

- **tail-enriched pairs** — large experimental Δ, where a flat prediction is unambiguous
- **near-zero pairs** — enantiomers with small measured Δ, where the model *should* be flat

Pairs must come from the **same assay**, not merely the same target and activity type: cross-lab
IC50 variance can exceed the stereochemical effect. Censored values (`>`, `<`) are excluded —
they are the likeliest source of the implausible fold-differences reported in archive-scale
chirality datasets.

Candidate upstream sources are recorded in [`data/external/SOURCES.md`](data/external/SOURCES.md).
Nothing is vendored yet.

## Install

```bash
uv venv --python 3.12 .venv && uv pip install --python .venv/bin/python -e '.[dev]'
```
```bash
.venv/bin/python -m pytest -q
```

Boltz itself is a separate install and is not a dependency of this package — the harness shells
out to the `boltz predict` CLI. See [jwohlwend/boltz](https://github.com/jwohlwend/boltz).

## Boltz-2 interface notes

Verified against the source at `v2.2.1` (`b1ebfc4`), not just the docs:

- Steering potentials are **`--use_potentials`, a flag that defaults to OFF**. There is no
  `--no_potentials`; the unsteered arm is the default arm. **But it does not reach the affinity
  stage**: that stage hard-disables `fk_steering`, `physical_guidance_update` and
  `contact_guidance_update` regardless of the flag (`main.py:1387-1389`). Steering can only
  change the structure handed to the affinity stage, never the geometry the affinity head scores.
- Affinity emits **one `affinity_<name>.json` per input**, not per diffusion sample.
  `affinity_pred_value` is log10(IC50) with IC50 in µM, so **lower is stronger**.
- `--diffusion_samples_affinity` (**default 3**, `main.py:1054`) is separate from
  `--diffusion_samples`, and `--sampling_steps_affinity` defaults to 200 (`main.py:1053`).
- Affinity takes **exactly one ligand chain**, ≤128 atoms, with ≤56 recommended — enforced in
  `parse/schema.py:1266-1272`, matching `boltz_io`'s limits.
- `affinity_pred_value` carries a **molecular-weight correction** applied after the ensemble
  average (`boltz2.py:688-697`). MW is identical between enantiomers, so it cancels in Δ — but
  it means the reported number is not the raw head output.
- Both members of a pair must be given the **same MSA**. The protein is identical between
  enantiomers; a differing MSA would inject a difference unrelated to stereochemistry.

## Source findings

Read at `jwohlwend/boltz` `v2.2.1` (`b1ebfc4`), 2026-11-25. Line numbers are that commit.
Both original open questions are answered, and a third channel turned up that the original
framing did not account for.

### 1. No chiral atom feature reaches the model — but chirality still arrives without the pose

The featuriser **does** compute a chiral tag per atom: `ref_chirality`, built from
`GetChiralTag()` over the reference mol (`featurizerv2.py:1227-1231, 1436`) and placed in the
feature dict (`featurizerv2.py:1549`). It is then **never read**. The only consumers of atom
reference features are the atom encoders, which take `ref_pos`, `ref_charge`, `ref_element`,
`ref_atom_name_chars` and `ref_space_uid` (`encodersv2.py:313-325`, `encoders.py:404-425`);
`ref_chirality` appears nowhere under `src/boltz/model/`. So there is no 2D chiral-tag channel.

**That is not the end of it.** `ref_pos` is the **RDKit ETKDGv3 reference conformer**, embedded
from the input SMILES with `AddHs` (`parse/schema.py:1242-1258`, `compute_3d_conformer` at
`:200-250`). ETKDG honours the chiral tags, so **two enantiomers get mirror-image reference
conformers**. The atom encoder embeds pairwise displacement vectors between reference-conformer
atoms through a plain `LinearNoBias(3, atom_z)` (`encodersv2.py:263, 354-381`), which is not
reflection-invariant. The augmentation applied beforehand is a **proper** rotation —
`randomly_rotate` → `random_rotations` (`modules/utils.py:68-76, 99-104`), det +1 — so it
scrambles orientation but preserves handedness.

Consequence for the experiment: stereochemistry can reach the affinity head **through the input
conformer**, on a path that never passes through the predicted pose. A non-zero Δ under the key
conditional therefore does not by itself prove the head is reading its scored coordinates — the
conformer is a live alternative explanation. Distinguishing them needs a third arm: hold the
input SMILES fixed and substitute the reference conformer, or supply both members with the same
conformer. This is a new requirement, not in the original design.

### 2. Confirmed — the affinity head scores its own, unwritten coordinates

The affinity stage is a **second, separate `trainer.predict` run** against a **different
checkpoint**, `boltz2_aff.ckpt` (`main.py:1370-1410`). It takes the rank-0 structure from the
first run as input — written as `pre_affinity_<id>.npz` (`write/writer.py:178-180`), loaded at
`module/inferencev2.py:62-64` — and then **runs its own trunk and its own diffusion** over it
(`recycling_steps=5`, `sampling_steps_affinity`, `diffusion_samples_affinity`;
`main.py:1372-1379`).

The head is then fed the best sample from *that* pass, ranked by ipTM:

```python
argsort = torch.argsort(dict_out["iptm"], descending=True)
best_idx = argsort[0].item()
coords_affinity = dict_out["sample_atom_coords"].detach()[best_idx][None, None]
```
`model/models/boltz2.py:622-626`

`BoltzAffinityWriter` writes **only the JSON summary** — no structure at all
(`write/writer.py:292-330`). So the stereochemistry in the output CIF is from the *structure*
model and is **not** what the affinity head saw. Checkpoint 3 is unreachable with a stock
install, exactly as suspected.

**Capture point for the patch:** `coords_affinity` at `boltz2.py:624`, or equivalently
`dict_out["sample_atom_coords"][best_idx]`. Atom identity and ordering come from the same
`feats`, so a dump there plus the existing structure npz is enough to rebuild the ligand.

### 3. Determinism

`--seed` exists (`main.py:919`) and defaults to None (`main.py:1101-1103`). Unseeded, two runs
of one input differ. See *The noise floor is not optional* for what that implies and for the two
distinct variance sources.

### Verification checkpoints

| Checkpoint | Reachable? |
|---|---|
| 0. Input conformer — what geometry did ETKDG hand the trunk? | yes, but only by re-embedding; not written out. **Newly relevant — see finding 1** |
| 1. Parsed input — did the requested R/S enter correctly? | yes, RDKit on the input SMILES |
| 2. Predicted complex — did structure generation preserve it? | yes, `AssignStereochemistryFrom3D` on the output CIF |
| 3. Affinity-stage complex — did the affinity pass preserve it? | **not without a patch** — confirmed never written; capture at `boltz2.py:624` |
| 4. Affinity output — right direction and magnitude? | yes, `affinity_<name>.json` |

## Prior work

- Goddard & Croll, UCSF ChimeraX analysis (2026) — flipped stereocenters in 27 of 52 HIV-protease
  ligands; 3 of 50 with steering. Establishes the *generation* error, not an affinity failure.
- Ujihara, *Pose-to-Affinity Information Flow in Boltz-2* (2026 preprint) — substituting
  crystallographic poses moved affinity <3-fold.
- Fried / Rowan (Nov 2026) — rescoring alternative poses with fixed upstream representations gave
  median affinity ranges of 1.16-fold with MSAs, 1.20-fold without. Explicitly did not evaluate
  experimental affinity accuracy.
- Schneider, Lewis, Fechner & Ertl, *Chiral Cliffs*, ChemMedChem 2018;13(13):1315–1324 — the
  origin of the question at archive scale.

This repo is the affinity-side analogue of what PoseBusters does for poses: if Boltz-2 nails
stereo, that is itself a notable result given the mutation-insensitivity literature.

## Licence

MIT (code). Any vendored dataset keeps its own licence — see `data/external/SOURCES.md`.
