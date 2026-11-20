# boltz2-stereoselectivity

**Does Boltz-2's affinity head respond to stereochemistry?**

Enantiomers are the cleanest possible probe of a structure-based affinity model.
Identical atoms, identical bonds, identical formula, identical 2D graph — mirror-image
geometry. Experimentally, one enantiomer frequently binds far more tightly than the other.
If Boltz-2 assigns them the same affinity, the head is reading chemistry, not shape.

There is a second, already-documented failure to attach this to: Boltz-2 often generates the
wrong enantiomer from a stereo-specified SMILES. That allows the paired question — *when the
geometry comes out mirrored, does the predicted affinity follow the geometry or the input
string?*

Existing probes of Boltz-2's coordinate sensitivity have all worked the protein side
(mutate the pocket, shuffle targets) or substituted whole poses. Nobody has held the 2D
graph fixed and changed only the 3D.

## Status

Early. The verification and evaluation core is implemented and tested; the curated pair set
is not built yet, and no inference has been run. See *Open questions* — two of them bear on
whether the experiment is interpretable at all, and both are read-only questions about the
Boltz source.

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

The headline result is a null, and a null is unreadable without a noise floor. Boltz-2
exposes **no seed flag**, so repeated runs of one input differ by some σ. If σ is 0.4 log
units and the pairs differ by 0.3, there is no finding in either direction. `noise_sigma`
measures it from replicate runs of identical inputs, and every threshold is derived from it
rather than chosen by hand. Measure σ before curating the full pair set — it sets the
replicate count for every later run.

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

From the Boltz prediction docs, since these shape the harness:

- Steering potentials are **`--use_potentials`, a flag that defaults to OFF**. There is no
  `--no_potentials`; the unsteered arm is the default arm.
- Affinity emits **one `affinity_<name>.json` per input**, not per diffusion sample.
  `affinity_pred_value` is log10(IC50) with IC50 in µM, so **lower is stronger**.
- `--diffusion_samples_affinity` (default 5) is separate from `--diffusion_samples` (default 1).
- Affinity takes **exactly one ligand chain**, ≤128 atoms, with ≤56 recommended.
- Both members of a pair must be given the **same MSA**. The protein is identical between
  enantiomers; a differing MSA would inject a difference unrelated to stereochemistry.

## Open questions

These are unresolved and both are answerable by reading the Boltz source. Neither has been
checked yet.

**1. Can the affinity head read stereochemistry without geometry?** Enantiomers have different
isomeric SMILES. If Boltz-2's ligand featurization carries a chiral tag or CIP code as an *atom
feature*, stereochemistry reaches the affinity head through a 2D channel that never passes
through coordinates — and a positive result would not demonstrate 3D sensitivity. If there is
no such feature, stereo information can only arrive via coordinates, which makes the
experiment clean.

**2. Does the affinity module use the structure you inspected?** The separate
`--diffusion_samples_affinity` suggests the affinity module runs its own diffusion pass rather
than consuming the ranked structure from the main run. If so, the stereochemistry in the output
CIF is *not* the geometry the affinity head saw, and the paired question needs the
affinity-stage coordinates instead — which the CLI does not write out.

### Verification checkpoints

| Checkpoint | Reachable? |
|---|---|
| 1. Parsed input — did the requested R/S enter correctly? | yes, RDKit on the input SMILES |
| 2. Predicted complex — did structure generation preserve it? | yes, `AssignStereochemistryFrom3D` on the output CIF |
| 3. Affinity-stage complex — did the affinity pass preserve it? | **no** — not written to disk; see open question 2 |
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
