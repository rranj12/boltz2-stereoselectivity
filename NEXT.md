# Where this stands and what to do next

Written 2026-11-25 overnight. Everything below runs without a GPU except where
it says otherwise.

## The one-line state

The mechanistic arm is built, tested, and **de-risked** — conformer
substitution reaches Boltz with no fork required, verified against a real
install. The pair set is in hand: 261 pairs, all 261 of which now build a
crossover. Nothing has been run through Boltz itself.

The Boltz install used for verification lives in this session's scratchpad
under `/private/tmp`, which is cleared on reboot. Reinstall it somewhere
durable before relying on it.

## What got answered

**Conformer substitution works without patching Boltz.** This was the biggest
open risk and it is closed. `scripts/verify_substitution.py` checks it against
a real Boltz install (`v2.2.1`), end to end: our pickled mol round-trips
through Boltz's own loader with atom names and conformer intact, the
coordinates that survive are ours rather than a re-embedded ETKDG conformer,
`parse_ccd_residue` accepts it, and native/reflected arms come back with
opposite handedness. All checks pass.

That means the crossover arm needs no fork. A fork is still required for
checkpoint 3 (the affinity-stage coordinates), but that is now a separate,
optional question rather than a blocker.

**The analysis can tell the hypotheses apart.** `scripts/analyze.py
--synthetic {conformer,smiles,flat}` fabricates outputs under a stated
hypothesis and runs the real analysis path over them. It correctly identifies
all three. This is a power check, not a demo: an analysis that could not
separate them on data where the answer is known would not separate them on
real data either.

## Run order when you have a GPU

```bash
.venv/bin/python -m pytest -q
```

```bash
PYTHONPATH=src /path/to/boltz-venv/bin/python scripts/verify_substitution.py
```

```bash
PYTHONPATH=src /path/to/boltz-venv/bin/python scripts/preflight.py
```

Then, smallest useful experiment first:

**1. Noise floor, before anything else.** It is the go/no-go gate. If σ is
0.5 pAct and pairs differ by 0.3, there is no finding in either direction and
no amount of curation rescues it. An hour of GPU.

```bash
.venv/bin/python scripts/build_noise_floor.py --msa-dir data/msa
cp runs/noise_floor/mols/*.pkl ~/.boltz/mols/
runs/noise_floor/run.sh
```

Two arms, deliberately. `substituted` pins the reference conformer, so only
Boltz's own stochasticity varies — that is the floor the crossover is measured
against. `native` lets Boltz embed its own conformer, which is additionally
unseeded. The gap between the two σ values *is* the contribution of the
unseeded conformer, which is the quantity the mechanistic arm is about.

**2. The crossover, on 3–5 pairs.** The headline experiment. Needs no
experimental data, so it does not wait on curation.

```bash
PYTHONPATH=src /path/to/boltz-venv/bin/python scripts/fetch_msas.py --max-targets 8
```

```bash
.venv/bin/python scripts/build_run.py --mode crossover --max-pairs 5 --replicates 3 --msa-dir data/msa
cp runs/pilot/mols/*.pkl ~/.boltz/mols/
runs/pilot/run.sh
```

`run.sh` is generated from the manifest, one `boltz predict` per input carrying
the seed the manifest records. Do not run the YAMLs by hand: `--seed` is a CLI
flag and not a YAML field, so by-hand runs leave the manifest's seeds fictional
and the replicates either unseeded or mislabelled. It skips inputs whose output
directory exists, so it resumes.

Then:

```bash
.venv/bin/python scripts/analyze.py
```

**3. The chemistry arm**, once `data/pairs/pilot.csv` is in hand.

## What each crossover outcome means

| Outcome | Reading |
|---|---|
| affinity follows the **conformer** | input channel is live; every paired comparison must hold conformers fixed |
| affinity follows the **SMILES** | something downstream of `ref_pos` carries stereochemistry |
| affinity follows **neither** | stereochemistry reaches the head only via the pose — the original design is clean as written |

The third is the one to hope for. The test's best case *restores* the original
experiment rather than complicating it.

## Known limitations to state in the writeup

These are worth naming explicitly rather than engineering around. For a class
project, naming an identification problem precisely is worth more than
half-solving it.

- **The crossover perturbs two things at once.** Changing `ref_pos` also
  changes the diffusion conditioning, so the pose moves too. Isolating the
  conformer cleanly means holding `x_pred` fixed, which needs the fork
  (capture point: `boltz2.py:624`).
- **Checkpoint 3 is unreachable.** The affinity head scores a re-diffused
  structure that is never written to disk. Any stereo verdict from an output
  CIF is checkpoint 2.
- **`pilot.csv` pairs are same-assay but not same-lab-protocol verified.**
  One assay id is a strong constraint, not a perfect one.
- **The MSA must be shared across a pair's arms.** With `--use_msa_server`,
  Boltz generates one per input file, so the four arms of a pair can be
  conditioned on different alignments — variance injected into exactly the
  comparison being measured, and invisible in the output. `fetch_msas.py`
  generates one per target and `build_run.py --msa-dir` points every arm at it.
  `build_run.py` warns when run without it.
- **ETKDG is unseeded inside Boltz.** `compute_3d_conformer` never sets
  `options.randomSeed`, so Boltz's own reference conformers vary run to run
  independently of `--seed`. Substituted mols remove this; native arms do not.

## Pipeline map

| Script | Does | Needs |
|---|---|---|
| `scripts/shortlist.py` | StereoAtlas → candidate pairs | the two CSVs (downloaded) |
| `scripts/requalify.py` | ChEMBL → SMILES, assay, qualifier → `pilot.csv` | network |
| `scripts/fetch_sequences.py` | ChEMBL target → UniProt sequence | network |
| `scripts/fetch_msas.py` | one MSA per target, cached | network + Boltz |
| `scripts/build_run.py` | pairs → YAMLs, substituted mols, manifest, `run.sh` | — |
| `scripts/build_noise_floor.py` | replicate inputs for σ, both arms | — |
| `scripts/preflight.py` | catches what fails 20 min into a GPU session | Boltz for the deep checks |
| `scripts/verify_substitution.py` | proves substitution reaches Boltz | Boltz installed |
| `scripts/analyze.py` | outputs → pilot table + metrics | a run, or `--synthetic` |

## New modules

| Module | Does |
|---|---|
| `stereobench.conformer` | builds Boltz's ETKDG reference conformer, reflects it, writes substitutable mols, assembles the 2×2 crossover |
| `stereobench.structure` | reads the predicted ligand out of a Boltz CIF by atom-name mapping onto a template; the missing link between a prediction and a stereo verdict |
| `stereobench.curate` | StereoAtlas joins, plausibility and enantiomer-candidate filters |
| `stereobench.chembl` | SMILES, same-assay lookup, qualifiers, target sequences |

## Two traps worth remembering

**Label rank spaces.** CIP labels are keyed by graph-canonical rank, and ranks
are computed over whatever graph they are given. Labels from a
hydrogen-bearing mol cannot be compared against labels from a stripped one —
they agree on the verdict and disagree on the key, which reads as "the centres
do not correspond". Use `structure.template_labels` /
`template_conformer_labels` to get both sides into one space. Two tests
already caught this; it will try again.

**Over-specifying stereochemistry.** `AssignStereochemistryFrom3D` tags every
atom whose neighbours are distinguishable in 3D, including atoms no
stereocentre supports — adamantyl bridgeheads are the common case. The mol is
still correct, but it no longer prints as the SMILES it came from, so an
identity check against the input string fails. This cost 10% of the pair set
before it was found, silently, at run-build time. `reflect_molecule` now
inverts the existing chiral tags instead, which keeps the specification exactly
as specific as the input was.

**One conformer per substituted mol.** The featuriser picks one at
`random.choice` (`featurizerv2.py:1193-1198`), so a multi-conformer mol
reintroduces the variance the controlled arm exists to remove.
`build_reference_conformer` raises if it ever produces more than one.
