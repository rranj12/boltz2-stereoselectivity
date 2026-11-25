"""Reference-conformer construction, reflection, and substitution.

Why this module exists
----------------------
Boltz-2 does not only see the predicted pose. Before any network runs, it
embeds a 3D reference conformer from the input SMILES with RDKit ETKDGv3 plus a
UFF cleanup (`boltz/data/parse/schema.py:200-254`, called at 1258). That
geometry is fed to the atom encoder as `ref_pos`, where the encoder forms
pairwise displacement vectors between reference atoms

    d = ref_pos[j] - ref_pos[i]

and embeds them through `LinearNoBias(3, atom_z)`
(`boltz/model/modules/encodersv2.py:263`, applied at 357 and 381). A linear map
of a displacement vector is *not* reflection-invariant, so the handedness of
the reference conformer is information the network can read -- on a path that
never passes through the predicted pose.

The only transform applied beforehand is `center_random_augmentation`
(`boltz/model/modules/utils.py:68-104`), a proper rotation drawn from SO(3)
plus a translation. Proper rotations preserve handedness, so the chirality of
the reference conformer survives into the features.

Running alongside it is a second channel, `embed_atompair_ref_dist`, which sees
only ``1 / (1 + |d|^2)``. Distances are reflection-invariant, so that channel is
blind to handedness by construction.

That asymmetry is what makes reflection the right intervention. Mirroring a
conformer preserves every interatomic distance exactly, so the distance channel
is bit-identical and *only* the displacement channel changes. Any response to a
reflected conformer localises to one linear layer.

What reflection is and is not
-----------------------------
`reflect_coords` negates a single Cartesian axis. That is an improper transform
(determinant -1), so it maps a conformer onto its mirror image. The result is a
chemically real geometry -- it is a valid conformer of the opposite enantiomer.

What is synthetic is the *pairing*: `reflect_conformer_only` leaves the
molecular graph untouched while mirroring the coordinates, producing a molecule
whose SMILES says one thing and whose reference geometry says the other. That
contradiction is deliberate. It is the crossover cell that reveals which input
the affinity head is reading.

Note that the graph's own chiral tags are *not* what reaches the model. Boltz
computes a per-atom `ref_chirality` feature (`featurizerv2.py:1227`) and ships
it in the feature dict (line 1549), but no module under `boltz/model/` ever
reads it. It is a dead feature. Stereochemistry reaches the network through
coordinates only -- the reference conformer here, and the predicted pose.
"""

from __future__ import annotations

import pickle
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from stereobench.stereo import Labels, labels_from_3d, labels_from_mol

# Boltz tags its generated conformer with this name and looks it up by it
# (schema.py:248 writes it, get_conformer at 278 reads it). A pickled mol
# missing the tag still loads via Boltz's positional fallback, but we set it so
# substitution follows the primary path rather than the fallback.
CONFORMER_NAME = "Computed"

# Mirroring negates this Cartesian axis. The choice is arbitrary -- any single
# axis gives an improper transform, and the result differs only by a rotation.
DEFAULT_REFLECT_AXIS = 0


def build_reference_conformer(
    smiles: str,
    *,
    version: str = "v3",
    seed: int | None = None,
    optimize: bool = True,
) -> Chem.Mol:
    """Build the reference conformer the way Boltz-2 builds it.

    Mirrors `boltz.data.parse.schema.compute_3d_conformer`, including the
    explicit hydrogens, the UFF cleanup, and the conformer name tag. Atom names
    follow Boltz's canonical-rank scheme (schema.py:1246-1257) so the resulting
    mol can be substituted into a Boltz run.

    `seed` is exposed because Boltz does not set one. ETKDG is stochastic, so
    two runs of one SMILES give slightly different reference geometries. Pin it
    when a pair must share an embedding run, and leave it None to sample the
    variation Boltz itself is subject to.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"unparseable SMILES: {smiles!r}")
    mol = Chem.AddHs(mol)

    # Boltz's atom naming, reproduced so substituted mols match its expectations.
    canonical_order = Chem.CanonicalRankAtoms(mol)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    for atom, can_idx in zip(mol.GetAtoms(), canonical_order):
        name = atom.GetSymbol().upper() + str(can_idx + 1)
        if len(name) > 4:
            raise ValueError(f"atom name longer than 4 characters: {name}")
        atom.SetProp("name", name)

    options = AllChem.ETKDGv3() if version == "v3" else AllChem.ETKDGv2()
    options.clearConfs = False
    if seed is not None:
        options.randomSeed = int(seed)

    conf_id = AllChem.EmbedMolecule(mol, options)
    if conf_id == -1:
        options.useRandomCoords = True
        conf_id = AllChem.EmbedMolecule(mol, options)
    if conf_id == -1:
        raise RuntimeError(f"ETKDG failed to embed a conformer for {smiles!r}")

    if optimize:
        try:
            AllChem.UFFOptimizeMolecule(mol, confId=conf_id, maxIters=1000)
        except (RuntimeError, ValueError):
            # Boltz swallows both of these and keeps the unoptimised conformer.
            pass

    conformer = mol.GetConformer(conf_id)
    conformer.SetProp("name", CONFORMER_NAME)
    conformer.SetProp("coord_generation", f"ETKDG{version}")

    # The featuriser picks a conformer with `random.choice` over everything the
    # mol carries (`featurizerv2.py:1193-1198`). A substituted mol holding more
    # than one would reintroduce exactly the run-to-run variance a controlled
    # arm exists to remove, so keep it to one and fail loudly if that changes.
    if mol.GetNumConformers() != 1:
        raise RuntimeError(
            f"expected exactly one conformer, got {mol.GetNumConformers()}; "
            "Boltz samples conformers at random and a multi-conformer mol "
            "would make the input nondeterministic"
        )
    return mol


def get_coords(mol: Chem.Mol, conf_id: int = -1) -> np.ndarray:
    """Conformer coordinates as an (n_atoms, 3) array."""
    if mol.GetNumConformers() == 0:
        raise ValueError("mol has no conformer")
    return np.array(mol.GetConformer(conf_id).GetPositions(), dtype=float)


def reflect_coords(coords: np.ndarray, axis: int = DEFAULT_REFLECT_AXIS) -> np.ndarray:
    """Mirror coordinates through a plane by negating one axis.

    Determinant -1, so handedness flips. Every pairwise distance is preserved
    exactly, which is the property that makes this a clean probe: Boltz's
    distance channel cannot see the change, and only its displacement channel
    can.
    """
    if axis not in (0, 1, 2):
        raise ValueError(f"axis must be 0, 1 or 2; got {axis}")
    out = np.array(coords, dtype=float, copy=True)
    out[:, axis] *= -1.0
    return out


def _set_coords(mol: Chem.Mol, coords: np.ndarray, conf_id: int = -1) -> None:
    conformer = mol.GetConformer(conf_id)
    for i in range(mol.GetNumAtoms()):
        x, y, z = (float(v) for v in coords[i])
        conformer.SetAtomPosition(i, Chem.rdGeometry.Point3D(x, y, z))


def reflect_conformer_only(
    mol: Chem.Mol, *, axis: int = DEFAULT_REFLECT_AXIS
) -> Chem.Mol:
    """Mirror the geometry and leave the molecular graph alone.

    The returned mol is deliberately self-contradictory: its graph still
    carries the original chiral tags, while its conformer now has the opposite
    handedness. This is the crossover input -- it changes exactly one of the two
    channels that can carry stereochemistry into Boltz.

    Use `reflect_molecule` instead if you want the genuine enantiomer.
    """
    out = Chem.Mol(mol)
    _set_coords(out, reflect_coords(get_coords(mol), axis=axis))
    return out


def reflect_molecule(mol: Chem.Mol, *, axis: int = DEFAULT_REFLECT_AXIS) -> Chem.Mol:
    """Mirror geometry *and* graph -- the honest enantiomer, no contradiction.

    Provided as the control for `reflect_conformer_only`: both produce the same
    coordinates, and they differ only in whether the graph was updated to
    match. Comparing them separates a response to the conformer from a response
    to the graph.
    """
    out = reflect_conformer_only(mol, axis=axis)
    Chem.AssignStereochemistryFrom3D(out)
    return out


def conformer_labels(mol: Chem.Mol) -> Labels:
    """CIP labels read off the geometry, keyed by graph-canonical rank."""
    return labels_from_3d(mol)


def graph_labels(mol: Chem.Mol) -> Labels:
    """CIP labels carried by the molecular graph, ignoring any conformer.

    Read off the graph in place rather than via a SMILES round trip. Canonical
    ranks are computed over whatever graph they are given, and re-parsing a
    SMILES renumbers the atoms, so a round trip would return labels keyed into
    a different space than `conformer_labels` uses and the two could not be
    compared.
    """
    stripped = Chem.Mol(mol)
    stripped.RemoveAllConformers()
    return labels_from_mol(stripped)


def max_distance_deviation(a: Chem.Mol, b: Chem.Mol) -> float:
    """Largest change in any interatomic distance between two conformers.

    The validity check for a reflection. A true mirror leaves every distance
    untouched, so this must come back at numerical noise. Anything larger means
    the geometry moved in some way beyond reflection, and the intervention is
    no longer clean -- a response could then be to shape rather than handedness.
    """
    ca, cb = get_coords(a), get_coords(b)
    if ca.shape != cb.shape:
        raise ValueError("conformers have different atom counts")
    da = np.linalg.norm(ca[:, None, :] - ca[None, :, :], axis=-1)
    db = np.linalg.norm(cb[:, None, :] - cb[None, :, :], axis=-1)
    return float(np.max(np.abs(da - db)))


@dataclass(frozen=True)
class CrossoverArm:
    """One cell of the SMILES x conformer crossover.

    `graph_source` / `conformer_source` name which member of the pair supplied
    each channel. When they differ, the input is a deliberate contradiction and
    the affinity the model returns says which channel it followed.
    """

    arm_id: str
    smiles: str
    graph_source: str
    conformer_source: str
    mol: Chem.Mol

    @property
    def is_crossed(self) -> bool:
        return self.graph_source != self.conformer_source

    def summary(self) -> dict:
        return {
            "arm_id": self.arm_id,
            "graph_source": self.graph_source,
            "conformer_source": self.conformer_source,
            "is_crossed": self.is_crossed,
            "graph_labels": graph_labels(self.mol),
            "conformer_labels": conformer_labels(self.mol),
        }


def build_crossover(
    smiles_a: str,
    smiles_b: str,
    *,
    seed: int | None = None,
    axis: int = DEFAULT_REFLECT_AXIS,
) -> list[CrossoverArm]:
    """The 2x2: each member's graph against each member's conformer handedness.

    Every arm descends from a SINGLE embedding of A. B's conformer is A's
    mirrored, not an independent ETKDG run.

    That is not a shortcut, it is the thing that makes the experiment valid.
    Embedding the two members separately -- even under one seed -- gives
    genuinely different geometries: ETKDG lands them in different torsional
    minima, and measured deviations run to well over 1 angstrom. The arms would
    then differ in shape as well as handedness, and any affinity response could
    be blamed on either. Deriving by reflection holds all 3N-6 internal
    coordinates fixed and varies only the sign, so the four arms differ in
    exactly the two bits under test.

    The result is two geometries crossed against two graphs. `a_crossed` and
    `b_native` share coordinates exactly, as do `b_crossed` and `a_native`;
    within each of those couples only the molecular graph differs.

    Raises if the two SMILES are not an enantiomeric pair, since reflecting A
    would not then produce B and the crossing would be meaningless.
    """
    mol_a = build_reference_conformer(smiles_a, seed=seed)
    mol_b = reflect_molecule(mol_a, axis=axis)

    produced = Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(mol_b)))
    expected = Chem.MolToSmiles(Chem.MolFromSmiles(smiles_b))
    if produced != expected:
        raise ValueError(
            "mirroring A did not reproduce B, so these are not enantiomers: "
            f"got {produced!r}, expected {expected!r}"
        )

    return [
        CrossoverArm("a_native", smiles_a, "a", "a", mol_a),
        CrossoverArm(
            "a_crossed", smiles_a, "a", "b", reflect_conformer_only(mol_a, axis=axis)
        ),
        CrossoverArm("b_native", smiles_b, "b", "b", mol_b),
        CrossoverArm(
            "b_crossed", smiles_b, "b", "a", reflect_conformer_only(mol_b, axis=axis)
        ),
    ]


def write_mol_pickle(mol: Chem.Mol, path: str | Path) -> Path:
    """Write a mol in the form Boltz's molecule cache expects.

    `boltz.data.mol.load_molecules` reads `<name>.pkl` from a directory and
    unpickles an RDKit Mol (`boltz/data/mol.py:16-36`). Dropping a file there
    and referencing its stem as a `ccd` code in the YAML substitutes this
    geometry for the one ETKDG would have produced -- which is how a crossover
    arm reaches the model without patching Boltz.

    RDKit drops atom and conformer properties when pickling unless the process
    is told otherwise, and Boltz needs those: `parse_ccd_residue` reads each
    atom's `name`, and `get_conformer` looks up the conformer's. Boltz sets the
    same global flag itself (`boltz/main.py:645` and `:1099`), so this matches
    its behaviour rather than imposing anything new -- but it IS process-global
    state, and it changes how every later RDKit pickle in this process behaves.

    Verify the route with `scripts/verify_substitution.py` before relying on it.
    """
    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        pickle.dump(mol, handle)
    return path
