"""Tests for reference-conformer construction and reflection.

The load-bearing test in this file is `test_reflection_preserves_all_distances`.
The whole reflection design rests on the claim that mirroring changes handedness
and nothing else -- if distances moved, Boltz's distance channel would see the
change too and a response could no longer be attributed to chirality.
"""

from __future__ import annotations

import pickle

import pytest
from rdkit import Chem

from stereobench.conformer import (
    CONFORMER_NAME,
    build_crossover,
    build_reference_conformer,
    conformer_labels,
    get_coords,
    graph_labels,
    max_distance_deviation,
    reflect_coords,
    reflect_conformer_only,
    reflect_molecule,
    write_mol_pickle,
)
from stereobench.stereo import classify_pair

# One stereocentre, small, embeds reliably.
R_PHENYLETHANOL = "C[C@H](O)c1ccccc1"
S_PHENYLETHANOL = "C[C@@H](O)c1ccccc1"
ACHIRAL = "CC(C)O"
# Meso tartaric acid: two stereocentres, but its own mirror image.
MESO_TARTARIC = "O[C@H](C(=O)O)[C@@H](O)C(=O)O"

SEED = 0xC0FFEE


@pytest.fixture(scope="module")
def mol_r():
    return build_reference_conformer(R_PHENYLETHANOL, seed=SEED)


def test_every_arm_carries_exactly_one_conformer():
    """Boltz picks a conformer at random (`featurizerv2.py:1193-1198`).

    A substituted mol holding several would make `ref_pos` nondeterministic,
    reintroducing the variance the crossover arm exists to eliminate.
    """
    for arm in build_crossover(R_PHENYLETHANOL, S_PHENYLETHANOL, seed=SEED):
        assert arm.mol.GetNumConformers() == 1, arm.arm_id


def test_builds_conformer_with_boltz_tags(mol_r):
    assert mol_r.GetNumConformers() == 1
    conf = mol_r.GetConformer()
    assert conf.GetProp("name") == CONFORMER_NAME
    # Boltz names every atom; a mol missing these cannot be substituted in.
    assert all(a.HasProp("name") for a in mol_r.GetAtoms())


def test_etkdg_honours_requested_chirality(mol_r):
    """The embedded geometry must match the SMILES it came from.

    If this failed, the reference conformer would carry no reliable
    stereochemical signal and the whole crossover design would be moot.
    """
    assert conformer_labels(mol_r) == graph_labels(mol_r)


def test_reflection_preserves_all_distances(mol_r):
    """A mirror leaves every interatomic distance exactly where it was.

    This is what makes reflection a clean probe rather than just a
    perturbation: Boltz's `embed_atompair_ref_dist` channel sees only
    distances, so it is provably blind to this intervention, and any response
    must come through the displacement channel.
    """
    mirrored = reflect_conformer_only(mol_r)
    assert max_distance_deviation(mol_r, mirrored) < 1e-9


def test_reflection_flips_geometric_handedness(mol_r):
    mirrored = reflect_conformer_only(mol_r)
    before, after = conformer_labels(mol_r), conformer_labels(mirrored)
    assert set(before) == set(after)
    assert all(before[k] != after[k] for k in before)


def test_reflect_conformer_only_leaves_the_graph_alone(mol_r):
    """The crossover input: graph says one thing, geometry says the other."""
    mirrored = reflect_conformer_only(mol_r)
    assert graph_labels(mirrored) == graph_labels(mol_r)
    assert conformer_labels(mirrored) != conformer_labels(mol_r)
    # The contradiction is the point.
    assert graph_labels(mirrored) != conformer_labels(mirrored)


def test_reflect_molecule_updates_both(mol_r):
    """The honest enantiomer: no contradiction between graph and geometry."""
    mirrored = reflect_molecule(mol_r)
    assert graph_labels(mirrored) == conformer_labels(mirrored)
    assert graph_labels(mirrored) != graph_labels(mol_r)


def test_reflect_molecule_gives_the_enantiomer(mol_r):
    mirrored = reflect_molecule(mol_r)
    a = Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(mol_r)))
    b = Chem.MolToSmiles(Chem.RemoveHs(mirrored))
    assert classify_pair(a, b) == "enantiomers"


def test_reflection_is_an_improper_transform():
    import numpy as np

    coords = np.array([[1.0, 2.0, 3.0], [-1.0, 0.5, 2.0]])
    for axis in (0, 1, 2):
        out = reflect_coords(coords, axis=axis)
        assert np.allclose(out[:, axis], -coords[:, axis])
        others = [i for i in (0, 1, 2) if i != axis]
        assert np.allclose(out[:, others], coords[:, others])


def test_reflect_coords_rejects_bad_axis():
    import numpy as np

    with pytest.raises(ValueError, match="axis must be"):
        reflect_coords(np.zeros((2, 3)), axis=3)


def test_achiral_molecule_has_no_labels_to_flip():
    mol = build_reference_conformer(ACHIRAL, seed=SEED)
    assert conformer_labels(mol) == {}
    assert conformer_labels(reflect_conformer_only(mol)) == {}


def test_meso_compound_is_its_own_mirror_image():
    """Reflecting a meso compound gives back the same molecule.

    The counterpart to `classify_pair`'s meso test: counting inverted centres
    would call this an enantiomeric pair, and it is not one.
    """
    mol = build_reference_conformer(MESO_TARTARIC, seed=SEED)
    mirrored = reflect_molecule(mol)
    a = Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(mol)))
    b = Chem.MolToSmiles(Chem.RemoveHs(mirrored))
    assert classify_pair(a, b) == "identical"


def test_crossover_has_four_arms_two_crossed():
    arms = build_crossover(R_PHENYLETHANOL, S_PHENYLETHANOL, seed=SEED)
    assert [a.arm_id for a in arms] == [
        "a_native",
        "a_crossed",
        "b_native",
        "b_crossed",
    ]
    assert [a.is_crossed for a in arms] == [False, True, False, True]


def test_crossover_native_arms_are_self_consistent():
    arms = {a.arm_id: a for a in build_crossover(R_PHENYLETHANOL, S_PHENYLETHANOL, seed=SEED)}
    for arm_id in ("a_native", "b_native"):
        arm = arms[arm_id]
        assert graph_labels(arm.mol) == conformer_labels(arm.mol)


def test_crossover_crossed_arms_contradict_themselves():
    arms = {a.arm_id: a for a in build_crossover(R_PHENYLETHANOL, S_PHENYLETHANOL, seed=SEED)}
    for arm_id in ("a_crossed", "b_crossed"):
        arm = arms[arm_id]
        assert graph_labels(arm.mol) != conformer_labels(arm.mol)


def test_crossed_arm_geometry_matches_the_other_member():
    """a_crossed must carry B's handedness while keeping A's graph.

    If this did not hold, the crossed arm would not be testing what the design
    claims it tests.
    """
    arms = {a.arm_id: a for a in build_crossover(R_PHENYLETHANOL, S_PHENYLETHANOL, seed=SEED)}
    assert conformer_labels(arms["a_crossed"].mol) == conformer_labels(arms["b_native"].mol)
    assert graph_labels(arms["a_crossed"].mol) == graph_labels(arms["a_native"].mol)


def test_crossover_keeps_distances_identical_across_all_arms():
    """All four arms share one shape; only handedness and graph differ."""
    arms = build_crossover(R_PHENYLETHANOL, S_PHENYLETHANOL, seed=SEED)
    reference = arms[0].mol
    for arm in arms[1:]:
        assert max_distance_deviation(reference, arm.mol) < 1e-6, arm.arm_id


def test_independent_embedding_would_NOT_give_a_mirror_image():
    """Why `build_crossover` reflects instead of embedding each member.

    Embedding the two enantiomers separately, even under one seed, lands them
    in different torsional minima. The geometries then differ in shape as well
    as handedness, which would confound the crossover entirely. This test pins
    that down so nobody 'simplifies' build_crossover back into two embeddings.
    """
    a = build_reference_conformer(R_PHENYLETHANOL, seed=SEED)
    b_embedded = build_reference_conformer(S_PHENYLETHANOL, seed=SEED)
    assert max_distance_deviation(a, b_embedded) > 0.1

    b_reflected = reflect_molecule(a)
    assert max_distance_deviation(a, b_reflected) < 1e-9


def test_crossover_pairs_share_coordinates_across_graphs():
    """Two geometries crossed against two graphs, not four separate shapes."""
    arms = {a.arm_id: a for a in build_crossover(R_PHENYLETHANOL, S_PHENYLETHANOL, seed=SEED)}
    import numpy as np

    from stereobench.conformer import get_coords

    assert np.allclose(get_coords(arms["a_crossed"].mol), get_coords(arms["b_native"].mol))
    assert np.allclose(get_coords(arms["b_crossed"].mol), get_coords(arms["a_native"].mol))


def test_crossover_rejects_a_non_enantiomeric_pair():
    with pytest.raises(ValueError, match="not enantiomers"):
        build_crossover(R_PHENYLETHANOL, "CC(C)O", seed=SEED)


def test_seed_makes_embedding_reproducible():
    a = build_reference_conformer(R_PHENYLETHANOL, seed=SEED)
    b = build_reference_conformer(R_PHENYLETHANOL, seed=SEED)
    assert max_distance_deviation(a, b) < 1e-9


def test_mol_pickle_roundtrip(tmp_path, mol_r):
    path = write_mol_pickle(mol_r, tmp_path / "LIG.pkl")
    assert path.name == "LIG.pkl"
    with path.open("rb") as handle:
        loaded = pickle.load(handle)
    assert loaded.GetNumAtoms() == mol_r.GetNumAtoms()
    assert conformer_labels(loaded) == conformer_labels(mol_r)
    assert loaded.GetConformer().GetProp("name") == CONFORMER_NAME


def test_unparseable_smiles_raises():
    with pytest.raises(ValueError, match="unparseable"):
        build_reference_conformer("not_a_molecule[")


def test_get_coords_requires_a_conformer():
    mol = Chem.MolFromSmiles(ACHIRAL)
    with pytest.raises(ValueError, match="no conformer"):
        get_coords(mol)


# --- over-specification regression -----------------------------------------
#
# `reflect_molecule` used to update the graph with
# `AssignStereochemistryFrom3D`, which tags every atom whose neighbours are
# distinguishable in 3D. That over-specifies relative to the input SMILES, and
# the extra tags reach the canonical SMILES, so `build_crossover`'s enantiomer
# check rejected molecules that were perfectly correct. It cost 10% of the
# pair set, silently, at run-build time.

ADAMANTYL_R = "C[C@@H](OC(=O)C12CC3CC(CC(C3)C1)C2)c1ccncc1"
ADAMANTYL_S = "C[C@H](OC(=O)C12CC3CC(CC(C3)C1)C2)c1ccncc1"


def test_reflecting_a_caged_molecule_reproduces_its_enantiomer():
    """The adamantyl bridgeheads must not acquire tags the input never had."""
    mol = build_reference_conformer(ADAMANTYL_R, seed=SEED)
    produced = Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(reflect_molecule(mol))))
    assert produced == Chem.MolToSmiles(Chem.MolFromSmiles(ADAMANTYL_S))


def test_crossover_accepts_a_caged_enantiomer_pair():
    arms = build_crossover(ADAMANTYL_R, ADAMANTYL_S, seed=SEED, axis=0)
    assert len(arms) == 4


def test_reflection_does_not_add_stereocentres_to_the_graph():
    """Mirroring changes handedness, never how specified the molecule is."""
    mol = build_reference_conformer(ADAMANTYL_R, seed=SEED)
    assert set(graph_labels(reflect_molecule(mol))) == set(graph_labels(mol))


def test_reflection_inverts_every_specified_centre():
    mol = build_reference_conformer(ADAMANTYL_R, seed=SEED)
    before, after = graph_labels(mol), graph_labels(reflect_molecule(mol))
    assert before and all(after[k] != v for k, v in before.items())


def test_the_enantiomer_guard_still_rejects_a_non_pair():
    """The fix must not have bought its coverage by weakening the check."""
    with pytest.raises(ValueError, match="not enantiomers"):
        build_crossover(ADAMANTYL_R, "C[C@H](O)c1ccccc1", seed=SEED, axis=0)
