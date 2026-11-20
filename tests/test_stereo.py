import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from stereobench import stereo

# (S)- and (R)-ibuprofen; (S) is the active eutomer.
S_IBU = "C[C@@H](c1ccc(CC(C)C)cc1)C(=O)O"
R_IBU = "C[C@H](c1ccc(CC(C)C)cc1)C(=O)O"
# Thalidomide, two stereocentre-bearing forms.
S_THAL = "O=C1CC[C@H](N2C(=O)c3ccccc3C2=O)C(=O)N1"
R_THAL = "O=C1CC[C@@H](N2C(=O)c3ccccc3C2=O)C(=O)N1"
# Threonine: two centres. Inverting ONE gives a diastereomer.
L_THR = "C[C@@H](O)[C@H](N)C(=O)O"
THR_DIA = "C[C@H](O)[C@H](N)C(=O)O"
# meso-tartaric acid: its own mirror image.
MESO_TARTARIC = "O[C@H](C(=O)O)[C@@H](O)C(=O)O"


def test_enantiomers_detected():
    assert stereo.classify_pair(S_IBU, R_IBU) == "enantiomers"
    assert stereo.classify_pair(S_THAL, R_THAL) == "enantiomers"


def test_identical_detected():
    assert stereo.classify_pair(S_IBU, S_IBU) == "identical"


def test_single_inversion_is_diastereomer_not_enantiomer():
    assert stereo.classify_pair(L_THR, THR_DIA) == "diastereomers"


def test_meso_compound_is_identical_to_its_mirror():
    # A naive "all centres inverted => enantiomers" rule gets this wrong.
    assert stereo.classify_pair(MESO_TARTARIC, stereo.mirror(MESO_TARTARIC)) == "identical"


def test_different_constitution_rejected():
    assert stereo.classify_pair(S_IBU, S_THAL) == "different_constitution"


def test_mirror_is_an_involution():
    assert stereo.mirror(stereo.mirror(S_IBU)) == stereo.canonical(S_IBU)


def test_enantiomers_share_flat_identity_but_differ_with_stereo():
    assert stereo.flat_canonical(S_IBU) == stereo.flat_canonical(R_IBU)
    assert stereo.canonical(S_IBU) != stereo.canonical(R_IBU)


def test_cip_labels_differ_between_enantiomers():
    a, b = stereo.labels_from_smiles(S_IBU), stereo.labels_from_smiles(R_IBU)
    assert set(a) == set(b)          # same centre, same canonical key
    assert a != b                    # opposite assignment
    assert set(a.values()) | set(b.values()) == {"R", "S"}


def _embed(smiles: str) -> Chem.Mol:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=0xF00D) == 0
    return Chem.RemoveHs(mol)


@pytest.mark.parametrize("smiles", [S_IBU, R_IBU, L_THR])
def test_labels_from_3d_recovers_requested_stereo(smiles):
    # RDKit's embedder enforces chirality, so geometry must agree with the
    # input SMILES. This is the same comparison we run on Boltz-2's output.
    assert stereo.labels_from_3d(_embed(smiles)) == stereo.labels_from_smiles(smiles)


def test_check_flags_correct_and_flipped():
    s = stereo.labels_from_smiles(S_IBU)
    r = stereo.labels_from_smiles(R_IBU)
    assert stereo.check(s, s).status == "correct"
    assert stereo.check(s, r).status == "flipped"
    assert stereo.check(s, r).n_inverted == 1


def test_check_flags_partial_inversion():
    req = stereo.labels_from_smiles(L_THR)
    got = stereo.labels_from_smiles(THR_DIA)
    assert stereo.check(req, got).status == "partial"


def test_check_is_ambiguous_when_a_centre_is_unreadable():
    req = stereo.labels_from_smiles(L_THR)
    assert stereo.check(req, {}).status == "ambiguous"
    assert stereo.check({}, {}).status == "ambiguous"
