import json

import pytest
from rdkit import Chem

from stereobench import boltz_io

SMILES = "C[C@@H](c1ccc(CC(C)C)cc1)C(=O)O"
SEQ = "MKVLATGW"


def test_yaml_requests_affinity_for_the_ligand_chain():
    y = boltz_io.build_yaml(SEQ, SMILES)
    assert "properties:" in y
    assert "- affinity:" in y
    assert "binder: B" in y
    assert f"sequence: {SEQ}" in y


def test_yaml_omits_msa_key_when_using_the_server():
    assert "msa:" not in boltz_io.build_yaml(SEQ, SMILES)
    assert "msa: ./x.a3m" in boltz_io.build_yaml(SEQ, SMILES, msa_path="./x.a3m")


def test_oversized_ligand_is_rejected():
    # Affinity output would be unreliable rather than absent, so fail loudly.
    with pytest.raises(ValueError, match="affinity limit"):
        boltz_io.build_yaml(SEQ, "C" * 200)


def test_atom_count_follows_the_removehs_convention():
    # Boltz counts "heavy atoms plus the hydrogens RDKit's RemoveHs keeps".
    # RemoveHs drops ordinary hydrogens, so for a normal ligand this is just
    # the heavy-atom count -- methane is 1 atom, not 5.
    assert boltz_io.affinity_atom_count("C") == 1
    assert boltz_io.affinity_atom_count(SMILES) == Chem.MolFromSmiles(SMILES).GetNumAtoms()


def test_parse_affinity_and_pic50_direction(tmp_path):
    p = tmp_path / "affinity_x.json"
    p.write_text(
        json.dumps(
            {
                "affinity_pred_value": -3.0,
                "affinity_probability_binary": 0.91,
                "affinity_pred_value1": -3.1,
                "affinity_pred_value2": -2.9,
            }
        )
    )
    a = boltz_io.parse_affinity(p)
    assert a.pred_value == -3.0
    # log10(IC50 uM) = -3 is IC50 = 1 nM, i.e. pIC50 9.
    assert a.pic50 == pytest.approx(9.0)
    assert a.ensemble_spread == pytest.approx(0.2)


def test_stronger_binder_has_lower_pred_value_but_higher_pic50():
    strong = boltz_io.Affinity(pred_value=-3.0, probability_binary=0.9)
    weak = boltz_io.Affinity(pred_value=2.0, probability_binary=0.2)
    assert strong.pred_value < weak.pred_value
    assert strong.pic50 > weak.pic50
