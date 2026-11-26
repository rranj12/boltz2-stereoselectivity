"""Tests for reading the predicted ligand back out of a CIF.

The theme here is that this module must fail loudly. A partial atom match
would leave some atoms at the origin, and those atoms would read as a
perfectly confident, perfectly wrong stereocentre -- a silent corruption of
the exact measurement the project rests on.
"""

from __future__ import annotations

import pytest
from rdkit import Chem

from stereobench.conformer import (
    build_reference_conformer,
    conformer_labels,
    get_coords,
    reflect_conformer_only,
)
from stereobench.structure import (
    ligand_from_cif,
    parse_atom_sites,
    realized_labels,
    template_conformer_labels,
    verify_prediction,
)

SMILES = "C[C@H](O)c1ccccc1"
SEED = 0xC0FFEE

CIF_HEADER = """data_test
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
"""


def write_cif(path, mol, *, comp_id="LIG", chain="B", skip_names=()):
    """Serialise a mol's heavy atoms as a minimal Boltz-shaped mmCIF."""
    heavy = Chem.RemoveHs(Chem.Mol(mol), sanitize=False)
    coords = get_coords(heavy)
    lines = [CIF_HEADER]
    serial = 1
    for atom, (x, y, z) in zip(heavy.GetAtoms(), coords):
        name = atom.GetProp("name")
        if name in skip_names:
            continue
        lines.append(
            f"ATOM {serial} {atom.GetSymbol()} {name} {comp_id} {chain} 1 "
            f"{x:.3f} {y:.3f} {z:.3f}"
        )
        serial += 1
    lines.append("#")
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture
def mol():
    return build_reference_conformer(SMILES, seed=SEED)


def test_parses_atom_sites(tmp_path, mol):
    atoms = parse_atom_sites(write_cif(tmp_path / "a.cif", mol))
    assert len(atoms) == Chem.RemoveHs(Chem.Mol(mol)).GetNumAtoms()
    assert {a.label_asym_id for a in atoms} == {"B"}


def test_column_order_is_read_not_assumed(tmp_path, mol):
    """Coordinates must be located by column name, not by position."""
    path = write_cif(tmp_path / "a.cif", mol)
    text = path.read_text()
    atoms = parse_atom_sites(path)
    # The x of the first atom must match the file, whatever the column layout.
    first_data = [l for l in text.splitlines() if l.startswith("ATOM")][0].split()
    assert atoms[0].x == pytest.approx(float(first_data[7]))


def test_roundtrip_recovers_the_geometry(tmp_path, mol):
    recovered = ligand_from_cif(write_cif(tmp_path / "a.cif", mol), mol)
    assert conformer_labels(recovered) == template_conformer_labels(mol)


def test_reads_the_handedness_that_was_written(tmp_path, mol):
    """A mirrored structure must come back mirrored, not silently corrected."""
    mirrored = reflect_conformer_only(mol)
    native_labels = realized_labels(write_cif(tmp_path / "n.cif", mol), mol)
    flipped_labels = realized_labels(write_cif(tmp_path / "m.cif", mirrored), mol)
    assert set(native_labels) == set(flipped_labels)
    assert all(native_labels[k] != flipped_labels[k] for k in native_labels)


def test_verify_prediction_calls_a_matching_structure_correct(tmp_path, mol):
    assert verify_prediction(write_cif(tmp_path / "a.cif", mol), mol).status == "correct"


def test_verify_prediction_calls_a_mirrored_structure_flipped(tmp_path, mol):
    mirrored = reflect_conformer_only(mol)
    assert verify_prediction(write_cif(tmp_path / "m.cif", mirrored), mol).status == "flipped"


def test_missing_heavy_atom_raises_rather_than_guessing(tmp_path, mol):
    """The failure that matters: a partial match must never score."""
    heavy = Chem.RemoveHs(Chem.Mol(mol), sanitize=False)
    dropped = next(a.GetProp("name") for a in heavy.GetAtoms() if a.GetSymbol() == "C")
    path = write_cif(tmp_path / "partial.cif", mol, skip_names=(dropped,))
    with pytest.raises(ValueError, match="absent from the CIF|mismatch"):
        ligand_from_cif(path, mol)


def test_no_ligand_chain_raises(tmp_path):
    path = tmp_path / "protein_only.cif"
    path.write_text(CIF_HEADER + "ATOM 1 C CA ALA A 1 0.000 0.000 0.000\n#\n")
    with pytest.raises(ValueError, match="no non-protein chain"):
        ligand_from_cif(path, build_reference_conformer(SMILES, seed=SEED))


def test_picks_the_largest_non_protein_chain(tmp_path, mol):
    """An ion or additive must not be mistaken for the ligand."""
    path = write_cif(tmp_path / "a.cif", mol)
    text = path.read_text().rstrip().removesuffix("#").rstrip()
    text += "\nATOM 999 ZN ZN ZN C 1 50.000 50.000 50.000\n#\n"
    path.write_text(text)
    recovered = ligand_from_cif(path, mol)
    assert conformer_labels(recovered) == template_conformer_labels(mol)


def test_label_spaces_must_not_be_compared_across_graphs(mol):
    """The hazard `heavy_atom_template` exists to prevent.

    Hydrogen-bearing and hydrogen-stripped graphs rank atoms differently, so
    labels from one are meaningless against the other. They agree on the
    verdict (S here) but not on the key, which is exactly the kind of mismatch
    that silently reads as "the centres do not correspond".
    """
    with_hydrogens = conformer_labels(mol)
    without = template_conformer_labels(mol)
    assert set(with_hydrogens.values()) == set(without.values())
    assert set(with_hydrogens) != set(without)
