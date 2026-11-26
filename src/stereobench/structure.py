"""Read the predicted ligand back out of a Boltz-2 CIF.

This is the step between "Boltz produced a structure" and "we can say what
stereochemistry it built". A CIF carries atoms and coordinates but no bond
orders, and RDKit cannot assign CIP labels without a bond graph. So rather than
inferring bonds from distances -- which is guesswork, and guesswork precisely
at the geometry we are trying to measure -- the coordinates are mapped onto a
template mol whose bonds are already known.

The template is the same mol that generated the input, so the mapping is by
atom NAME, not by index. Boltz names ligand atoms by canonical rank
(`parse/schema.py:1246-1257`) and preserves those names into the CIF, which
makes the join exact and order-independent.

What this cannot do
-------------------
It reads the *structure* model's output. Per the README's finding 2, the
affinity head scores a different, re-diffused structure that is never written
to disk. So a stereo verdict from here is checkpoint 2, not checkpoint 3, and
the distinction matters for any claim about what the head saw.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from rdkit import Chem

from stereobench.stereo import Labels, StereoCheck, check, labels_from_3d, labels_from_mol


@dataclass(frozen=True)
class CifAtom:
    label_asym_id: str
    label_comp_id: str
    label_atom_id: str
    element: str
    x: float
    y: float
    z: float


def parse_atom_sites(path: str | Path) -> list[CifAtom]:
    """Pull the `_atom_site` loop out of an mmCIF.

    Deliberately minimal: Boltz writes a regular, single-block mmCIF and a full
    parser would be a dependency for no gain. Reads the loop header to find
    column positions rather than assuming a fixed order, since that ordering is
    not guaranteed.
    """
    lines = Path(path).read_text().splitlines()
    columns: list[str] = []
    atoms: list[CifAtom] = []
    in_loop = False

    for line in lines:
        stripped = line.strip()
        if stripped.startswith("_atom_site."):
            columns.append(stripped.split(".", 1)[1].split()[0])
            in_loop = True
            continue
        if not in_loop:
            continue
        if not stripped or stripped.startswith(("#", "loop_", "_", "data_")):
            if columns and atoms:
                break
            continue

        fields = stripped.split()
        if len(fields) < len(columns):
            continue
        row = dict(zip(columns, fields))
        try:
            atoms.append(
                CifAtom(
                    label_asym_id=row.get("label_asym_id", ""),
                    label_comp_id=row.get("label_comp_id", ""),
                    label_atom_id=row.get("label_atom_id", "").strip('"'),
                    element=row.get("type_symbol", ""),
                    x=float(row["Cartn_x"]),
                    y=float(row["Cartn_y"]),
                    z=float(row["Cartn_z"]),
                )
            )
        except (KeyError, ValueError):
            continue
    return atoms


def ligand_atoms(atoms: list[CifAtom], *, chain_id: str | None = None) -> list[CifAtom]:
    """Atoms belonging to the ligand chain.

    With `chain_id` given, selects it directly. Otherwise picks the chain whose
    residue is not a standard amino acid -- which is how the ligand is
    distinguished when the chain labelling in the output does not match the
    input ids.
    """
    if chain_id is not None:
        return [a for a in atoms if a.label_asym_id == chain_id]

    standard = set(Chem.rdchem.__dict__.get("_dummy", [])) or {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "SEC", "PYL", "UNK", "HOH",
    }
    candidates = [a for a in atoms if a.label_comp_id.upper() not in standard]
    if not candidates:
        raise ValueError("no non-protein chain found in CIF")

    by_chain: dict[str, list[CifAtom]] = {}
    for atom in candidates:
        by_chain.setdefault(atom.label_asym_id, []).append(atom)
    # The largest non-protein chain is the ligand; anything smaller is an ion
    # or a crystallisation additive.
    return max(by_chain.values(), key=len)


def ligand_from_cif(
    cif_path: str | Path, template: Chem.Mol, *, chain_id: str | None = None
) -> Chem.Mol:
    """Template mol carrying the CIF's coordinates, matched by atom name.

    Raises rather than guessing if the names do not line up. A silent partial
    match would produce a mol with some atoms at the origin, and those atoms
    would read as a perfectly confident, perfectly wrong stereocentre.
    """
    atoms = ligand_atoms(parse_atom_sites(cif_path), chain_id=chain_id)
    coords_by_name = {a.label_atom_id: (a.x, a.y, a.z) for a in atoms}

    out = Chem.Mol(template)
    out.RemoveAllConformers()
    conformer = Chem.Conformer(out.GetNumAtoms())

    missing = []
    for atom in out.GetAtoms():
        if not atom.HasProp("name"):
            missing.append(f"template atom {atom.GetIdx()} has no name")
            continue
        name = atom.GetProp("name")
        if name not in coords_by_name:
            # Hydrogens are routinely absent from predicted output; they carry
            # no CIP information, so their loss is harmless.
            if atom.GetAtomicNum() == 1:
                continue
            missing.append(name)
            continue
        x, y, z = coords_by_name[name]
        conformer.SetAtomPosition(atom.GetIdx(), Chem.rdGeometry.Point3D(x, y, z))

    if missing:
        raise ValueError(
            f"{len(missing)} template atoms absent from the CIF "
            f"(first few: {missing[:5]}); refusing to assign stereo from a "
            "partial match"
        )

    heavy_in_cif = sum(1 for a in atoms if a.element.upper() != "H")
    heavy_in_template = sum(1 for a in out.GetAtoms() if a.GetAtomicNum() != 1)
    if heavy_in_cif != heavy_in_template:
        raise ValueError(
            f"heavy-atom count mismatch: CIF {heavy_in_cif}, "
            f"template {heavy_in_template}"
        )

    # Drop hydrogens before reading stereo: predicted H positions are the least
    # reliable coordinates in the file, and CIP does not need them.
    out.AddConformer(conformer, assignId=True)
    return Chem.RemoveHs(out, sanitize=False)


def heavy_atom_template(template: Chem.Mol) -> Chem.Mol:
    """The template reduced to the graph `ligand_from_cif` returns.

    Canonical ranks are computed over whatever graph they are given, so labels
    from a hydrogen-bearing mol are keyed into a different space than labels
    from a hydrogen-stripped one and the two CANNOT be compared. Everything
    that compares predicted against requested stereochemistry must pass through
    here first, or it will compare rank 18 against rank 8 and conclude the
    centres do not match.
    """
    return Chem.RemoveHs(Chem.Mol(template), sanitize=False)


def template_labels(template: Chem.Mol) -> Labels:
    """Requested CIP labels, keyed to match `realized_labels` output."""
    graph_only = heavy_atom_template(template)
    graph_only.RemoveAllConformers()
    return labels_from_mol(graph_only)


def template_conformer_labels(template: Chem.Mol) -> Labels:
    """Reference-conformer CIP labels, keyed to match `realized_labels`.

    Distinct from `template_labels` whenever the template carries a
    deliberately mismatched conformer -- which is the whole point of a crossed
    arm.
    """
    return labels_from_3d(heavy_atom_template(template))


def realized_labels(cif_path: str | Path, template: Chem.Mol, **kwargs) -> Labels:
    """CIP labels of the geometry Boltz actually built.

    Keyed in the heavy-atom graph's rank space. Compare against
    `template_labels` or `template_conformer_labels`, never against
    `conformer.conformer_labels` on a hydrogen-bearing mol.
    """
    return labels_from_3d(ligand_from_cif(cif_path, template, **kwargs))


def verify_prediction(
    cif_path: str | Path, template: Chem.Mol, **kwargs
) -> StereoCheck:
    """Requested vs realized, as `stereo.check` classifies it.

    `requested` is read off the template's graph, not its conformer, so this
    answers "did Boltz build what the SMILES asked for" even when the template
    carries a deliberately mismatched reference conformer.
    """
    return check(
        template_labels(template), realized_labels(cif_path, template, **kwargs)
    )


def rmsd_to_template(cif_path: str | Path, template: Chem.Mol, **kwargs) -> float:
    """Coordinate RMSD between the predicted ligand and the template conformer.

    Not a quality metric -- the two are in different frames. It exists to catch
    the degenerate case where a predicted pose is essentially the input
    conformer unchanged.
    """
    predicted = ligand_from_cif(cif_path, template, **kwargs)
    reference = Chem.RemoveHs(Chem.Mol(template), sanitize=False)
    a = np.array(predicted.GetConformer().GetPositions())
    b = np.array(reference.GetConformer().GetPositions())
    a = a - a.mean(axis=0)
    b = b - b.mean(axis=0)
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))
