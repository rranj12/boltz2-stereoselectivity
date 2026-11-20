"""Stereo assignment and pair classification.

Two jobs:
  1. What stereochemistry did we ask for, and what did Boltz-2 actually build?
  2. Are two molecules enantiomers, diastereomers, or neither?

Atom indices are not comparable between two RDKit mols built from different
SMILES, so CIP labels are keyed by graph-canonical rank (computed with
includeChirality=False, so the key is stable across stereoisomers of one
constitution).
"""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

Labels = dict[int, str]  # canonical rank -> "R" | "S"


def _canon_ranks(mol: Chem.Mol) -> list[int]:
    return list(Chem.CanonicalRankAtoms(mol, includeChirality=False, breakTies=True))


def labels_from_mol(mol: Chem.Mol) -> Labels:
    """CIP labels keyed by graph-canonical rank. Uses the accurate CIP labeller."""
    mol = Chem.Mol(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    rdCIPLabeler.AssignCIPLabels(mol)
    ranks = _canon_ranks(mol)
    out: Labels = {}
    for atom in mol.GetAtoms():
        if atom.HasProp("_CIPCode"):
            out[ranks[atom.GetIdx()]] = atom.GetProp("_CIPCode")
    return out


def labels_from_smiles(smiles: str) -> Labels:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"unparseable SMILES: {smiles!r}")
    return labels_from_mol(mol)


def labels_from_3d(mol: Chem.Mol) -> Labels:
    """CIP labels read off a conformer, ignoring any stereo already on the graph.

    `mol` must carry a 3D conformer. This is how we find out what geometry
    Boltz-2 actually produced, as opposed to what the input SMILES requested.
    """
    if mol.GetNumConformers() == 0:
        raise ValueError("mol has no conformer")
    mol = Chem.Mol(mol)
    Chem.RemoveStereochemistry(mol)
    Chem.AssignStereochemistryFrom3D(mol)
    return labels_from_mol(mol)


def mirror(smiles: str) -> str:
    """Canonical SMILES of the mirror image: every tetrahedral centre inverted.

    Reflection does not affect double-bond (E/Z) stereo, so only the
    tetrahedral tags are flipped.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"unparseable SMILES: {smiles!r}")
    flip = {
        Chem.ChiralType.CHI_TETRAHEDRAL_CW: Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
        Chem.ChiralType.CHI_TETRAHEDRAL_CCW: Chem.ChiralType.CHI_TETRAHEDRAL_CW,
    }
    for atom in mol.GetAtoms():
        tag = atom.GetChiralTag()
        if tag in flip:
            atom.SetChiralTag(flip[tag])
    return Chem.MolToSmiles(mol)


def canonical(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"unparseable SMILES: {smiles!r}")
    return Chem.MolToSmiles(mol)


def flat_canonical(smiles: str) -> str:
    """Canonical SMILES with all stereochemistry stripped -- the 2D identity."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"unparseable SMILES: {smiles!r}")
    Chem.RemoveStereochemistry(mol)
    return Chem.MolToSmiles(mol)


def classify_pair(smiles_a: str, smiles_b: str) -> str:
    """One of: identical, enantiomers, diastereomers, different_constitution.

    Enantiomers are defined exactly -- the mirror image of A is B -- rather than
    by counting inverted centres, so meso compounds fall out as `identical`
    instead of being miscalled. StereoAtlas-style archives group by stripped
    SMILES only, which mixes enantiomers and diastereomers together, so this
    check has to run over any candidate set regardless of its source.
    """
    ca, cb = canonical(smiles_a), canonical(smiles_b)
    if ca == cb:
        return "identical"
    if flat_canonical(smiles_a) != flat_canonical(smiles_b):
        return "different_constitution"
    if mirror(smiles_a) == cb:
        return "enantiomers"
    return "diastereomers"


@dataclass(frozen=True)
class StereoCheck:
    """Comparison of requested vs realized stereochemistry at every centre."""

    status: str  # correct | flipped | partial | ambiguous
    requested: Labels
    realized: Labels

    @property
    def n_inverted(self) -> int:
        return sum(
            1 for k, v in self.requested.items() if k in self.realized and self.realized[k] != v
        )


def check(requested: Labels, realized: Labels) -> StereoCheck:
    """Classify a prediction as correct / flipped / partial / ambiguous.

    `ambiguous` means we could not read every requested centre off the geometry,
    so the pair must be excluded rather than scored -- an unreadable centre is
    not evidence either way.
    """
    if not requested or set(requested) != set(realized):
        return StereoCheck("ambiguous", requested, realized)
    inverted = [k for k, v in requested.items() if realized[k] != v]
    if not inverted:
        status = "correct"
    elif len(inverted) == len(requested):
        status = "flipped"
    else:
        status = "partial"
    return StereoCheck(status, requested, realized)
