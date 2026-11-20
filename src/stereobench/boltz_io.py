"""Build Boltz-2 inputs and read its affinity outputs.

Schema and output layout follow the Boltz prediction docs. Two constraints from
those docs are enforced here because violating them yields numbers that look
fine and mean nothing: affinity takes exactly one ligand chain, and the ligand
must be small.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from rdkit import Chem

# Boltz-2 affinity limits, per the prediction docs.
AFFINITY_MAX_ATOMS = 128
AFFINITY_RECOMMENDED_MAX_ATOMS = 56

# RT*ln(10) at 298 K, kcal/mol -- the docs' pIC50 -> kcal/mol factor.
_RT_LN10 = 1.364


def affinity_atom_count(smiles: str) -> int:
    """Atoms as Boltz counts them: heavy atoms plus the Hs RDKit's RemoveHs keeps."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"unparseable SMILES: {smiles!r}")
    return Chem.RemoveHs(Chem.AddHs(mol)).GetNumAtoms()


def build_yaml(
    protein_sequence: str,
    ligand_smiles: str,
    *,
    protein_id: str = "A",
    ligand_id: str = "B",
    msa_path: str | None = None,
) -> str:
    """A Boltz-2 YAML for one protein + one ligand, requesting affinity.

    `msa_path` of None means the caller intends to pass --use_msa_server. Both
    members of a pair must be given the SAME msa_path: the protein is identical
    between enantiomers, so a differing MSA would introduce a difference that
    has nothing to do with stereochemistry.
    """
    n_atoms = affinity_atom_count(ligand_smiles)
    if n_atoms > AFFINITY_MAX_ATOMS:
        raise ValueError(
            f"ligand has {n_atoms} atoms, over Boltz-2's affinity limit of "
            f"{AFFINITY_MAX_ATOMS}"
        )
    lines = [
        "version: 1",
        "sequences:",
        "  - protein:",
        f"      id: {protein_id}",
        f"      sequence: {protein_sequence}",
    ]
    if msa_path is not None:
        lines.append(f"      msa: {msa_path}")
    lines += [
        "  - ligand:",
        f"      id: {ligand_id}",
        f"      smiles: '{ligand_smiles}'",
        "properties:",
        "  - affinity:",
        f"      binder: {ligand_id}",
    ]
    return "\n".join(lines) + "\n"


@dataclass(frozen=True)
class Affinity:
    """One affinity_<name>.json.

    `pred_value` is log10(IC50) with IC50 in micromolar, so LOWER is stronger.
    Boltz emits one of these per input, not per diffusion sample; value1/value2
    are the two ensemble members.
    """

    pred_value: float
    probability_binary: float
    pred_value1: float | None = None
    pred_value2: float | None = None

    @property
    def pic50(self) -> float:
        """pIC50. Comparable to experimental pKi/pIC50: HIGHER is stronger."""
        return 6.0 - self.pred_value

    @property
    def kcal_per_mol(self) -> float:
        return self.pic50 * _RT_LN10

    @property
    def ensemble_spread(self) -> float | None:
        """|value1 - value2|. A cheap check that the two heads agree at all."""
        if self.pred_value1 is None or self.pred_value2 is None:
            return None
        return abs(self.pred_value1 - self.pred_value2)


def parse_affinity(path: str | Path) -> Affinity:
    data = json.loads(Path(path).read_text())
    return Affinity(
        pred_value=float(data["affinity_pred_value"]),
        probability_binary=float(data["affinity_probability_binary"]),
        pred_value1=_opt(data, "affinity_pred_value1"),
        pred_value2=_opt(data, "affinity_pred_value2"),
    )


def _opt(data: dict, key: str) -> float | None:
    return float(data[key]) if key in data and data[key] is not None else None


def find_affinity_json(out_dir: str | Path, name: str) -> Path:
    """Locate affinity_<name>.json under a --out_dir."""
    path = Path(out_dir) / "predictions" / name / f"affinity_{name}.json"
    if not path.exists():
        raise FileNotFoundError(path)
    return path
