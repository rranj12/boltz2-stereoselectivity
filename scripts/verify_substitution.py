"""Does conformer substitution actually reach Boltz? Run before trusting it.

The crossover arm depends on one claim: that a mol we pickle into Boltz's
molecule cache is the mol Boltz featurises, with OUR coordinates rather than a
conformer it embeds itself. That claim is cheap to test and expensive to get
wrong -- a silent fallback to ETKDG would leave every crossed arm identical to
its native counterpart and the experiment would quietly measure nothing.

Checks, in order of how badly each would break the design:

1. Boltz's own `load_molecules` reads our pickle back with atom names and the
   conformer intact. RDKit drops properties on pickle unless told otherwise.
2. The coordinates that survive are ours, to numerical precision.
3. `parse_ccd_residue` accepts it and carries those coordinates into the
   parsed residue -- the structure the featuriser reads `ref_pos` from.
4. Native and reflected mols come back with OPPOSITE handedness. If this fails
   the two arms are the same input and nothing is being tested.

Run with the Boltz environment, pointing PYTHONPATH at this package's src:

    PYTHONPATH=src path/to/boltz-venv/bin/python scripts/verify_substitution.py
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

from stereobench.conformer import (
    build_reference_conformer,
    conformer_labels,
    get_coords,
    reflect_conformer_only,
    write_mol_pickle,
)

SMILES = "C[C@H](O)c1ccccc1"
CODE_NATIVE, CODE_REFLECTED = "SB900", "SB901"

PASS, FAIL = "  PASS", "  FAIL"


def main() -> int:
    try:
        from boltz.data.mol import load_molecules
        from boltz.data.parse.schema import get_conformer, parse_ccd_residue
    except ImportError:
        print("boltz is not importable in this interpreter.", file=sys.stderr)
        print(__doc__.strip().splitlines()[-1], file=sys.stderr)
        return 2

    failures = 0

    def report(label: str, ok: bool, detail: str = "") -> None:
        nonlocal failures
        print(f"{PASS if ok else FAIL}  {label}{(' -- ' + detail) if detail else ''}")
        if not ok:
            failures += 1

    native = build_reference_conformer(SMILES, seed=0xC0FFEE)
    reflected = reflect_conformer_only(native)

    with tempfile.TemporaryDirectory() as tmp:
        mol_dir = Path(tmp)
        write_mol_pickle(native, mol_dir / f"{CODE_NATIVE}.pkl")
        write_mol_pickle(reflected, mol_dir / f"{CODE_REFLECTED}.pkl")

        print("1. round trip through Boltz's own loader")
        loaded = load_molecules(str(mol_dir), [CODE_NATIVE, CODE_REFLECTED])
        report("both codes load", set(loaded) == {CODE_NATIVE, CODE_REFLECTED})

        loaded_native = loaded[CODE_NATIVE]
        report(
            "atom names survive pickling",
            all(a.HasProp("name") for a in loaded_native.GetAtoms()),
            "Boltz matches atoms by name; without these, parsing fails",
        )
        report(
            "exactly one conformer",
            loaded_native.GetNumConformers() == 1,
            "the featuriser picks one at random",
        )

        print("\n2. coordinates are ours, not re-embedded")
        report(
            "native coordinates preserved",
            np.allclose(get_coords(loaded_native), get_coords(native), atol=1e-6),
        )
        report(
            "reflected coordinates preserved",
            np.allclose(
                get_coords(loaded[CODE_REFLECTED]), get_coords(reflected), atol=1e-6
            ),
        )

        print("\n3. the parser carries them into the residue")
        parsed_ok = True
        try:
            residue = parse_ccd_residue(CODE_NATIVE, loaded_native, 0)
            parsed_ok = residue is not None
        except Exception as exc:  # noqa: BLE001 - reporting, not handling
            report("parse_ccd_residue accepts the mol", False, repr(exc))
            parsed_ok = False
        else:
            report("parse_ccd_residue accepts the mol", parsed_ok)

        if parsed_ok:
            conformer = get_conformer(loaded_native)
            report(
                "conformer retrievable by Boltz's own lookup",
                conformer is not None,
            )
            parsed_coords = np.array(
                [[a.conformer[0], a.conformer[1], a.conformer[2]] for a in residue.atoms]
            )
            ours = get_coords(loaded_native)
            # parse_ccd_residue drops leaving/implicit atoms, so compare the
            # overlap by position count rather than demanding equal shapes.
            n = min(len(parsed_coords), len(ours))
            report(
                "parsed residue carries our coordinates",
                bool(np.allclose(np.sort(parsed_coords[:n], axis=0),
                                 np.sort(ours[:n], axis=0), atol=1e-3)),
                f"compared {n} atoms",
            )

        print("\n4. the two arms are genuinely different inputs")
        native_labels = conformer_labels(loaded[CODE_NATIVE])
        reflected_labels = conformer_labels(loaded[CODE_REFLECTED])
        report("native has readable stereocentres", bool(native_labels))
        report(
            "handedness is opposite between arms",
            bool(native_labels)
            and set(native_labels) == set(reflected_labels)
            and all(native_labels[k] != reflected_labels[k] for k in native_labels),
            f"{native_labels} vs {reflected_labels}",
        )

    print(
        f"\n{'ALL CHECKS PASSED' if not failures else f'{failures} CHECK(S) FAILED'}"
    )
    if failures:
        print(
            "Substitution is not safe to rely on. The crossover arm needs a "
            "patched Boltz instead -- capture point is boltz2.py:624."
        )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
