"""Attach a protein sequence to every target in the pair set.

Boltz needs a sequence, StereoAtlas has none, and ChEMBL target ids are not
UniProt accessions. Two hops: ChEMBL target -> UniProt accession -> sequence.

Targets that are not a SINGLE PROTEIN are dropped rather than guessed at --
complexes and protein families cannot be modelled as one chain, and a pair
whose target is a family is not a pair against a defined binding site.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from stereobench.chembl import fetch_target_uniprot, fetch_uniprot_sequence


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pairs", default="data/pairs/pilot.csv", type=Path)
    ap.add_argument("--out", default="data/pairs/target_sequences.json", type=Path)
    ap.add_argument("--max-length", type=int, default=2000,
                    help="skip very long chains; they dominate runtime")
    args = ap.parse_args()

    targets = sorted(pd.read_csv(args.pairs).target_chembl_id.unique())
    print(f"{len(targets)} distinct targets")

    sequences, skipped = {}, {}
    for i, target in enumerate(targets, 1):
        accession = fetch_target_uniprot(target)
        if accession is None:
            skipped[target] = "not_single_protein"
            continue
        sequence = fetch_uniprot_sequence(accession)
        if not sequence:
            skipped[target] = "no_sequence"
            continue
        if len(sequence) > args.max_length:
            skipped[target] = f"too_long:{len(sequence)}"
            continue
        sequences[target] = sequence
        print(f"  [{i}/{len(targets)}] {target} -> {accession} ({len(sequence)} aa)")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(sequences, indent=2))
    print(f"\nresolved {len(sequences)}, skipped {len(skipped)}")
    for target, reason in list(skipped.items())[:10]:
        print(f"  {target}: {reason}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
