"""Build the candidate pair table from StereoAtlas. No network, no GPU."""

from __future__ import annotations

import argparse
from pathlib import Path

from stereobench.curate import (
    build_candidates,
    load_fold_table,
    load_member_table,
    rank_for_pilot,
)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--atlas", default="data/external/stereoatlas", type=Path)
    ap.add_argument("--out", default="data/pairs/candidates.csv", type=Path)
    args = ap.parse_args()

    fold = load_fold_table(args.atlas / "FAMILY_MINMAX_FOLD.csv")
    members = load_member_table(args.atlas / "TRUE_ISOMER_LITERATURE_CLEAN.csv")

    candidates = rank_for_pilot(build_candidates(fold, members))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    candidates.to_csv(args.out, index=False)

    print(f"fold rows in              : {len(fold):,}")
    print(f"candidate pairs out       : {len(candidates):,}")
    print(f"  tail (Strong/Extreme)   : {(candidates.group == 'tail').sum():,}")
    print(f"  near_zero (<2x)         : {(candidates.group == 'near_zero').sum():,}")
    print(f"  distinct targets        : {candidates.target_chembl_id.nunique():,}")
    print(f"  human targets           : {candidates.is_human.sum():,}")
    print(f"  single stereocentre     : {(candidates.n_stereocenters == 1).sum():,}")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
