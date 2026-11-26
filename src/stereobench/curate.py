"""Turn StereoAtlas rows into candidate pairs.

StereoAtlas is organised family-first. `FAMILY_MINMAX_FOLD.csv` gives, per
family x target x metric, the least and most potent stereoform -- but it names
them only as stereo-descriptor strings (`R,R,S`), not as compounds.
`TRUE_ISOMER_LITERATURE_CLEAN.csv` is the per-member table, so joining on
`stereochemical_family_id` + `stereochemical_variant` turns each descriptor
string back into a `compound_chembl_id`.

Three things this module will not do, by design.

It does not decide that a pair is enantiomeric. The variant strings can only
say that every centre differs, which is a *candidate* test -- meso compounds
satisfy it while being their own mirror image. Confirmation needs real SMILES
and `stereo.classify_pair`, which needs the ChEMBL pass.

It does not enforce the same-assay requirement, because no StereoAtlas file
carries an `assay_chembl_id`. `mean_activity_value` is, as the name says, a
mean -- potentially across labs and protocols, where cross-lab IC50 variance
can exceed the stereochemical effect being measured. Until the ChEMBL pass
attaches real assay ids, every row here is a candidate and nothing more.

It does not filter censored measurements, for the same reason: the qualifier
(`>`, `<`) lives in ChEMBL, not here. Censored values are the likeliest source
of the implausible fold-differences in this archive, so `plausible_value`
catches the worst of them by range alone as a stopgap.
"""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

# Values outside this range in nM are not real measurements. The archive
# contains fold-differences up to 1e20, which arise from censored or
# mis-parsed values rather than from chemistry.
MIN_PLAUSIBLE_NM = 1e-4
MAX_PLAUSIBLE_NM = 1e7

# SCHEMA.md restricts pairs to binding/potency metrics with a defined direction.
ALLOWED_ACTIVITY_TYPES = ("Ki", "IC50", "Kd")

# StereoAtlas tiers. A blank classification is the "None" tier (<2x), which is
# the near-zero arm -- the pairs where the model is *supposed* to be flat.
TAIL_TIERS = ("Strong", "Extreme")
NEAR_ZERO_TIERS = ("None",)


def pact_from_nm(value_nm: float) -> float:
    """nM to pAct (-log10 molar). Higher = stronger, matching SCHEMA.md."""
    return 9.0 - math.log10(value_nm)


def plausible_value(value: float) -> bool:
    return (
        value is not None
        and not pd.isna(value)
        and MIN_PLAUSIBLE_NM <= value <= MAX_PLAUSIBLE_NM
    )


def parse_variant(variant: str) -> list[str] | None:
    """Split a descriptor string like 'R,R,S' into centres.

    Returns None if any centre is not R or S -- undefined or mixed centres
    cannot support an inversion test.
    """
    if not isinstance(variant, str) or not variant.strip():
        return None
    centres = [c.strip() for c in variant.split(",")]
    if not centres or any(c not in ("R", "S") for c in centres):
        return None
    return centres


def is_enantiomer_candidate(variant_a: str, variant_b: str) -> bool:
    """True when every stereocentre differs between the two descriptors.

    A *candidate* test only. Necessary but not sufficient: a meso compound has
    all centres inverted relative to itself and is nonetheless identical to its
    mirror image. `stereo.classify_pair` is the real decision, and it needs
    SMILES this archive does not carry.
    """
    a, b = parse_variant(variant_a), parse_variant(variant_b)
    if a is None or b is None or len(a) != len(b):
        return False
    return all(x != y for x, y in zip(a, b))


def load_fold_table(path: str | Path) -> pd.DataFrame:
    """Family-level min/max table, with the blank tier named explicitly."""
    df = pd.read_csv(path, low_memory=False)
    # A blank classification is the "None" (<2x) tier, not missing data. Left
    # as NaN it silently drops out of every tier filter, which would quietly
    # delete the entire near-zero arm.
    df["classification"] = df["classification"].fillna("None")
    return df


def load_member_table(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, low_memory=False)


def build_variant_index(members: pd.DataFrame) -> pd.DataFrame:
    """Map (family, variant) -> compound, dropping keys that resolve to several.

    A descriptor string is *usually* unique within a family, but not always:
    about 237 of 76,302 keys name more than one compound (distinguished
    upstream only by a `stereo_form` tag the fold table does not reference).
    Those keys cannot be resolved without guessing, so they are dropped rather
    than arbitrarily picked.
    """
    index = (
        members.groupby(["stereochemical_family_id", "stereochemical_variant"])[
            "compound_chembl_id"
        ]
        .agg(["nunique", "first"])
        .reset_index()
    )
    unambiguous = index[index["nunique"] == 1].copy()
    return unambiguous.rename(columns={"first": "compound_chembl_id"})[
        ["stereochemical_family_id", "stereochemical_variant", "compound_chembl_id"]
    ]


def resolve_pairs(fold: pd.DataFrame, members: pd.DataFrame) -> pd.DataFrame:
    """Attach a ChEMBL id to each side of every fold row.

    The eutomer is the more potent member. Activity values are nM, so the
    *lower* value is the stronger binder -- `min_value` is the eutomer and
    `max_value` the distomer, which is the opposite of what the column names
    suggest in pAct space.
    """
    index = build_variant_index(members)

    out = fold.merge(
        index.rename(
            columns={
                "stereochemical_variant": "min_variant",
                "compound_chembl_id": "eutomer_chembl_id",
            }
        ),
        on=["stereochemical_family_id", "min_variant"],
        how="left",
    ).merge(
        index.rename(
            columns={
                "stereochemical_variant": "max_variant",
                "compound_chembl_id": "distomer_chembl_id",
            }
        ),
        on=["stereochemical_family_id", "max_variant"],
        how="left",
    )

    resolved = out.dropna(subset=["eutomer_chembl_id", "distomer_chembl_id"])
    # A family representative can appear on both sides when the descriptor
    # strings collide; that is not a pair.
    return resolved[
        resolved["eutomer_chembl_id"] != resolved["distomer_chembl_id"]
    ].copy()


def assign_group(classification: str) -> str | None:
    if classification in TAIL_TIERS:
        return "tail"
    if classification in NEAR_ZERO_TIERS:
        return "near_zero"
    return None


def build_candidates(fold: pd.DataFrame, members: pd.DataFrame) -> pd.DataFrame:
    """Full candidate table: resolved, plausible, enantiomer-candidate pairs.

    Everything here still needs the ChEMBL pass before it can be called a pair:
    SMILES for `classify_pair`, an assay id for the same-assay rule, and a
    qualifier to exclude censored values.
    """
    pairs = resolve_pairs(fold, members)

    pairs = pairs[pairs["activity_type"].isin(ALLOWED_ACTIVITY_TYPES)]
    pairs = pairs[
        pairs["min_value"].apply(plausible_value)
        & pairs["max_value"].apply(plausible_value)
    ]
    pairs = pairs[
        pairs.apply(
            lambda r: is_enantiomer_candidate(r["min_variant"], r["max_variant"]),
            axis=1,
        )
    ]

    pairs["group"] = pairs["classification"].apply(assign_group)
    pairs = pairs.dropna(subset=["group"])

    pairs["eutomer_pact"] = pairs["min_value"].apply(pact_from_nm)
    pairs["distomer_pact"] = pairs["max_value"].apply(pact_from_nm)
    pairs["delta_pact"] = pairs["eutomer_pact"] - pairs["distomer_pact"]
    pairs["n_stereocenters"] = pairs["min_variant"].apply(
        lambda v: len(parse_variant(v) or [])
    )

    columns = [
        "stereochemical_family_id",
        "target_chembl_id",
        "target_name",
        "organism",
        "activity_type",
        "eutomer_chembl_id",
        "distomer_chembl_id",
        "min_variant",
        "max_variant",
        "n_stereocenters",
        "eutomer_pact",
        "distomer_pact",
        "delta_pact",
        "fold_difference",
        "classification",
        "group",
    ]
    return pairs[columns].reset_index(drop=True)


def rank_for_pilot(candidates: pd.DataFrame) -> pd.DataFrame:
    """Order candidates by how easy they are to verify by hand.

    A pilot wants pairs a human can check, not a representative sample. Fewer
    stereocentres means an unambiguous CIP reading off the predicted geometry
    and fewer `partial` verdicts; a human target means a sequence is easy to
    find.

    **The two arms must sort delta in OPPOSITE directions**, and getting this
    wrong is silent. The tail arm wants the largest separations, where a flat
    prediction is unambiguous. The near-zero arm wants the SMALLEST ones --
    its whole purpose is pairs where the model is supposed to be flat. Sorting
    both descending fills the near-zero arm with pairs pressed against the 2x
    tier ceiling, which is the least near-zero sample available and tests
    nothing. A first version of this function did exactly that, and the
    symptom was a third of the pair set landing within 0.01 of log10(2).
    """
    ranked = candidates.copy()
    ranked["is_human"] = ranked["organism"].eq("Homo sapiens")
    # Sort key that is descending for the tail arm and ascending for near-zero.
    ranked["_delta_rank"] = ranked.apply(
        lambda r: -r["delta_pact"] if r["group"] == "tail" else r["delta_pact"],
        axis=1,
    )
    ranked = ranked.sort_values(
        ["group", "is_human", "n_stereocenters", "_delta_rank"],
        ascending=[True, False, True, True],
    )
    return ranked.drop(columns=["_delta_rank"]).reset_index(drop=True)


# Concentrations that recur as "highest tested" rather than as measurements.
# A compound reported inactive at the top of a dilution series often lands in
# ChEMBL with relation "=" at the cutoff, so the qualifier filter cannot see
# it. These are flagged rather than dropped: the value may be genuine, and
# that is a judgement for whoever inspects the pair.
SUSPICIOUS_CUTOFFS_NM = (
    1e4, 2e4, 2.5e4, 3e4, 5e4, 1e5, 2e5, 2.5e5, 3e5, 5e5, 1e6, 1e7,
)


def looks_like_a_cutoff(value: float, *, tolerance: float = 1e-6) -> bool:
    """True for values sitting exactly on a common assay ceiling."""
    if value is None or pd.isna(value):
        return False
    return any(
        abs(value - cutoff) <= tolerance * max(cutoff, 1.0)
        for cutoff in SUSPICIOUS_CUTOFFS_NM
    )
