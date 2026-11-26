"""Tests for the StereoAtlas -> candidate-pair path.

The sharp edge here is `rank_for_pilot`. Its two arms must sort delta in
opposite directions, and getting it wrong produces a plausible-looking table
that silently tests nothing.
"""

from __future__ import annotations

import pandas as pd
import pytest

from stereobench.curate import (
    assign_group,
    build_variant_index,
    is_enantiomer_candidate,
    looks_like_a_cutoff,
    pact_from_nm,
    parse_variant,
    plausible_value,
    rank_for_pilot,
)


def test_pact_conversion_is_nanomolar():
    assert pact_from_nm(1.0) == pytest.approx(9.0)
    assert pact_from_nm(1000.0) == pytest.approx(6.0)


def test_plausibility_rejects_the_archive_absurdities():
    assert plausible_value(50.0)
    assert not plausible_value(1.9e19)
    assert not plausible_value(None)


def test_parse_variant_rejects_undefined_centres():
    assert parse_variant("R,S,R") == ["R", "S", "R"]
    assert parse_variant("R,X") is None
    assert parse_variant("") is None


def test_enantiomer_candidate_needs_every_centre_inverted():
    assert is_enantiomer_candidate("R,R", "S,S")
    assert not is_enantiomer_candidate("R,R", "S,R")
    assert not is_enantiomer_candidate("R,R", "S,S,S")


def test_meso_style_descriptor_still_passes_the_candidate_test():
    """Why this is a candidate test and not a decision.

    A meso compound has every centre inverted relative to itself and is
    nonetheless its own mirror image. Only `stereo.classify_pair` on real
    SMILES can tell the difference.
    """
    assert is_enantiomer_candidate("R,S", "S,R")


def test_blank_tier_is_the_near_zero_arm():
    assert assign_group("None") == "near_zero"
    assert assign_group("Extreme") == "tail"
    assert assign_group("Mild") is None


def test_ambiguous_variant_keys_are_dropped_not_guessed():
    members = pd.DataFrame(
        {
            "stereochemical_family_id": ["F1", "F1", "F2"],
            "stereochemical_variant": ["R,S", "R,S", "S,R"],
            "compound_chembl_id": ["C1", "C2", "C3"],
        }
    )
    index = build_variant_index(members)
    assert list(index.compound_chembl_id) == ["C3"]


def test_cutoff_detection_catches_assay_ceilings():
    assert looks_like_a_cutoff(300000.0)
    assert looks_like_a_cutoff(1e4)
    assert not looks_like_a_cutoff(5.5)
    assert not looks_like_a_cutoff(None)


def _candidates():
    return pd.DataFrame(
        {
            "group": ["tail"] * 3 + ["near_zero"] * 3,
            "delta_pact": [1.0, 5.0, 3.0, 0.29, 0.001, 0.15],
            "organism": ["Homo sapiens"] * 6,
            "n_stereocenters": [1] * 6,
        }
    )


def test_tail_arm_takes_the_largest_deltas():
    ranked = rank_for_pilot(_candidates())
    tail = ranked[ranked.group == "tail"]
    assert list(tail.delta_pact) == [5.0, 3.0, 1.0]


def test_near_zero_arm_takes_the_SMALLEST_deltas():
    """The regression that matters.

    Sorting both arms descending fills the near-zero arm with pairs pressed
    against the 2x tier ceiling -- the least near-zero sample available. The
    first version of this did exactly that, and a third of the resulting pair
    set landed within 0.01 of log10(2).
    """
    ranked = rank_for_pilot(_candidates())
    near_zero = ranked[ranked.group == "near_zero"]
    assert list(near_zero.delta_pact) == [0.001, 0.15, 0.29]
    # The arm's top pick must be near zero, not near the tier boundary.
    assert near_zero.delta_pact.iloc[0] < 0.01


def test_ranking_does_not_leak_its_sort_key():
    assert "_delta_rank" not in rank_for_pilot(_candidates()).columns
