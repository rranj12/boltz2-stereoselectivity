import pytest

from stereobench.metrics import Pair, flat_fraction, noise_sigma, sign_accuracy, summarize


def mk(pid, de, dp, se="correct", sd="flipped", ce=None, cd=None):
    """Pair with experimental delta `de` and predicted delta `dp`.

    `ce`/`cd` are the reference-conformer tags. Left unset they take the
    dataclass defaults, which differ -- the untreated state, where each member
    was embedded from its own SMILES.
    """
    kwargs = {}
    if ce is not None:
        kwargs["conformer_eutomer"] = ce
    if cd is not None:
        kwargs["conformer_distomer"] = cd
    return Pair(pid, exp_eutomer=de, exp_distomer=0.0, pred_eutomer=dp,
                pred_distomer=0.0, stereo_eutomer=se, stereo_distomer=sd, **kwargs)


def test_sign_accuracy_counts_direction_only():
    pairs = [mk("a", 2.0, 0.1), mk("b", 2.0, -0.1)]
    assert sign_accuracy(pairs) == pytest.approx(0.5)


def test_sign_accuracy_is_none_without_scorable_pairs():
    assert sign_accuracy([]) is None
    assert sign_accuracy([mk("a", 2.0, 0.0)]) is None


def test_flat_fraction_uses_the_threshold():
    pairs = [mk("a", 2.0, 0.05), mk("b", 2.0, 1.5)]
    assert flat_fraction(pairs, threshold=0.1) == pytest.approx(0.5)
    assert flat_fraction(pairs, threshold=2.0) == pytest.approx(1.0)


def test_noise_sigma_needs_replicates():
    assert noise_sigma([1.0, 1.0, 1.0]) == 0.0
    with pytest.raises(ValueError):
        noise_sigma([1.0])


def test_ambiguous_pairs_are_excluded_not_scored():
    pairs = [mk("a", 2.0, 1.0), mk("b", 2.0, 1.0, se="ambiguous")]
    s = summarize(pairs, threshold=0.1)
    assert s["n_pairs"] == 2
    assert s["n_usable"] == 1
    assert s["n_ambiguous"] == 1


def test_same_realized_geometry_requires_matching_status():
    assert mk("a", 2.0, 1.0, se="correct", sd="correct").same_realized_geometry
    assert mk("b", 2.0, 1.0, se="flipped", sd="flipped").same_realized_geometry
    # One correct, one flipped => the two members have DIFFERENT geometry.
    assert not mk("c", 2.0, 1.0, se="correct", sd="flipped").same_realized_geometry
    assert not mk("d", 2.0, 1.0, se="correct", sd="ambiguous").same_realized_geometry


def test_matched_geometry_alone_is_not_enough():
    """Matched poses still leave the reference-conformer channel open.

    Both members built identically, but each was embedded from its own SMILES,
    so their reference conformers are mirror images. A non-zero delta here has
    two explanations and the pair cannot settle the question.
    """
    pair = mk("a", 2.0, 0.8, se="correct", sd="correct")
    assert pair.same_realized_geometry
    assert not pair.same_reference_conformer
    assert not pair.pose_channel_isolated
    assert pair.conditional_state == "confounded_by_conformer"


def test_pose_channel_isolated_requires_both_channels_matched():
    pair = mk("a", 2.0, 0.8, se="correct", sd="correct", ce="shared", cd="shared")
    assert pair.pose_channel_isolated
    assert pair.conditional_state == "pose_channel_isolated"


def test_differing_geometry_is_its_own_state():
    pair = mk("a", 2.0, 1.0, se="correct", sd="flipped", ce="shared", cd="shared")
    assert not pair.same_realized_geometry
    assert pair.conditional_state == "geometry_differs"


def test_summarize_splits_the_three_conditional_states():
    pairs = [
        mk("differs", 2.0, 1.0, se="correct", sd="flipped"),
        mk("confounded", 2.0, 0.8, se="correct", sd="correct"),
        mk("isolated", 2.0, 0.4, se="correct", sd="correct", ce="shared", cd="shared"),
    ]
    s = summarize(pairs, threshold=0.1)
    assert s["n_by_conditional_state"] == {
        "geometry_differs": 1,
        "confounded_by_conformer": 1,
        "pose_channel_isolated": 1,
    }
    assert s["mean_abs_delta_pred_by_state"]["confounded_by_conformer"] == pytest.approx(0.8)
    assert s["n_pose_channel_isolated"] == 1
    assert s["mean_abs_delta_pred_isolated"] == pytest.approx(0.4)


def test_isolated_delta_is_none_without_a_substitution_arm():
    """The default for an ordinary run, and it must not read as 'no effect'."""
    s = summarize([mk("a", 2.0, 1.0, se="correct", sd="correct")], threshold=0.1)
    assert s["n_pose_channel_isolated"] == 0
    assert s["mean_abs_delta_pred_isolated"] is None
