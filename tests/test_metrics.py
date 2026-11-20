import pytest

from stereobench.metrics import Pair, flat_fraction, noise_sigma, sign_accuracy, summarize


def mk(pid, de, dp, se="correct", sd="flipped"):
    """Pair with experimental delta `de` and predicted delta `dp`."""
    return Pair(pid, exp_eutomer=de, exp_distomer=0.0, pred_eutomer=dp,
                pred_distomer=0.0, stereo_eutomer=se, stereo_distomer=sd)


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


def test_key_conditional_reported_for_same_geometry_pairs():
    # Both members built identically => a geometry-driven head must be flat.
    pairs = [mk("a", 2.0, 0.8, se="correct", sd="correct")]
    s = summarize(pairs, threshold=0.1)
    assert s["n_same_realized_geometry"] == 1
    assert s["mean_abs_delta_pred_same_geometry"] == pytest.approx(0.8)


def test_conditional_is_none_when_no_pair_shares_geometry():
    s = summarize([mk("a", 2.0, 1.0, se="correct", sd="flipped")], threshold=0.1)
    assert s["n_same_realized_geometry"] == 0
    assert s["mean_abs_delta_pred_same_geometry"] is None
