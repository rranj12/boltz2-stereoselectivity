"""Tests for the assay-selection rule.

`find_shared_assay` is the only place in the pipeline that makes a *choice*
between equally valid measurements, which makes it the only place a selection
bias can enter. These tests are pure: they build Measurement lists directly and
never touch the network.
"""

from __future__ import annotations

import math

import pytest

from stereobench.chembl import Measurement, find_shared_assay


def m(assay: str, value: float, *, relation: str | None = "=") -> Measurement:
    return Measurement(
        assay_chembl_id=assay,
        standard_type="Ki",
        standard_value=value,
        standard_units="nM",
        standard_relation=relation,
        document_chembl_id="CHEMBL_DOC",
    )


def test_no_shared_assay_returns_none():
    assert find_shared_assay([m("A1", 10)], [m("A2", 100)]) is None


def test_shared_assay_is_found():
    shared = find_shared_assay([m("A1", 10)], [m("A1", 100)])
    assert shared is not None
    assert shared.eutomer.assay_chembl_id == "A1"
    assert shared.n_shared_assays == 1
    assert shared.n_measurements == 2


def test_censored_measurements_cannot_form_a_pair():
    """A `>` reading is a bound, not a value, and must not be collapsed in."""
    assert find_shared_assay([m("A1", 10)], [m("A1", 10000, relation=">")]) is None


def test_nonpositive_values_are_dropped():
    """log() of a zero or negative value is not a number we can average."""
    assert find_shared_assay([m("A1", 10)], [m("A1", 0)]) is None


def test_replicates_collapse_geometrically():
    """Arithmetic mean would be dragged upward by the weakest replicate."""
    shared = find_shared_assay([m("A1", 10), m("A1", 1000)], [m("A1", 50)])
    assert shared.eutomer.standard_value == pytest.approx(100.0)  # not 505
    assert shared.n_measurements == 3


def test_collapsed_replicates_are_marked_uncensored():
    shared = find_shared_assay([m("A1", 10), m("A1", 1000)], [m("A1", 50)])
    assert not shared.eutomer.is_censored


def test_picks_the_BEST_SUPPORTED_assay_not_the_most_separated():
    """The regression that matters.

    A1 is backed by four measurements and shows a 2x separation. A2 is backed
    by two and shows 1000x. The old rule took A2 -- the assay that most
    flattered the hypothesis under test. Replicate count is independent of the
    outcome; separation is the outcome.
    """
    eutomer = [m("A1", 100), m("A1", 100), m("A2", 10)]
    distomer = [m("A1", 200), m("A1", 200), m("A2", 10000)]

    shared = find_shared_assay(eutomer, distomer)

    assert shared.eutomer.assay_chembl_id == "A1"
    assert shared.n_shared_assays == 2
    assert shared.n_measurements == 4
    separation = abs(
        math.log10(shared.distomer.standard_value)
        - math.log10(shared.eutomer.standard_value)
    )
    assert separation == pytest.approx(math.log10(2), abs=1e-9)


def test_ties_break_deterministically_on_assay_id():
    """Set iteration order varies between runs; the pair set must not."""
    eutomer = [m("B9", 10), m("A1", 10)]
    distomer = [m("B9", 20), m("A1", 20)]
    picks = {
        find_shared_assay(eutomer, distomer).eutomer.assay_chembl_id
        for _ in range(20)
    }
    assert picks == {"A1"}


def test_support_counts_both_sides():
    shared = find_shared_assay(
        [m("A1", 10), m("A1", 20), m("A1", 30)], [m("A1", 40)]
    )
    assert shared.n_measurements == 4
