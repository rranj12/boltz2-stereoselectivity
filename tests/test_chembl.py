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


# --- transient empty responses ---------------------------------------------


def test_empty_response_is_retried_when_the_page_claims_rows(monkeypatch):
    """ChEMBL sometimes answers 200 with no rows but a nonzero total_count.

    Taking that at face value silently shrinks the pair set: the pair is
    recorded as having no shared assay, which looks like a finding rather than
    a dropped request.
    """
    from stereobench import chembl

    responses = [
        {"activities": [], "page_meta": {"total_count": 2}},
        {
            "activities": [
                {
                    "assay_chembl_id": "A1",
                    "standard_type": "Ki",
                    "standard_value": "10",
                    "standard_units": "nM",
                    "standard_relation": "=",
                    "document_chembl_id": "D1",
                }
            ],
            "page_meta": {"total_count": 2},
        },
    ]
    calls = []

    def fake_get(path, params):
        calls.append(params)
        return responses[min(len(calls) - 1, len(responses) - 1)]

    monkeypatch.setattr(chembl, "_get", fake_get)
    monkeypatch.setattr(chembl.time, "sleep", lambda _: None)

    out = chembl.fetch_activities("CHEMBL1", "CHEMBL2", "Ki")

    assert len(calls) == 2
    assert len(out) == 1


def test_genuinely_empty_response_is_not_retried(monkeypatch):
    """A compound with no measurements must not cost a second request."""
    from stereobench import chembl

    calls = []

    def fake_get(path, params):
        calls.append(params)
        return {"activities": [], "page_meta": {"total_count": 0}}

    monkeypatch.setattr(chembl, "_get", fake_get)
    out = chembl.fetch_activities("CHEMBL1", "CHEMBL2", "Ki")

    assert len(calls) == 1
    assert out == []


def test_non_nanomolar_units_are_skipped(monkeypatch):
    from stereobench import chembl

    monkeypatch.setattr(
        chembl,
        "_get",
        lambda path, params: {
            "activities": [
                {
                    "assay_chembl_id": "A1",
                    "standard_type": "Ki",
                    "standard_value": "1",
                    "standard_units": "ug.mL-1",
                    "standard_relation": "=",
                    "document_chembl_id": "D1",
                }
            ],
            "page_meta": {"total_count": 1},
        },
    )
    assert chembl.fetch_activities("CHEMBL1", "CHEMBL2", "Ki") == []
