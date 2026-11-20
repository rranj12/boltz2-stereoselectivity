"""Paired metrics.

Everything is in pAct space (higher = stronger) so predicted and experimental
differences carry the same sign convention. Boltz's pred_value runs the other
way, so convert with Affinity.pic50 before coming in here.

The headline question is a null -- does |predicted delta| collapse to zero --
and a null needs a noise floor or it cannot be read. Boltz-2 exposes no seed
flag, so the floor is measured by rerunning identical inputs: see noise_sigma.
"""

from __future__ import annotations

import statistics
from dataclasses import dataclass


@dataclass(frozen=True)
class Pair:
    """One eutomer/distomer pair. pAct values: higher = stronger."""

    pair_id: str
    exp_eutomer: float
    exp_distomer: float
    pred_eutomer: float
    pred_distomer: float
    # From stereo.check(): correct | flipped | partial | ambiguous
    stereo_eutomer: str = "ambiguous"
    stereo_distomer: str = "ambiguous"

    @property
    def delta_exp(self) -> float:
        return self.exp_eutomer - self.exp_distomer

    @property
    def delta_pred(self) -> float:
        return self.pred_eutomer - self.pred_distomer

    @property
    def same_realized_geometry(self) -> bool:
        """True when Boltz built the same stereochemistry for both members.

        This is the conditional that matters: if the realized geometries are
        identical, a geometry-driven affinity head MUST give |delta_pred| ~ 0.
        A non-zero delta here means the number came from somewhere other than
        the coordinates.
        """
        return self.stereo_eutomer == self.stereo_distomer and self.stereo_eutomer in (
            "correct",
            "flipped",
        )

    @property
    def usable(self) -> bool:
        return "ambiguous" not in (self.stereo_eutomer, self.stereo_distomer)


def noise_sigma(replicate_preds: list[float]) -> float:
    """Stdev of repeated runs of ONE input. The floor every delta is read against."""
    if len(replicate_preds) < 2:
        raise ValueError("need at least 2 replicates")
    return statistics.stdev(replicate_preds)


def sign_accuracy(pairs: list[Pair]) -> float | None:
    """Fraction of pairs whose predicted direction matches experiment."""
    scored = [p for p in pairs if p.delta_exp != 0 and p.delta_pred != 0]
    if not scored:
        return None
    hits = sum(1 for p in scored if (p.delta_pred > 0) == (p.delta_exp > 0))
    return hits / len(scored)


def flat_fraction(pairs: list[Pair], threshold: float) -> float | None:
    """Fraction of pairs the model calls effectively identical.

    `threshold` should come from noise_sigma, not be picked by hand.
    """
    if not pairs:
        return None
    return sum(1 for p in pairs if abs(p.delta_pred) <= threshold) / len(pairs)


def summarize(pairs: list[Pair], threshold: float) -> dict:
    usable = [p for p in pairs if p.usable]
    same_geom = [p for p in usable if p.same_realized_geometry]
    return {
        "n_pairs": len(pairs),
        "n_usable": len(usable),
        "n_ambiguous": len(pairs) - len(usable),
        "sign_accuracy": sign_accuracy(usable),
        "flat_fraction": flat_fraction(usable, threshold),
        "threshold": threshold,
        # The key conditional.
        "n_same_realized_geometry": len(same_geom),
        "mean_abs_delta_pred_same_geometry": (
            statistics.fmean(abs(p.delta_pred) for p in same_geom) if same_geom else None
        ),
        "mean_abs_delta_exp": (
            statistics.fmean(abs(p.delta_exp) for p in usable) if usable else None
        ),
        "mean_abs_delta_pred": (
            statistics.fmean(abs(p.delta_pred) for p in usable) if usable else None
        ),
    }
