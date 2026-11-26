"""Paired metrics.

Everything is in pAct space (higher = stronger) so predicted and experimental
differences carry the same sign convention. Boltz's pred_value runs the other
way, so convert with Affinity.pic50 before coming in here.

The headline question is a null -- does |predicted delta| collapse to zero --
and a null needs a noise floor or it cannot be read. Boltz-2 exposes no seed
flag, so the floor is measured by rerunning identical inputs: see noise_sigma.

A correction to the original design
-----------------------------------
This module once claimed that a non-zero delta under matched realized geometry
proved the affinity number came from somewhere other than the coordinates.
That inference does not hold, and the reason is in `stereobench.conformer`:
stereochemistry reaches the model through TWO coordinate channels, not one.
Besides the predicted pose there is `ref_pos`, the ETKDG reference conformer
embedded from the input SMILES, which enters the atom encoder as pairwise
displacement vectors through a layer that is not reflection-invariant.

So when Boltz builds the same geometry for both members, their reference
conformers still differ -- they came from different SMILES. A non-zero delta
there has two live explanations and the original conditional cannot separate
them. `Pair` therefore tracks both channels, and the conditional is three-way:

- geometry differs            -> delta is expected; compare against experiment
- geometry matched, conformers differ -> AMBIGUOUS; the conformer may explain it
- both matched               -> a non-zero delta really does come from elsewhere

The third state is only reachable by substituting conformers deliberately (see
`conformer.build_crossover`), which is why the crossover arm is not optional.
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
    # Handedness of the ETKDG reference conformer each member was given, as a
    # tag such as "native" or "reflected". Defaults differ because, left alone,
    # Boltz embeds each member from its own SMILES and the two conformers are
    # mirror images -- the confound in its untreated state.
    conformer_eutomer: str = "native_eutomer"
    conformer_distomer: str = "native_distomer"

    @property
    def delta_exp(self) -> float:
        return self.exp_eutomer - self.exp_distomer

    @property
    def delta_pred(self) -> float:
        return self.pred_eutomer - self.pred_distomer

    @property
    def same_realized_geometry(self) -> bool:
        """True when Boltz built the same stereochemistry for both members.

        Necessary for the key conditional but NOT sufficient on its own. The
        reference conformers can still differ here, so a non-zero delta is not
        yet attributable to the pose. See `pose_channel_isolated`.
        """
        return self.stereo_eutomer == self.stereo_distomer and self.stereo_eutomer in (
            "correct",
            "flipped",
        )

    @property
    def same_reference_conformer(self) -> bool:
        """True when both members were given the same reference handedness.

        False for any ordinary run: Boltz embeds each member from its own
        SMILES, so their conformers are mirror images. Only a deliberate
        substitution makes this true.
        """
        return self.conformer_eutomer == self.conformer_distomer

    @property
    def pose_channel_isolated(self) -> bool:
        """Both stereo channels matched -- the conditional with teeth.

        When the realized geometries agree AND both members saw the same
        reference conformer, every route by which stereochemistry could reach
        the affinity head has been equalised. A geometry-driven head must give
        |delta_pred| ~ 0 here. A non-zero delta is then genuinely unexplained
        by coordinates of any kind.
        """
        return self.same_realized_geometry and self.same_reference_conformer

    @property
    def conditional_state(self) -> str:
        """Which of the three interpretive regimes this pair falls into."""
        if not self.same_realized_geometry:
            return "geometry_differs"
        if not self.same_reference_conformer:
            return "confounded_by_conformer"
        return "pose_channel_isolated"

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


def _mean_abs_delta_pred(pairs: list[Pair]) -> float | None:
    return statistics.fmean(abs(p.delta_pred) for p in pairs) if pairs else None


def summarize(pairs: list[Pair], threshold: float) -> dict:
    usable = [p for p in pairs if p.usable]
    by_state = {
        state: [p for p in usable if p.conditional_state == state]
        for state in (
            "geometry_differs",
            "confounded_by_conformer",
            "pose_channel_isolated",
        )
    }
    isolated = by_state["pose_channel_isolated"]
    return {
        "n_pairs": len(pairs),
        "n_usable": len(usable),
        "n_ambiguous": len(pairs) - len(usable),
        "sign_accuracy": sign_accuracy(usable),
        "flat_fraction": flat_fraction(usable, threshold),
        "threshold": threshold,
        # The three-way conditional. Report all of it: the counts are what let
        # a reader see how much of the set actually supports an inference.
        "n_by_conditional_state": {k: len(v) for k, v in by_state.items()},
        "mean_abs_delta_pred_by_state": {
            k: _mean_abs_delta_pred(v) for k, v in by_state.items()
        },
        # The conditional with teeth, kept at top level because it is the
        # headline. None means no pair reached the isolated state -- which is
        # the default outcome without a conformer-substitution arm, and must
        # not be read as "no effect".
        "n_pose_channel_isolated": len(isolated),
        "mean_abs_delta_pred_isolated": _mean_abs_delta_pred(isolated),
        "mean_abs_delta_exp": (
            statistics.fmean(abs(p.delta_exp) for p in usable) if usable else None
        ),
        "mean_abs_delta_pred": _mean_abs_delta_pred(usable),
    }
