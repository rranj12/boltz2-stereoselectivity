"""Collect a Boltz-2 run into the pilot table and the summary metrics.

Produces the chain the project is after, one row per pair per replicate:

    requested stereo -> reference conformer -> realized pose stereo
    -> predicted affinity -> experimental affinity

Checkpoint 3 (the affinity-stage geometry) is absent by construction: Boltz
never writes it. See the README's finding 2.

The `--synthetic` mode fabricates outputs under an explicit hypothesis so the
whole path can be exercised without a GPU. That is not a toy: it is a power
check. If the analysis cannot separate `conformer` from `smiles` on fabricated
data where the answer is known, it will not separate them on real data either,
and that is worth finding out before spending GPU hours.
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

import pandas as pd

from stereobench.boltz_io import Affinity, parse_affinity
from stereobench.metrics import Pair, noise_sigma, summarize

SYNTHETIC_MODES = ("flat", "conformer", "smiles", "noise_only")


def synthesize(manifest: pd.DataFrame, mode: str, sigma: float, seed: int) -> dict:
    """Fabricate one affinity per run under a stated hypothesis.

    - flat        : the head ignores stereochemistry entirely
    - conformer   : the head follows the reference conformer's handedness
    - smiles      : the head follows the molecular graph
    - noise_only  : no signal at all, only the noise floor
    """
    rng = random.Random(seed)
    out = {}
    for row in manifest.itertuples():
        base = 1.0
        if mode == "conformer":
            base += 0.8 if row.conformer_source == "a" else -0.8
        elif mode == "smiles":
            base += 0.8 if row.graph_source == "a" else -0.8
        elif mode == "flat":
            base += 0.05 * (1 if row.graph_source == "a" else -1)
        out[row.run_name] = base + rng.gauss(0.0, sigma)
    return out


def load_predictions(manifest: pd.DataFrame, out_dir: Path) -> dict:
    """Read affinity_<name>.json for every run that produced one."""
    predictions = {}
    for row in manifest.itertuples():
        path = out_dir / "predictions" / row.run_name / f"affinity_{row.run_name}.json"
        if path.exists():
            predictions[row.run_name] = parse_affinity(path).pic50
    return predictions


def build_pairs(table: pd.DataFrame) -> list[Pair]:
    """One `metrics.Pair` per (pair_id, replicate), eutomer vs distomer.

    Within a replicate the two members are the two NATIVE arms -- that is the
    chemistry comparison. The crossed arms are the mechanistic control and are
    analysed separately by `crossover_effect`.
    """
    pairs = []
    for (pair_id, replicate), group in table.groupby(["pair_id", "replicate"]):
        native = group[~group.is_crossed]
        by_source = {r.graph_source: r for r in native.itertuples()}
        if set(by_source) != {"a", "b"}:
            continue
        a, b = by_source["a"], by_source["b"]
        pairs.append(
            Pair(
                pair_id=f"{pair_id}__r{replicate}",
                exp_eutomer=a.eutomer_pact,
                exp_distomer=a.distomer_pact,
                pred_eutomer=a.pred_pic50,
                pred_distomer=b.pred_pic50,
                stereo_eutomer=getattr(a, "realized_stereo", "ambiguous"),
                stereo_distomer=getattr(b, "realized_stereo", "ambiguous"),
                conformer_eutomer=a.conformer_source,
                conformer_distomer=b.conformer_source,
            )
        )
    return pairs


def crossover_effect(table: pd.DataFrame) -> dict:
    """Which input does the affinity follow -- the graph or the conformer?

    Holds one channel fixed and varies the other, which is the only way to
    attribute a response. Reported as mean |delta| in pAct units, directly
    comparable to the noise floor.
    """
    rows = []
    for pair_id, group in table.groupby("pair_id"):
        by_arm = group.groupby("arm_id").pred_pic50.mean().to_dict()
        if not {"a_native", "a_crossed", "b_native", "b_crossed"} <= set(by_arm):
            continue
        rows.append(
            {
                "pair_id": pair_id,
                # Same graph, conformer flipped.
                "conformer_effect_a": by_arm["a_crossed"] - by_arm["a_native"],
                "conformer_effect_b": by_arm["b_crossed"] - by_arm["b_native"],
                # Same conformer, graph flipped. a_crossed and b_native share
                # coordinates exactly, so this isolates the graph.
                "graph_effect": by_arm["a_crossed"] - by_arm["b_native"],
            }
        )
    if not rows:
        return {"n_pairs": 0}

    effects = pd.DataFrame(rows)
    conformer = pd.concat(
        [effects.conformer_effect_a, effects.conformer_effect_b]
    ).abs()
    return {
        "n_pairs": len(effects),
        "mean_abs_conformer_effect": float(conformer.mean()),
        "mean_abs_graph_effect": float(effects.graph_effect.abs().mean()),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", default="runs/pilot/manifest.csv", type=Path)
    ap.add_argument("--out-dir", default="runs/pilot/out", type=Path)
    ap.add_argument("--table", default="runs/pilot/pilot_table.csv", type=Path)
    ap.add_argument("--synthetic", choices=SYNTHETIC_MODES)
    ap.add_argument("--synthetic-sigma", type=float, default=0.15)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    manifest = pd.read_csv(args.manifest)

    if args.synthetic:
        print(f"SYNTHETIC DATA, hypothesis = {args.synthetic!r}. Not a result.\n")
        predictions = synthesize(
            manifest, args.synthetic, args.synthetic_sigma, args.seed
        )
    else:
        predictions = load_predictions(manifest, args.out_dir)

    table = manifest.copy()
    table["pred_pic50"] = table.run_name.map(predictions)
    missing = table.pred_pic50.isna().sum()
    table = table.dropna(subset=["pred_pic50"])
    print(f"runs with a prediction: {len(table)} (missing {missing})")
    if table.empty:
        print("nothing to analyse")
        return

    args.table.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.table, index=False)

    # The noise floor, measured from replicates of identical inputs. Every
    # threshold below is derived from it rather than chosen.
    sigmas = [
        noise_sigma(list(group.pred_pic50))
        for _, group in table.groupby(["pair_id", "arm_id"])
        if len(group) >= 2
    ]
    sigma = float(pd.Series(sigmas).mean()) if sigmas else 0.0
    threshold = 2 * sigma

    print(f"\nnoise floor sigma : {sigma:.3f} pAct  (from {len(sigmas)} input groups)")
    print(f"flat threshold 2s : {threshold:.3f}")

    pairs = build_pairs(table)
    if pairs:
        print("\n--- paired metrics (native arms) ---")
        for key, value in summarize(pairs, threshold).items():
            print(f"{key:35s} {value}")

    print("\n--- crossover: which input drives affinity? ---")
    effect = crossover_effect(table)
    for key, value in effect.items():
        print(f"{key:35s} {value}")

    if effect.get("n_pairs"):
        conformer_effect = effect["mean_abs_conformer_effect"]
        graph_effect = effect["mean_abs_graph_effect"]
        print("\nreading:")
        if max(conformer_effect, graph_effect) < threshold:
            print("  Neither channel moves affinity beyond the noise floor.")
            print("  Stereochemistry has no route to the head except the pose,")
            print("  and the original key conditional is clean as written.")
        elif conformer_effect > graph_effect:
            print("  Affinity tracks the REFERENCE CONFORMER more than the graph.")
            print("  The key conditional is confounded; conformers must be held")
            print("  fixed in every paired comparison.")
        else:
            print("  Affinity tracks the MOLECULAR GRAPH more than the conformer.")
            print("  Something downstream of ref_pos carries stereochemistry.")

    print(f"\nwrote {args.table}")


if __name__ == "__main__":
    main()
