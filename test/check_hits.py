#!/usr/bin/env python3
"""Check a hit dump from dump_hits.py: sanity requirements, then agreement with a reference summary.

Usage: check_hits.py <hits.csv> [--reference ref.json] [--write-reference ref.json]
                     [--layers PREFIX=ID,ID,... ...]

The reference holds the summary statistics of each sub-detector. A change that is meant to alter
the output regenerates it with --write-reference, so the change shows up in the diff of the PR.
"""
import argparse
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from digi_common import layer_of  # noqa: E402
from hit_stats import load_hits, summarise  # noqa: E402

# Sanity requirements on every configuration
# reconstructed hits / SimTrackerHits, with a margin below the values of master (the outer
# tracker endcap loses more hits, e.g. on layers outside LayerIDs)
MIN_EFFICIENCY = {"VXDBarrel": 0.95, "VXDEndcap": 0.95, "ITBarrel": 0.95, "ITEndcap": 0.95,
                  "OTBarrel": 0.9, "OTEndcap": 0.7}
MAX_FAR_FRACTION = 0.01  # reconstructed hits more than 1 mm from their SimTrackerHit
MAX_MEDIAN_RESIDUAL_UM = {"VXD": 15., "IT": 30., "OT": 30.}

# Agreement with the reference: relative tolerance, absolute tolerance. Small enough to catch any
# intended change of the output, large enough for floating-point differences between platforms.
TOLERANCES = {
    "hits": (0.01, 2),
    "efficiency": (0.01, 0.),
    "mean_pixels": (0.02, 0.),
    "mean_charge_keV": (0.02, 0.),
    "median_residual_um": (0.05, 0.2),
    "time_rms_ps": (0.05, 0.5),
}

parser = argparse.ArgumentParser()
parser.add_argument("hits")
parser.add_argument("--reference")
parser.add_argument("--write-reference")
parser.add_argument("--layers", action="append", default=[], metavar="PREFIX=ID,ID,...",
                    help="the only layers (cell-ID layer field) expected to have hits in sub-detector "
                         "PREFIX: hits on other layers fail, and the efficiency counts these layers only")
args = parser.parse_args()

hits, meta = load_hits(args.hits)
layers = {p: {int(i) for i in ids.split(",")} for p, ids in (s.split("=", 1) for s in args.layers)}
failures = []
summary = {}
for subdet, rows in sorted(hits.items()):
    nsim = meta["nsim"].get(subdet, 0)
    if subdet in layers:
        by_layer = meta.get("nsim_by_layer", {}).get(subdet, {})
        nsim = sum(n for layer, n in by_layer.items() if int(layer) in layers[subdet])
        stray = sorted({layer_of(int(r["cellID"]) & 0xffffffff) for r in rows} - layers[subdet])
        if stray:
            failures.append(f"{subdet}: hits on layers {stray}, expected only on {sorted(layers[subdet])}")
    summary[subdet] = summarise(rows, nsim)


def fail(message):
    failures.append(message)


if not summary:
    fail("no reconstructed hits at all")
for subdet, s in summary.items():
    print(f"{subdet:10} " + "  ".join(f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}" for k, v in s.items()))
    if s["non_finite"]:
        fail(f"{subdet}: {s['non_finite']} hits with non-finite position, time or charge")
    if s["empty_clusters"]:
        fail(f"{subdet}: {s['empty_clusters']} hits without fired pixels")
    if s["non_positive_charge"]:
        fail(f"{subdet}: {s['non_positive_charge']} hits with charge <= 0")
    if not s["efficiency"] >= MIN_EFFICIENCY[subdet]:
        fail(f"{subdet}: efficiency {s['efficiency']:.3f} < {MIN_EFFICIENCY[subdet]}")
    if s["far_from_sim"] > MAX_FAR_FRACTION * s["hits"]:
        fail(f"{subdet}: {s['far_from_sim']} hits more than 1 mm from their SimTrackerHit")
    limit = MAX_MEDIAN_RESIDUAL_UM[subdet[:3] if subdet.startswith("VXD") else subdet[:2]]
    if not s["median_residual_um"] <= limit:
        fail(f"{subdet}: median residual {s['median_residual_um']:.1f} um > {limit} um")

if args.reference and not os.path.exists(args.reference):
    fail(f"reference {args.reference} missing: run ctest, then "
         "'cmake --build <build> --target update_test_references'")
elif args.reference:
    with open(args.reference) as f:
        reference = json.load(f)
    if sorted(reference["subdetectors"]) != sorted(summary):
        fail(f"sub-detectors {sorted(summary)} differ from the reference {sorted(reference['subdetectors'])}")
    for subdet, ref in reference["subdetectors"].items():
        if subdet not in summary:
            continue
        for key, (rel, abs_) in TOLERANCES.items():
            new, old = summary[subdet][key], ref[key]
            if math.isnan(old) and math.isnan(new):
                continue
            if not abs(new - old) <= rel * abs(old) + abs_:
                fail(f"{subdet}: {key} = {new:.4g}, reference {old:.4g} (tolerance {rel:.0%} + {abs_})")

if args.write_reference:
    os.makedirs(os.path.dirname(os.path.abspath(args.write_reference)), exist_ok=True)
    with open(args.write_reference, "w") as f:
        json.dump({"events": meta["events"], "subdetectors": summary}, f, indent=1, sort_keys=True)
        f.write("\n")
    print(f"Reference written to {args.write_reference}")

for message in failures:
    print("FAILED:", message)
sys.exit(1 if failures else 0)
