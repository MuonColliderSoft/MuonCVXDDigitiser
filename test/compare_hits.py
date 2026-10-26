#!/usr/bin/env python3
"""Compare two hit dumps from dump_hits.py hit by hit and in summary.

Usage: compare_hits.py <a.csv> <b.csv> [--label TEXT] [--markdown FILE] [--expect identical]
                       [--max-rel-diff METRIC=VALUE ...]

Hits are matched by event and SimTrackerHit cell ID and position. With --expect identical the
script fails unless every hit is the same in both dumps; --max-rel-diff bounds the relative
difference of a summary metric (e.g. mean_pixels=0.05) in every sub-detector.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hit_stats import hit_key, load_hits, summarise  # noqa: E402

FIELDS = ("x", "y", "z", "t", "edep", "npix", "pixels")

parser = argparse.ArgumentParser()
parser.add_argument("a")
parser.add_argument("b")
parser.add_argument("--label", default="")
parser.add_argument("--markdown", help="append the comparison as a Markdown table to this file")
parser.add_argument("--expect", choices=["identical"])
parser.add_argument("--max-rel-diff", action="append", default=[], metavar="METRIC=VALUE")
args = parser.parse_args()

hits_a, meta_a = load_hits(args.a)
hits_b, meta_b = load_hits(args.b)
limits = {k: float(v) for k, v in (s.split("=", 1) for s in args.max_rel_diff)}

header = ["sub-detector", "hits A", "hits B", "matched", "identical",
          "pixels A", "pixels B", "charge A (keV)", "charge B (keV)", "residual A (um)", "residual B (um)",
          "time RMS A (ps)", "time RMS B (ps)"]
rows, failures, all_identical = [], [], True
for subdet in sorted(set(hits_a) | set(hits_b)):
    ra = {hit_key(r): r for r in hits_a.get(subdet, [])}
    rb = {hit_key(r): r for r in hits_b.get(subdet, [])}
    common = set(ra) & set(rb)
    identical = sum(all(ra[k][f] == rb[k][f] for f in FIELDS) for k in common)
    if not (identical == len(ra) == len(rb)):
        all_identical = False
    sa = summarise(list(ra.values()), meta_a["nsim"].get(subdet, 0))
    sb = summarise(list(rb.values()), meta_b["nsim"].get(subdet, 0))
    rows.append([subdet, len(ra), len(rb), len(common), identical,
                 f"{sa['mean_pixels']:.3f}", f"{sb['mean_pixels']:.3f}",
                 f"{sa['mean_charge_keV']:.2f}", f"{sb['mean_charge_keV']:.2f}",
                 f"{sa['median_residual_um']:.2f}", f"{sb['median_residual_um']:.2f}",
                 f"{sa['time_rms_ps']:.1f}", f"{sb['time_rms_ps']:.1f}"])
    for metric, limit in limits.items():
        old, new = sa[metric], sb[metric]
        if not abs(new - old) <= limit * abs(old):
            failures.append(f"{subdet}: {metric} {old:.4g} -> {new:.4g}, more than {limit:.0%} apart")

if args.expect == "identical" and not all_identical:
    failures.append("the dumps are not identical")

title = f"{args.label}: {os.path.basename(args.a)} (A) vs {os.path.basename(args.b)} (B)"
print(f"== {title}")
widths = [max(len(str(x)) for x in col) for col in zip(header, *rows)]
for row in [header] + rows:
    print("  ".join(str(x).rjust(w) for x, w in zip(row, widths)))
print("identical" if all_identical else "DIFFERENT")

if args.markdown:
    with open(args.markdown, "a") as f:
        f.write(f"\n#### {title}: {'identical' if all_identical else '**different**'}\n\n")
        f.write("| " + " | ".join(header) + " |\n|" + "---|" * len(header) + "\n")
        for row in rows:
            f.write("| " + " | ".join(str(x) for x in row) + " |\n")

for message in failures:
    print("FAILED:", message)
sys.exit(1 if failures else 0)
