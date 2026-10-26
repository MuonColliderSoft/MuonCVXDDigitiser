"""Summary statistics of a hit dump from dump_hits.py, shared by the check and compare scripts."""
import csv
import json
import math
from collections import defaultdict


def load_hits(path):
    """Hits of a dump grouped by sub-detector, and the dump's metadata (events, SimTrackerHits)."""
    hits = defaultdict(list)
    with open(path) as f:
        for row in csv.DictReader(f):
            hits[row["subdet"]].append(row)
    with open(path + ".json") as f:
        meta = json.load(f)
    return hits, meta


def hit_key(row):
    return (row["event"], row["cellID"], row["sim_x"], row["sim_y"], row["sim_z"])


def residual(row):
    """Distance between the reconstructed and the simulated hit, in mm."""
    return math.dist([float(row[c]) for c in "xyz"], [float(row["sim_" + c]) for c in "xyz"])


def _median(values):
    s = sorted(values)
    return s[len(s) // 2] if s else float("nan")


def _mean(values):
    return sum(values) / len(values) if values else float("nan")


def _rms(values):
    if not values:
        return float("nan")
    m = _mean(values)
    return math.sqrt(sum((v - m) ** 2 for v in values) / len(values))


def summarise(rows, nsim):
    """Statistics of the hits of one sub-detector."""
    finite = [r for r in rows if all(math.isfinite(float(r[c])) for c in ("x", "y", "z", "t", "edep"))]
    dr = [residual(r) for r in finite]
    return {
        "hits": len(rows),
        "efficiency": len(rows) / nsim if nsim else float("nan"),
        "non_finite": len(rows) - len(finite),
        "empty_clusters": sum(int(r["npix"]) < 1 for r in rows),
        "non_positive_charge": sum(float(r["edep"]) <= 0 for r in finite),
        "far_from_sim": sum(d > 1.0 for d in dr),
        "mean_pixels": _mean([int(r["npix"]) for r in rows]),
        "mean_charge_keV": _mean([float(r["edep"]) * 1e6 for r in finite]),
        "median_residual_um": _median(dr) * 1e3,
        "time_rms_ps": _rms([float(r["t"]) - float(r["sim_t"]) for r in finite]) * 1e3,
    }


def summarise_dump(path):
    hits, meta = load_hits(path)
    return {subdet: summarise(rows, meta["nsim"].get(subdet, 0)) for subdet, rows in sorted(hits.items())}, meta
