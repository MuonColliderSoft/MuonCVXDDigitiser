#!/usr/bin/env python3
# Sanity checks of the MuonCVXDRealDigitiser output of runMuonCVXDRealDigitiser.py
#
# Usage: check_real_digi.py digi_real.edm4hep.root
import argparse
import math
import sys
from collections import Counter

from podio.root_io import Reader

argparser = argparse.ArgumentParser()
argparser.add_argument("file")
argparser.add_argument("--sim", default="VertexBarrelCollection")
argparser.add_argument("--prefix", default="VXDBarrelReal")
args = argparser.parse_args()

# MAIA GlobalTrackerReadoutID: system:5,side:-2,layer:6,module:11,sensor:8
def layer(cellID):
    return (cellID >> 7) & 0x3F

nsim = nreco = nlinks = bad_cellid = 0
res = []
layers_sim, layers_reco = Counter(), Counter()
linked = set()
for ev in Reader(args.file).get("events"):
    sims = ev.get(args.sim)
    hits = ev.get(args.prefix + "Hits")
    rels = ev.get(args.prefix + "HitsRelations")
    nsim += len(sims)
    nreco += len(hits)
    nlinks += len(rels)
    layers_sim.update(layer(h.getCellID()) for h in sims)
    layers_reco.update(layer(h.getCellID()) for h in hits)
    for link in rels:
        hit, sim = link.getFrom(), link.getTo()
        linked.add(hit.id().index)
        if hit.getCellID() != sim.getCellID():
            bad_cellid += 1
        hp, sp = hit.getPosition(), sim.getPosition()
        res.append(math.dist((hp.x, hp.y, hp.z), (sp.x, sp.y, sp.z)))

res.sort()
median = res[len(res) // 2] if res else float("nan")
print(f"sim hits {nsim}, reco hits {nreco}, links {nlinks}, "
      f"median |dr| {median:.4f} mm, max |dr| {max(res, default=float('nan')):.4f} mm, cellID mismatches {bad_cellid}")
print("layer IDs sim :", dict(sorted(layers_sim.items())))
print("layer IDs reco:", dict(sorted(layers_reco.items())))

failures = []
if nreco < 0.8 * nsim: failures.append("efficiency below 80%")
if nlinks < nreco: failures.append("reco hits without sim links")
if bad_cellid: failures.append("reco hit and linked sim hit on different sensors")
if not res or res[-1] > 0.1: failures.append("residual above 100 um")
if set(layers_reco) != set(layers_sim): failures.append("some layers produced no hits")
for f in failures: print("FAILED:", f)
sys.exit(1 if failures else 0)
