#!/usr/bin/env python3
# Sanity checks of the MuonCVXDDigitiser output of runMuonCVXDDigitiser.py
#
# Usage: check_digi.py digi.edm4hep.root [--pixels]
import argparse, math, sys
from podio.root_io import Reader

argparser = argparse.ArgumentParser()
argparser.add_argument("file")
argparser.add_argument("--pixels", action="store_true", help="Fired pixels were stored")
args = argparser.parse_args()

prefixes = {"VXDBarrel": "VertexBarrelCollection", "VXDEndcap": "VertexEndcapCollection",
            "ITBarrel": "InnerTrackerBarrelCollection", "ITEndcap": "InnerTrackerEndcapCollection",
            "OTBarrel": "OuterTrackerBarrelCollection", "OTEndcap": "OuterTrackerEndcapCollection"}
stats = {p: dict(sim=0, reco=0, links=0, pixels=0, rawlinks=0, res=[], dt=[], edep=[], bad=0) for p in prefixes}
for ev in Reader(args.file).get("events"):
    for p, simname in prefixes.items():
        s = stats[p]
        sims = ev.get(simname); hits = ev.get(p + "Hits"); rels = ev.get(p + "HitsRelations")
        pix = ev.get(p + "Pixels"); raw = ev.get(p + "RawHitRelations")
        s["sim"] += len(sims); s["reco"] += len(hits); s["links"] += len(rels)
        s["pixels"] += len(pix); s["rawlinks"] += len(raw)
        for l in rels:
            h, sh = l.getFrom(), l.getTo()
            if h.getCellID() != sh.getCellID(): s["bad"] += 1
            hp, sp = h.getPosition(), sh.getPosition()
            s["res"].append(math.dist((hp.x, hp.y, hp.z), (sp.x, sp.y, sp.z)))
            s["dt"].append(h.getTime() - sh.getTime())
            s["edep"].append(h.getEDep() / sh.getEDep() if sh.getEDep() > 0 else float("nan"))
        wsum = {}
        for l in raw:
            key = l.getFrom().id().index
            wsum[key] = wsum.get(key, 0.) + l.getWeight()
            if l.getTo().getCellID() != l.getFrom().getCellID(): s["bad"] += 1
        if any(abs(w - 1) > 1e-4 for w in wsum.values()): s["bad"] += 1

def med(v):
    v = sorted(x for x in v if x == x)
    return v[len(v)//2] if v else float("nan")
print(f"{'':10s} {'sim':>5s} {'reco':>5s} {'links':>5s} {'pixels':>6s} {'rawlnk':>6s} {'med|dr|mm':>9s} {'max|dr|mm':>9s} {'med dt ns':>9s} {'medE/Esim':>9s} bad")
for p, s in stats.items():
    print(f"{p:10s} {s['sim']:5d} {s['reco']:5d} {s['links']:5d} {s['pixels']:6d} {s['rawlinks']:6d} "
          f"{med(s['res']):9.4f} {max(s['res'], default=float('nan')):9.4f} {med(s['dt']):9.4f} {med(s['edep']):9.3f} {s['bad']}")

failures = []
for p, s in stats.items():
    if s["reco"] == 0: failures.append(f"{p}: no hits")
    if s["reco"] < 0.8 * s["sim"]: failures.append(f"{p}: efficiency below 80%")
    if s["links"] != s["reco"]: failures.append(f"{p}: not one sim link per hit")
    if s["bad"]: failures.append(f"{p}: {s['bad']} inconsistent links")
    if max(s["res"], default=0) > 0.1: failures.append(f"{p}: residual above 100 um")
    if args.pixels and (s["pixels"] == 0 or s["rawlinks"] != s["pixels"]): failures.append(f"{p}: pixels or raw links missing")
    if not args.pixels and (s["pixels"] or s["rawlinks"]): failures.append(f"{p}: pixels stored although not requested")
for f in failures: print("FAILED:", f)
sys.exit(1 if failures else 0)
