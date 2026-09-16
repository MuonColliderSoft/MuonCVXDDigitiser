#!/usr/bin/env python3
# Check that two runMuonCVXDDigitiser.py outputs contain identical hits and pixels
#
# Usage: compare_digi.py a.edm4hep.root b.edm4hep.root
import sys
from podio.root_io import Reader
def dump(fname):
    out = {}
    for ev in Reader(fname).get("events"):
        evn = ev.get("EventHeader")[0].getEventNumber()
        rows = []
        for p in ["VXDBarrel", "VXDEndcap", "ITBarrel", "ITEndcap", "OTBarrel", "OTEndcap"]:
            for h in ev.get(p + "Hits"):
                pos = h.getPosition()
                rows.append((p, h.getCellID(), pos.x, pos.y, pos.z, h.getTime(), h.getEDep()))
            for h in ev.get(p + "Pixels"):
                pos = h.getPosition()
                rows.append((p + "Pix", h.getCellID(), pos.x, pos.y, h.getTime(), h.getEDep()))
        out[evn] = rows
    return out
a, b = dump(sys.argv[1]), dump(sys.argv[2])
nrows = sum(len(v) for v in a.values())
same = a == b
print(f"{sys.argv[1]} vs {sys.argv[2]}: events {len(a)}/{len(b)}, rows {nrows}, identical = {same}")
if not same:
    for k in a:
        if a.get(k) != b.get(k): print("  differs in event", k, len(a[k]), len(b.get(k, [])))
sys.exit(0 if same else 1)
