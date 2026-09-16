#!/usr/bin/env python3
# Check that two digitisation outputs contain identical hits, pixels and links
#
# Every TrackerHitPlane, SimTrackerHit and link collection of the files is compared, event by event.
#
# Usage: compare_digi.py a.edm4hep.root b.edm4hep.root
import sys

from podio.root_io import Reader


def row(obj):
    if hasattr(obj, "getFrom"):
        return (obj.getFrom().id().index, obj.getTo().id().index, obj.getWeight())
    pos = obj.getPosition()
    return (obj.getCellID(), pos.x, pos.y, pos.z, obj.getTime(), obj.getEDep())


def dump(fname):
    out = {}
    for ev in Reader(fname).get("events"):
        evn = ev.get("EventHeader")[0].getEventNumber()
        collections = {}
        for name in ev.getAvailableCollections():
            coll = ev.get(name)
            if any(t in str(coll.getValueTypeName()) for t in ("TrackerHitPlane", "SimTrackerHit", "Link")):
                collections[str(name)] = [row(o) for o in coll]
        out[evn] = collections
    return out


a, b = dump(sys.argv[1]), dump(sys.argv[2])
nrows = sum(len(rows) for ev in a.values() for rows in ev.values())
same = a == b
print(f"{sys.argv[1]} vs {sys.argv[2]}: events {len(a)}/{len(b)}, rows {nrows}, identical = {same}")
if not same:
    for evn in sorted(set(a) | set(b)):
        for name in sorted(set(a.get(evn, {})) | set(b.get(evn, {}))):
            if a.get(evn, {}).get(name) != b.get(evn, {}).get(name):
                print(f"  event {evn}: {name} differs")
sys.exit(0 if same else 1)
