#!/usr/bin/env python3
"""Write a copy of a sample with the deposited energy and the path length of every tracker
SimTrackerHit scaled, everything else unchanged.

Usage: scale_sim_hits.py <input.slcio> <output.slcio> [--edep FACTOR] [--path-length FACTOR]
"""
import argparse
import os
import sys

from pyLCIO import EVENT, IOIMPL

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from digi_common import SUBDETECTORS  # noqa: E402

parser = argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("output")
parser.add_argument("--edep", type=float, default=1.)
parser.add_argument("--path-length", type=float, default=1.)
args = parser.parse_args()

sim_collections = {collection for _, _, collection, _, _ in SUBDETECTORS}
factory = IOIMPL.LCFactory.getInstance()
reader, writer = factory.createLCReader(), factory.createLCWriter()
writer.open(args.output, EVENT.LCIO.WRITE_NEW)

# Run headers first, then the events opened for update so that the hits can be changed
reader.open(args.input)
while (run := reader.readNextRunHeader()):
    writer.writeRunHeader(run)
reader.close()

reader.open(args.input)
nhits = 0
while (event := reader.readNextEvent(EVENT.LCIO.UPDATE)):
    for name in set(event.getCollectionNames()) & sim_collections:
        hits = event.getCollection(name)
        for i in range(hits.getNumberOfElements()):
            hit = hits.getElementAt(i)
            hit.setEDep(hit.getEDep() * args.edep)
            hit.setPathLength(hit.getPathLength() * args.path_length)
            nhits += 1
    writer.writeEvent(event)
reader.close()
writer.close()
if nhits == 0:
    sys.exit("no tracker SimTrackerHits found")
print(f"{nhits} SimTrackerHits scaled: EDep x {args.edep}, path length x {args.path_length}")
