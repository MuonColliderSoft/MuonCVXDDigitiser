#!/usr/bin/env python3
"""Dump the hits of the Gaudi MuonCVXDDigitiser, with their SimTrackerHit and fired pixels, to CSV.

With --sim, also dump the SimTrackerHits of the input, to compute efficiencies.

Usage: dump_gaudi.py <digi.edm4hep.root> <hits.csv> [--sim sim_hits.csv]
"""
import argparse
import csv
import os
import sys
from collections import defaultdict

from podio.root_io import Reader

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import HIT_COLUMNS, SUBDETECTORS, pixel_string  # noqa: E402

parser = argparse.ArgumentParser()
parser.add_argument("input")
parser.add_argument("output")
parser.add_argument("--sim", help="Also write the SimTrackerHits to this file")
args = parser.parse_args()

sim_file = open(args.sim, "w", newline="") if args.sim else None
sim_writer = csv.writer(sim_file) if sim_file else None
if sim_writer:
    sim_writer.writerow(["subdet", "event", "cellID", "has_momentum"])

with open(args.output, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(HIT_COLUMNS)
    for frame in Reader(args.input).get("events"):
        event = frame.get("EventHeader")[0].getEventNumber()
        for prefix, _, sim_collection, _, _ in SUBDETECTORS:
            if sim_writer:
                for sim in frame.get(sim_collection):
                    m = sim.getMomentum()
                    has_momentum = (m.x != 0 or m.y != 0 or m.z != 0) and sim.getPathLength() > 0
                    sim_writer.writerow([prefix, event, sim.getCellID(), int(has_momentum)])
            hits, pixels = frame.get(prefix + "Hits"), frame.get(prefix + "Pixels")
            fired = defaultdict(list)
            for link in frame.get(prefix + "RawHitRelations"):
                fired[link.getFrom().id().index].append(pixels[link.getTo().id().index])
            for link in frame.get(prefix + "HitsRelations"):
                index = link.getFrom().id().index
                reco, sim = hits[index], link.getTo()
                rp, sp, sm = reco.getPosition(), sim.getPosition(), sim.getMomentum()
                writer.writerow([prefix, event, sim.getCellID(), sp.x, sp.y, sp.z, sim.getTime(), sim.getEDep(),
                                 sm.x, sm.y, sm.z, sim.getPathLength(), rp.x, rp.y, rp.z, reco.getTime(),
                                 reco.getEDep(), len(fired[index]),
                                 pixel_string((p.getPosition().x, p.getPosition().y, p.getEDep())
                                              for p in fired[index])])

if sim_file:
    sim_file.close()
