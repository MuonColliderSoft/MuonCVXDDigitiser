#!/usr/bin/env python3
"""Dump the reconstructed hits of MuonCVXDDigitiser, with their SimTrackerHit and fired pixels, to CSV.

Also writes <hits.csv>.json with the number of events and of SimTrackerHits per sub-detector.

Usage: dump_hits.py <digi.slcio> <hits.csv>
"""
import csv
import json
import os
import sys

from pyLCIO import IOIMPL

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from digi_common import HIT_COLUMNS, SUBDETECTORS, pixel_string  # noqa: E402


def cellid(hit):
    return (hit.getCellID0() & 0xffffffff) | ((hit.getCellID1() & 0xffffffff) << 32)


def vec3(v):
    # Index explicitly: iterating the returned pointer has no length and can run off the array
    return [v[0], v[1], v[2]]


input_file, output_file = sys.argv[1:3]
nsim = {prefix: 0 for prefix, *_ in SUBDETECTORS}
nevents = 0
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(input_file)
with open(output_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(HIT_COLUMNS)
    for event in reader:
        nevents += 1
        names = set(event.getCollectionNames())
        for prefix, _, sim_collection, _, pixel_collection in SUBDETECTORS:
            if prefix + "HitsRelations" not in names:
                continue
            if sim_collection in names:
                nsim[prefix] += event.getCollection(sim_collection).getNumberOfElements()
            pixels = event.getCollection(pixel_collection) if pixel_collection in names else None
            relations = event.getCollection(prefix + "HitsRelations")
            ipix = 0
            for i in range(relations.getNumberOfElements()):
                relation = relations.getElementAt(i)
                reco, sim = relation.getFrom(), relation.getTo()
                # The fired pixels of each hit are stored consecutively, in the order of the hits
                npix = reco.getRawHits().size()
                fired = [pixels.getElementAt(ipix + k) for k in range(npix)] if pixels is not None else []
                ipix += npix
                rp, sp, sm = vec3(reco.getPosition()), vec3(sim.getPosition()), vec3(sim.getMomentum())
                writer.writerow([prefix, event.getEventNumber(), cellid(sim), sp[0], sp[1], sp[2], sim.getTime(),
                                 sim.getEDep(), sm[0], sm[1], sm[2], sim.getPathLength(), rp[0], rp[1], rp[2],
                                 reco.getTime(), reco.getEDep(), npix,
                                 pixel_string((p.getPosition()[0], p.getPosition()[1], p.getEDep()) for p in fired)])
            if pixels is not None and ipix != pixels.getNumberOfElements():
                sys.exit(f"{prefix}: {ipix} raw hits but {pixels.getNumberOfElements()} fired pixels")
reader.close()

with open(output_file + ".json", "w") as f:
    json.dump({"events": nevents, "nsim": {k: v for k, v in nsim.items() if v}}, f, indent=1)
