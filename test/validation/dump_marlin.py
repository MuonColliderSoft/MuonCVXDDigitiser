#!/usr/bin/env python3
"""Dump the hits of the Marlin MuonCVXDDigitiser, with their SimTrackerHit and fired pixels, to CSV.

Usage: dump_marlin.py <digi.slcio> <hits.csv>
"""
import csv
import os
import sys

from pyLCIO import IOIMPL

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import HIT_COLUMNS, SUBDETECTORS, pixel_string  # noqa: E402


def cellid(hit):
    return (hit.getCellID0() & 0xffffffff) | ((hit.getCellID1() & 0xffffffff) << 32)


input_file, output_file = sys.argv[1:3]
reader = IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(input_file)
with open(output_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(HIT_COLUMNS)
    for event in reader:
        names = set(event.getCollectionNames())
        for prefix, _, _, _, pixel_collection in SUBDETECTORS:
            pixels = event.getCollection(pixel_collection) if pixel_collection in names else None
            ipix = 0
            for relation in event.getCollection(prefix + "HitsRelations"):
                reco, sim = relation.getFrom(), relation.getTo()
                # The fired pixels of each hit are stored consecutively, in the order of the hits
                npix = reco.getRawHits().size()
                fired = [pixels.getElementAt(ipix + k) for k in range(npix)]
                ipix += npix
                rp, sp, sm = reco.getPosition(), sim.getPosition(), sim.getMomentum()
                writer.writerow([prefix, event.getEventNumber(), cellid(sim), sp[0], sp[1], sp[2], sim.getTime(),
                                 sim.getEDep(), sm[0], sm[1], sm[2], sim.getPathLength(), rp[0], rp[1], rp[2],
                                 reco.getTime(), reco.getEDep(), npix,
                                 pixel_string((p.getPosition()[0], p.getPosition()[1], p.getEDep()) for p in fired)])
            if pixels is not None and ipix != pixels.getNumberOfElements():
                sys.exit(f"{prefix}: {ipix} raw hits but {pixels.getNumberOfElements()} fired pixels")
reader.close()
