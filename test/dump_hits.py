#!/usr/bin/env python3
"""Dump the reconstructed hits of MuonCVXDDigitiser, with their SimTrackerHit and fired pixels, to CSV.

One row per link of <prefix>HitsRelations, i.e. per reconstructed hit and SimTrackerHit pair. npix is
the number of fired pixels of the hit (<prefix>RawHitRelations), or, without fired pixels, the
number of SimTrackerHits linked to it (MuonCVXDRealDigitiser).

Also writes <hits.csv>.json with the number of events and of SimTrackerHits per sub-detector
and per layer.

Usage: dump_hits.py <digi.edm4hep.root> <hits.csv>
"""
import csv
import json
import os
import sys
from collections import defaultdict

from podio.root_io import Reader

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from digi_common import HIT_COLUMNS, SUBDETECTORS, layer_of, pixel_string  # noqa: E402


def vec3(v):
    return [v.x, v.y, v.z]


input_file, output_file = sys.argv[1:3]
nsim = {prefix: 0 for prefix, *_ in SUBDETECTORS}
nsim_by_layer = {prefix: {} for prefix, *_ in SUBDETECTORS}
nevents = 0
with open(output_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(HIT_COLUMNS)
    for event in Reader(input_file).get("events"):
        nevents += 1
        names = set(str(n) for n in event.getAvailableCollections())
        event_number = event.get("EventHeader")[0].getEventNumber()
        for prefix, _, sim_collection, _ in SUBDETECTORS:
            if prefix + "HitsRelations" not in names:
                continue
            if sim_collection in names:
                for sim in event.get(sim_collection):
                    nsim[prefix] += 1
                    layer = str(layer_of(sim.getCellID()))
                    nsim_by_layer[prefix][layer] = nsim_by_layer[prefix].get(layer, 0) + 1
            relations = event.get(prefix + "HitsRelations")

            # Fired pixels of each reconstructed hit, in the order of their links
            pixels = None
            if prefix + "RawHitRelations" in names:
                pixels = defaultdict(list)
                raw_links = event.get(prefix + "RawHitRelations")
                for link in raw_links:
                    pixels[link.getFrom().id().index].append(link.getTo())
                if prefix + "Pixels" in names and len(raw_links) != len(event.get(prefix + "Pixels")):
                    sys.exit(f"{prefix}: {len(raw_links)} raw hit links but "
                             f"{len(event.get(prefix + 'Pixels'))} fired pixels")
            sims_per_hit = defaultdict(int)
            for relation in relations:
                sims_per_hit[relation.getFrom().id().index] += 1

            for relation in relations:
                reco, sim = relation.getFrom(), relation.getTo()
                fired = pixels.get(reco.id().index, []) if pixels is not None else []
                npix = len(fired) if pixels is not None else sims_per_hit[reco.id().index]
                rp, sp, sm = vec3(reco.getPosition()), vec3(sim.getPosition()), vec3(sim.getMomentum())
                writer.writerow([prefix, event_number, sim.getCellID(), sp[0], sp[1], sp[2], sim.getTime(),
                                 sim.getEDep(), sm[0], sm[1], sm[2], sim.getPathLength(), rp[0], rp[1], rp[2],
                                 reco.getTime(), reco.getEDep(), npix,
                                 pixel_string((p.getPosition().x, p.getPosition().y, p.getEDep(), p.getTime())
                                              for p in fired)])

with open(output_file + ".json", "w") as f:
    json.dump({"events": nevents, "nsim": {k: v for k, v in nsim.items() if v},
               "nsim_by_layer": {k: v for k, v in nsim_by_layer.items() if v}}, f, indent=1)
