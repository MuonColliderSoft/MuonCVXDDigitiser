#!/usr/bin/env python3
"""Dump the tracker hits of a MAIAConfig digitisation output, with their SimTrackerHit and,
if they were stored, their fired pixels, to CSV.

Every TrackerHitPlane field is written, so that the parametric and realistic digitisers can be
compared field by field. With a third argument, the SimTrackerHits are also written, to compute
efficiencies.

Usage: dump_hits.py <digi.edm4hep.root> <hits.csv> [<sim_hits.csv>]
"""
import csv
import sys

from podio.root_io import Reader

# collection prefix: SimTrackerHit collection
SUBDETS = {"VXDBarrel": "VertexBarrelCollection", "VXDEndcap": "VertexEndcapCollection",
           "ITBarrel": "InnerTrackerBarrelCollection", "ITEndcap": "InnerTrackerEndcapCollection",
           "OTBarrel": "OuterTrackerBarrelCollection", "OTEndcap": "OuterTrackerEndcapCollection"}

input_file, output_file = sys.argv[1:3]
sim_file = sys.argv[3] if len(sys.argv) > 3 else None

fs = open(sim_file, "w", newline="") if sim_file else None
ws = csv.writer(fs) if fs else None
if ws:
    ws.writerow(["subdet", "event", "cellID", "x", "y", "z", "t", "edep", "px", "py", "pz", "pathLength"])

with open(output_file, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["subdet", "event", "cellID", "sim_x", "sim_y", "sim_z", "sim_t", "sim_edep", "px", "py", "pz",
                "pathLength", "reco_cellID", "x", "y", "z", "t", "edep", "edepError", "du", "dv",
                "u_theta", "u_phi", "v_theta", "v_phi", "cov0", "cov1", "cov2", "cov3", "cov4", "cov5",
                "type", "quality", "nlinks_sim", "pixels"])
    for frame in Reader(input_file).get("events"):
        ev = frame.get("EventHeader")[0].getEventNumber()
        names = set(str(n) for n in frame.getAvailableCollections())
        for prefix, sim_collection in SUBDETS.items():
            if ws:
                for s in frame.get(sim_collection):
                    p, m = s.getPosition(), s.getMomentum()
                    ws.writerow([prefix, ev, s.getCellID(), p.x, p.y, p.z, s.getTime(), s.getEDep(),
                                 m.x, m.y, m.z, s.getPathLength()])
            if prefix + "HitsRelations" not in names:
                continue
            hits = frame.get(prefix + "Hits")
            links = frame.get(prefix + "HitsRelations")

            # Fired pixels, if stored, as "u:v:charge:time;..." with u, v in pixel units and charge in electrons.
            # MAIAConfig names them <prefix>Pixels; older versions named them after the subdetector.
            pixels = {}
            pixel_names = [prefix + "Pixels", sim_collection.replace("Collection", "")]
            pixel_name = next((n for n in pixel_names if n in names), None)
            if pixel_name and prefix + "RawHitRelations" in names:
                pixel_collection = frame.get(pixel_name)
                for raw in frame.get(prefix + "RawHitRelations"):
                    px = pixel_collection[raw.getTo().id().index]
                    pp = px.getPosition()
                    pixels.setdefault(raw.getFrom().id().index, []).append(
                        f"{pp.x:.1f}:{pp.y:.1f}:{px.getEDep():.0f}:{px.getTime():.4f}")

            nlinks = {}
            for link in links:
                nlinks[link.getTo().id().index] = nlinks.get(link.getTo().id().index, 0) + 1
            for link in links:
                # the link points to the TrackerHit interface: read the TrackerHitPlane from its collection
                r, s = hits[link.getFrom().id().index], link.getTo()
                rp, sp, sm = r.getPosition(), s.getPosition(), s.getMomentum()
                cov = r.getCovMatrix()
                u, v = r.getU(), r.getV()
                w.writerow([prefix, ev, s.getCellID(), sp.x, sp.y, sp.z, s.getTime(), s.getEDep(), sm.x, sm.y, sm.z,
                            s.getPathLength(), r.getCellID(), rp.x, rp.y, rp.z, r.getTime(), r.getEDep(),
                            r.getEDepError(), r.getDu(), r.getDv(), u.a, u.b, v.a, v.b,
                            *[cov[i] for i in range(6)], r.getType(), r.getQuality(),
                            nlinks[s.id().index], ";".join(pixels.get(link.getFrom().id().index, []))])
if fs:
    fs.close()
