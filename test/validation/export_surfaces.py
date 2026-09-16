#!/usr/bin/env python3
"""Export the sensitive surfaces of the tracker sub-detectors (cellID -> origin, u, v, normal) to CSV.

Lengths in mm.

Usage: export_surfaces.py <compact.xml> <surfaces.csv>
"""
import csv
import os
import sys

import cppyy
import DDG4  # noqa: F401  (loads the DD4hep dictionaries)
from ROOT import dd4hep

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import SUBDETECTORS  # noqa: E402

compact, output_file = sys.argv[1:3]
detector = dd4hep.Detector.getInstance()
detector.fromCompact(compact)
cppyy.include("DDRec/SurfaceManager.h")
surface_manager = detector.extension[dd4hep.rec.SurfaceManager]()

cm = 10.  # DD4hep lengths are in cm
with open(output_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["subdet", "cellID", "ox", "oy", "oz", "ux", "uy", "uz", "vx", "vy", "vz", "nx", "ny", "nz",
                     "len_u", "len_v"])
    for prefix, subdet, *_ in SUBDETECTORS:
        for item in surface_manager.map(subdet):
            s = item.second
            o, u, v, n = s.origin(), s.u(), s.v(), s.normal()
            writer.writerow([prefix, item.first, cm * o[0], cm * o[1], cm * o[2], u[0], u[1], u[2], v[0], v[1], v[2],
                             n[0], n[1], n[2], cm * s.length_along_u(), cm * s.length_along_v()])
