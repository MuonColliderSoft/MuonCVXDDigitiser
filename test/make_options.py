#!/usr/bin/env python3
"""Write a k4run options file running MuonCVXDDigitiser (or MuonCVXDRealDigitiser) on the MAIA tracker
sub-detectors.

Usage: make_options.py <config> <compact.xml> <input.edm4hep.root> <output.edm4hep.root> <options.py>
                       [--algorithm TYPE] [--subdet PREFIX ...] [--set [PREFIX.]NAME=VALUE ...] [--threads N]
"""
import argparse
import ast
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from digi_common import CONFIGS, SUBDETECTORS  # noqa: E402

parser = argparse.ArgumentParser()
parser.add_argument("config", choices=sorted(CONFIGS))
parser.add_argument("compact")
parser.add_argument("input")
parser.add_argument("output")
parser.add_argument("options")
parser.add_argument("--algorithm", default="MuonCVXDDigitiser", help="algorithm type")
parser.add_argument("--subdet", action="append", help="only these sub-detectors (output prefix)")
parser.add_argument("--set", action="append", default=[], metavar="[PREFIX.]NAME=VALUE",
                    help="extra algorithm property, applied to every digitiser, or only to the one "
                         "of sub-detector PREFIX (e.g. VXDBarrel.LayerIDs=0 1 2)")
parser.add_argument("--threads", type=int, default=1, help="number of threads and event slots")
args = parser.parse_args()


def value(text):
    """Property value from the command line: a Python literal, a list for space-separated values
    (e.g. LayerIDs), a string otherwise."""
    items = text.split()
    if len(items) > 1:
        return [value(item) for item in items]
    try:
        return ast.literal_eval(text)
    except (ValueError, SyntaxError):
        return text


extra = dict(s.split("=", 1) for s in args.set)
# PREFIX.NAME applies to one sub-detector only
common = {k: value(v) for k, v in extra.items() if "." not in k}
own = {}
for k, v in extra.items():
    if "." in k:
        prefix, name = k.split(".", 1)
        own.setdefault(prefix, {})[name] = value(v)
subdets = [s for s in SUBDETECTORS if not args.subdet or s[0] in args.subdet]

algorithms = []
for prefix, subdet, collection, layer_ids in subdets:
    properties = {
        "SubDetectorName": subdet,
        "LayerIDs": layer_ids,
        "CollectionName": [collection],
        "EventHeader": ["EventHeader"],
        "OutputCollectionName": [prefix + "Hits"],
        "RelationColName": [prefix + "HitsRelations"],
    }
    if args.algorithm == "MuonCVXDDigitiser":
        properties.update({
            "SimHitLocCollectionName": [prefix + "Pixels"],
            "RawHitsLinkColName": [prefix + "RawHitRelations"],
            "StoreFiredPixels": True,
        })
    properties.update({**CONFIGS[args.config], **common, **own.get(prefix, {})})
    body = "".join(f"    {k}={v!r},\n" for k, v in properties.items())
    algorithms.append(f"{args.algorithm}(\n    {prefix + 'Digitiser'!r},\n{body})")

if args.threads > 1:
    event_loop = f"""
extSvc.append(HiveWhiteBoard("EventDataSvc", EventSlots={args.threads}, ForceLeaves=True))
AvalancheSchedulerSvc(ThreadPoolSize={args.threads}, OutputLevel=WARNING)
mt = dict(EventLoop=HiveSlimEventLoopMgr(SchedulerName="AvalancheSchedulerSvc", OutputLevel=WARNING),
          MessageSvcType="InertMessageSvc")
"""
else:
    event_loop = """
extSvc.append(EventDataSvc("EventDataSvc"))
mt = {}
"""

with open(args.options, "w") as f:
    f.write(f"""# Written by make_options.py ({args.config})
from Gaudi.Configuration import INFO, WARNING
from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EventDataSvc, GeoSvc, UniqueIDGenSvc
from Configurables import HiveSlimEventLoopMgr, HiveWhiteBoard, AvalancheSchedulerSvc
from Configurables import {args.algorithm}

GeoSvc("GeoSvc", detectors=[{args.compact!r}], EnableGeant4Geo=False, OutputLevel=WARNING)
iosvc = IOSvc()
iosvc.Input = {args.input!r}
iosvc.Output = {args.output!r}

algorithms = [
{(","+chr(10)).join(algorithms)},
]

extSvc = [UniqueIDGenSvc("UniqueIDGenSvc")]
{event_loop}
ApplicationMgr(TopAlg=algorithms, EvtSel="NONE", EvtMax=-1, ExtSvc=extSvc, OutputLevel=INFO, **mt)
""")
