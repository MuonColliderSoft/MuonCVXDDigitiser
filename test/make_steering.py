#!/usr/bin/env python3
"""Write a Marlin steering file running MuonCVXDDigitiser (or MuonCVXDRealDigitiser) on the MAIA tracker
sub-detectors.

Usage: make_steering.py <config> <compact.xml> <input.slcio> <output.slcio> <steering.xml>
                        [--processor TYPE] [--subdet PREFIX ...] [--set NAME=VALUE ...]
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from digi_common import CONFIGS, SUBDETECTORS  # noqa: E402

parser = argparse.ArgumentParser()
parser.add_argument("config", choices=sorted(CONFIGS))
parser.add_argument("compact")
parser.add_argument("input")
parser.add_argument("output")
parser.add_argument("steering")
parser.add_argument("--processor", default="MuonCVXDDigitiser", help="processor type")
parser.add_argument("--subdet", action="append", help="only these sub-detectors (output prefix)")
parser.add_argument("--set", action="append", default=[], metavar="[PREFIX.]NAME=VALUE",
                    help="extra processor parameter, applied to every digitiser, or only to the one "
                         "of sub-detector PREFIX (e.g. VXDBarrel.LayerIDs=0 1 2)")
args = parser.parse_args()

extra = dict(s.split("=", 1) for s in args.set)
# PREFIX.NAME applies to one sub-detector only
common = {k: v for k, v in extra.items() if "." not in k}
own = {}
for k, v in extra.items():
    if "." in k:
        prefix, name = k.split(".", 1)
        own.setdefault(prefix, {})[name] = v
subdets = [s for s in SUBDETECTORS if not args.subdet or s[0] in args.subdet]

processors = []
for prefix, subdet, collection, layer_ids, _ in subdets:
    params = {
        "CollectionName": collection,
        "OutputCollectionName": prefix + "Hits",
        "RelationColName": prefix + "HitsRelations",
        "SubDetectorName": subdet,
        "LayerIDs": " ".join(map(str, layer_ids)),
        "StoreFiredPixels": 1,
        **CONFIGS[args.config],
        **common,
        **own.get(prefix, {}),
    }
    body = "\n".join(f'    <parameter name="{k}">{v}</parameter>' for k, v in params.items())
    processors.append(f'  <processor name="{prefix}Digitiser" type="{args.processor}">\n{body}\n  </processor>')

execute = "\n".join(f'    <processor name="{s[0]}Digitiser"/>' for s in subdets)
nl = "\n"
with open(args.steering, "w") as f:
    f.write(f"""<marlin>
  <execute>
    <processor name="InitDD4hep"/>
{execute}
    <processor name="Output"/>
  </execute>
  <global>
    <parameter name="LCIOInputFiles">{args.input}</parameter>
    <parameter name="MaxRecordNumber" value="-1"/>
    <parameter name="Verbosity">MESSAGE</parameter>
  </global>
  <processor name="InitDD4hep" type="InitializeDD4hep">
    <parameter name="DD4hepXMLFile">{args.compact}</parameter>
  </processor>
{nl.join(processors)}
  <processor name="Output" type="LCIOOutputProcessor">
    <parameter name="LCIOOutputFile">{args.output}</parameter>
    <parameter name="LCIOWriteMode">WRITE_NEW</parameter>
  </processor>
</marlin>
""")
