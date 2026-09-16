#!/usr/bin/env python3
"""Write the Marlin steering file digitising all the MAIA tracker sub-detectors with MuonCVXDDigitiser.

Usage: make_marlin_steering.py <config> <compact.xml> <input.slcio> <output.slcio> <steering.xml>
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import CONFIGS, SUBDETECTORS  # noqa: E402

config, compact, input_file, output_file, steering = sys.argv[1:6]


def marlin_value(value):
    return str(int(value)) if isinstance(value, bool) else str(value)


processors = []
for prefix, subdet, collection, layer_ids, _ in SUBDETECTORS:
    params = {
        "CollectionName": collection,
        "OutputCollectionName": prefix + "Hits",
        "RelationColName": prefix + "HitsRelations",
        "SubDetectorName": subdet,
        "LayerIDs": " ".join(map(str, layer_ids)),
        "StoreFiredPixels": 1,
        **CONFIGS[config],
    }
    body = "\n".join(f'    <parameter name="{k}">{marlin_value(v)}</parameter>' for k, v in params.items())
    processors.append(f'  <processor name="{prefix}Digitiser" type="MuonCVXDDigitiser">\n{body}\n  </processor>')

execute = "\n".join(f'    <processor name="{prefix}Digitiser"/>' for prefix, *_ in SUBDETECTORS)
newline = "\n"
with open(steering, "w") as f:
    f.write(f"""<marlin>
  <execute>
    <processor name="InitDD4hep"/>
{execute}
    <processor name="Output"/>
  </execute>
  <global>
    <parameter name="LCIOInputFiles">{input_file}</parameter>
    <parameter name="MaxRecordNumber" value="-1"/>
    <parameter name="Verbosity">MESSAGE</parameter>
  </global>
  <processor name="InitDD4hep" type="InitializeDD4hep">
    <parameter name="DD4hepXMLFile">{compact}</parameter>
  </processor>
{newline.join(processors)}
  <processor name="Output" type="LCIOOutputProcessor">
    <parameter name="LCIOOutputFile">{output_file}</parameter>
    <parameter name="LCIOWriteMode">WRITE_NEW</parameter>
  </processor>
</marlin>
""")
