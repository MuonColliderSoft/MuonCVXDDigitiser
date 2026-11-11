#
# k4run options digitising all the MAIA tracker sub-detectors with MuonCVXDDigitiser, fired pixels stored
#
# k4run runValidationDigi.py --compact MAIA_v0.xml --config {default,nosmear} \
#       --IOSvc.Input sim.edm4hep.root --IOSvc.Output digi.edm4hep.root
#
import os
import sys

from Gaudi.Configuration import INFO
from k4FWCore import ApplicationMgr, IOSvc
from k4FWCore.parseArgs import parser
from Configurables import EventDataSvc, GeoSvc, UniqueIDGenSvc
from Configurables import Gaudi__Monitoring__MessageSvcSink as MessageSvcSink
from Configurables import MuonCVXDDigitiser

# common.py is found through PYTHONPATH (set by run_validation.sh) or next to this file
if "__file__" in globals():
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import CONFIGS, SUBDETECTORS  # noqa: E402

parser.add_argument("--compact", required=True, help="Compact file of the detector geometry")
parser.add_argument("--config", choices=list(CONFIGS), default="default")
args = parser.parse_known_args()[0]

GeoSvc("GeoSvc", detectors=[args.compact], EnableGeant4Geo=False, OutputLevel=INFO)
IOSvc()

digitisers = [
    MuonCVXDDigitiser(
        prefix + "Digitiser",
        SubDetectorName=subdet,
        LayerIDs=layer_ids,
        CollectionName=[collection],
        EventHeader=["EventHeader"],
        OutputCollectionName=[prefix + "Hits"],
        RelationColName=[prefix + "HitsRelations"],
        SimHitLocCollectionName=[prefix + "Pixels"],
        RawHitsLinkColName=[prefix + "RawHitRelations"],
        StoreFiredPixels=True,
        **CONFIGS[args.config],
    )
    for prefix, subdet, collection, layer_ids, _ in SUBDETECTORS
]

ApplicationMgr(
    TopAlg=digitisers,
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc"), UniqueIDGenSvc("UniqueIDGenSvc"), MessageSvcSink()],
    OutputLevel=INFO,
)
