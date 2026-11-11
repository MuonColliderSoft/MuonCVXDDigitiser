#
# Digitise the MAIA vertex and tracker SimTrackerHits with MuonCVXDDigitiser
#
# k4run runMuonCVXDDigitiser.py --IOSvc.Input sim.edm4hep.root --IOSvc.Output digi.edm4hep.root [--threads N] [--compact MAIA_v0.xml]
#
import os

from Gaudi.Configuration import INFO, WARNING
from k4FWCore import ApplicationMgr, IOSvc
from k4FWCore.parseArgs import parser
from Configurables import EventDataSvc, GeoSvc, UniqueIDGenSvc
from Configurables import HiveSlimEventLoopMgr, HiveWhiteBoard, AvalancheSchedulerSvc
from Configurables import Gaudi__Monitoring__MessageSvcSink as MessageSvcSink
from Configurables import MuonCVXDDigitiser

parser.add_argument("--threads", type=int, default=1, help="Number of threads (and event slots)")
parser.add_argument("--compact", default=os.environ["K4GEO"] + "/MuColl/MAIA/compact/MAIA_v0/MAIA_v0.xml",
                    help="Compact file of the detector geometry")
args = parser.parse_known_args()[0]

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [args.compact]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

iosvc = IOSvc()
iosvc.Input = "sim.edm4hep.root"
iosvc.Output = "digi.edm4hep.root"

# name, SubDetectorName, input collection, output prefix, LayerIDs
subdetectors = [
    ("VXDBarrelDigitiser", "VertexBarrel", "VertexBarrelCollection", "VXDBarrel", [0, 1, 2, 4, 6]),
    ("VXDEndcapDigitiser", "VertexEndcap", "VertexEndcapCollection", "VXDEndcap", [0, 1, 2, 3, 4, 5, 6, 7]),
    ("ITBarrelDigitiser", "InnerTrackerBarrel", "InnerTrackerBarrelCollection", "ITBarrel", [0, 1, 2]),
    ("ITEndcapDigitiser", "InnerTrackerEndcap", "InnerTrackerEndcapCollection", "ITEndcap", [0, 1, 2, 3, 4, 5, 6]),
    ("OTBarrelDigitiser", "OuterTrackerBarrel", "OuterTrackerBarrelCollection", "OTBarrel", [0, 1, 2]),
    ("OTEndcapDigitiser", "OuterTrackerEndcap", "OuterTrackerEndcapCollection", "OTEndcap", [0, 1, 2, 3]),
]

digitisers = [
    MuonCVXDDigitiser(
        name,
        SubDetectorName=subdet,
        LayerIDs=layer_ids,
        CollectionName=[input_collection],
        EventHeader=["EventHeader"],
        SimHitLocCollectionName=[prefix + "Pixels"],
        OutputCollectionName=[prefix + "Hits"],
        RelationColName=[prefix + "HitsRelations"],
        RawHitsLinkColName=[prefix + "RawHitRelations"],
        OutputLevel=INFO,
    )
    for name, subdet, input_collection, prefix, layer_ids in subdetectors
]

extSvc = [UniqueIDGenSvc("UniqueIDGenSvc"), MessageSvcSink()]
mt = {}
if args.threads > 1:
    extSvc.append(HiveWhiteBoard("EventDataSvc", EventSlots=args.threads, ForceLeaves=True))
    mt = dict(
        EventLoop=HiveSlimEventLoopMgr(SchedulerName="AvalancheSchedulerSvc", OutputLevel=WARNING),
        MessageSvcType="InertMessageSvc",
    )
    AvalancheSchedulerSvc(ThreadPoolSize=args.threads, OutputLevel=WARNING)
else:
    extSvc.append(EventDataSvc("EventDataSvc"))

ApplicationMgr(
    TopAlg=digitisers,
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=extSvc,
    OutputLevel=INFO,
    **mt,
)
