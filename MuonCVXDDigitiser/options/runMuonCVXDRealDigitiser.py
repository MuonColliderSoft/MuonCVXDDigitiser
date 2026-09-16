#
# Digitise the MAIA vertex barrel SimTrackerHits with MuonCVXDRealDigitiser
#
# k4run runMuonCVXDRealDigitiser.py --IOSvc.Input sim.edm4hep.root --IOSvc.Output digi_real.edm4hep.root \
#       [--sensor-type {0,1}] [--stats-file stats.root] [--threads N]
#
import os

from Gaudi.Configuration import INFO, WARNING
from k4FWCore import ApplicationMgr, IOSvc
from k4FWCore.parseArgs import parser
from Configurables import EventDataSvc, GeoSvc, UniqueIDGenSvc
from Configurables import HiveSlimEventLoopMgr, HiveWhiteBoard, AvalancheSchedulerSvc
from Configurables import Gaudi__Monitoring__MessageSvcSink as MessageSvcSink
from Configurables import Gaudi__Histograming__Sink__Root as RootHistoSink
from Configurables import MuonCVXDRealDigitiser

parser.add_argument("--threads", type=int, default=1, help="Number of threads (and event slots)")
parser.add_argument("--sensor-type", type=int, default=1, choices=[0, 1],
                    help="Sensor model: 0 = chip RD53A, 1 = trivial")
parser.add_argument("--stats-file", default="", help="Fill the cluster statistics histograms and write them to this file")
args = parser.parse_known_args()[0]

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [os.environ["K4GEO"] + "/MuColl/MAIA/compact/MAIA_v0/MAIA_v0.xml"]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

iosvc = IOSvc()
iosvc.Input = "sim.edm4hep.root"
iosvc.Output = "digi_real.edm4hep.root"

digitiser = MuonCVXDRealDigitiser(
    "VXDBarrelRealDigitiser",
    SubDetectorName="VertexBarrel",
    LayerIDs=[0, 1, 2, 4, 6],
    SensorType=args.sensor_type,
    CreateStats=bool(args.stats_file),
    CollectionName=["VertexBarrelCollection"],
    EventHeader=["EventHeader"],
    OutputCollectionName=["VXDBarrelRealHits"],
    RelationColName=["VXDBarrelRealHitsRelations"],
    OutputLevel=INFO,
)

extSvc = [UniqueIDGenSvc("UniqueIDGenSvc"), MessageSvcSink()]
if args.stats_file:
    extSvc.append(RootHistoSink("RootHistoSink", FileName=args.stats_file))

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
    TopAlg=[digitiser],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=extSvc,
    OutputLevel=INFO,
    **mt,
)
