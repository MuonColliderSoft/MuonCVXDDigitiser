#
# Tracker digitisation only, with the MAIAConfig configurations of DDPlanarDigi (parametric, default)
# or MuonCVXDDigitiser (--realistic), as digi_steer.py sets them up.
#
# k4run tracker_digi.py --maiaconfig <MAIAConfig>/MAIAConfig [--realistic] [--store-pixels]
#       [--tracker-pixel-v 0.025] [--fix-time-window] [--compact MAIA_v0.xml]
#       --IOSvc.Input sim.edm4hep.root --IOSvc.Output digi.edm4hep.root
#
import argparse
import os
import sys

from Gaudi.Configuration import INFO, WARNING
from k4FWCore import ApplicationMgr, IOSvc
from k4FWCore.parseArgs import parser
from Configurables import EventDataSvc, GeoSvc, UniqueIDGenSvc
from Configurables import Gaudi__Monitoring__MessageSvcSink as MessageSvcSink

parser.add_argument("--maiaconfig", default=os.environ.get("MAIACONFIG", ""),
                    help="Directory of the MAIAConfig steering files (<MAIAConfig>/MAIAConfig)")
parser.add_argument("--compact", default=os.environ.get("K4GEO", "") + "/MuColl/MAIA/compact/MAIA_v0/MAIA_v0.xml",
                    help="Compact file of the detector geometry")
parser.add_argument("--realistic", action="store_true", help="Use the realistic MuonCVXDDigitiser configurations")
parser.add_argument("--store-pixels", action="store_true", help="MuonCVXDDigitiser: store the fired pixels")
parser.add_argument("--tracker-pixel-v", type=float, default=None,
                    help="MuonCVXDDigitiser: override the PixelSizeY (mm) of the inner and outer trackers")
parser.add_argument("--fix-time-window", action="store_true",
                    help="DDPlanarDigi: also correct the time window cut for the propagation time "
                         "(CorrectTimeWindowForPropagation)")
args = parser.parse_known_args()[0]

if not args.maiaconfig:
    sys.exit("Set --maiaconfig or MAIACONFIG to the MAIAConfig/MAIAConfig directory")
sys.path.insert(0, args.maiaconfig)

if args.realistic:
    from TrackerDigi.realistic_vertex import VXDBarrel_cfg, VXDEndcap_cfg
    from TrackerDigi.realistic_inner import ITBarrel_cfg, ITEndcap_cfg
    from TrackerDigi.realistic_outer import OTBarrel_cfg, OTEndcap_cfg
else:
    from TrackerDigi.tracking_vertex import VXDBarrel_cfg, VXDEndcap_cfg
    from TrackerDigi.tracking_inner import ITBarrel_cfg, ITEndcap_cfg
    from TrackerDigi.tracking_outer import OTBarrel_cfg, OTEndcap_cfg

# The arguments the MAIAConfig configuration functions read
cfg_args = argparse.Namespace(doOverlayFull=False, doOverlayIP=False, doTrackerConing=False,
                              storeFiredPixels=args.store_pixels)
algs = [cfg(cfg_args) for cfg in (VXDBarrel_cfg, VXDEndcap_cfg, ITBarrel_cfg, ITEndcap_cfg, OTBarrel_cfg, OTEndcap_cfg)]

if args.realistic:
    for alg in algs:
        # older MAIAConfig versions do not pass storeFiredPixels to the digitisers
        alg.StoreFiredPixels = args.store_pixels
    if args.tracker_pixel_v is not None:
        for alg in algs[2:]:
            alg.PixelSizeY = args.tracker_pixel_v
elif args.fix_time_window:
    for alg in algs:
        alg.CorrectTimeWindowForPropagation = True

GeoSvc("GeoSvc", detectors=[args.compact], EnableGeant4Geo=False, OutputLevel=WARNING)
IOSvc()
ApplicationMgr(TopAlg=algs, EvtSel="NONE", EvtMax=-1, OutputLevel=INFO,
               ExtSvc=[EventDataSvc("EventDataSvc"), UniqueIDGenSvc("UniqueIDGenSvc", Seed=42), MessageSvcSink()])
