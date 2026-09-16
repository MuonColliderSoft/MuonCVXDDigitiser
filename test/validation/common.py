"""Definitions shared by the validation scripts."""

# prefix used for the output collections, SubDetectorName, SimTrackerHit collection,
# LayerIDs for MAIA_v0, name of the fired pixel collection of the Marlin processor
SUBDETECTORS = [
    ("VXDBarrel", "VertexBarrel", "VertexBarrelCollection", [0, 1, 2, 4, 6], "VBPixels"),
    ("VXDEndcap", "VertexEndcap", "VertexEndcapCollection", [0, 1, 2, 3, 4, 5, 6, 7], "VEPixels"),
    ("ITBarrel", "InnerTrackerBarrel", "InnerTrackerBarrelCollection", [0, 1, 2], "IBPixels"),
    ("ITEndcap", "InnerTrackerEndcap", "InnerTrackerEndcapCollection", [0, 1, 2, 3, 4, 5, 6], "IEPixels"),
    ("OTBarrel", "OuterTrackerBarrel", "OuterTrackerBarrelCollection", [0, 1, 2], "OBPixels"),
    ("OTEndcap", "OuterTrackerEndcap", "OuterTrackerEndcapCollection", [0, 1, 2, 3], "OEPixels"),
]

# Settings of the two digitisation configurations, as MuonCVXDDigitiser properties
CONFIGS = {
    "default": {},
    "nosmear": {"PoissonSmearing": False, "ElectronicEffects": False, "ThresholdSmearSigma": 0, "TimeSmearingModel": 0},
}

# Implementations, in the order they are shown
IMPLEMENTATIONS = {
    "marlin": "Marlin (spg-berkeleylab master)",
    "port": "Gaudi port",
    "revert_cutondeltarays": "Gaudi port, CutOnDeltaRays fix reverted",
    "revert_all": "Gaudi port, all fixes reverted",
}

HIT_COLUMNS = ["subdet", "event", "cellID", "sim_x", "sim_y", "sim_z", "sim_t", "sim_edep", "px", "py", "pz",
               "pathLength", "x", "y", "z", "t", "edep", "npix", "pixels"]


def layer_of(cellid):
    """Layer field of the MAIA GlobalTrackerReadoutID (system:5,side:-2,layer:6,module:11,sensor:8)."""
    return (cellid >> 7) & 0x3F


def pixel_string(pixels):
    """Fired pixels as 'x:y:charge;...', with x and y in pixel units and charge in electrons."""
    return ";".join(f"{x:.1f}:{y:.1f}:{q:.0f}" for x, y, q in pixels)
