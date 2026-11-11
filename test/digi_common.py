"""Definitions shared by the test scripts."""

# prefix of the output collections, SubDetectorName, SimTrackerHit collection, LayerIDs for MAIA_v0.
# The outputs of a digitiser are <prefix>Hits, <prefix>HitsRelations (links to the SimTrackerHits),
# and with StoreFiredPixels <prefix>Pixels and <prefix>RawHitRelations (links to the fired pixels).
SUBDETECTORS = [
    ("VXDBarrel", "VertexBarrel", "VertexBarrelCollection", [0, 1, 2, 4, 6]),
    ("VXDEndcap", "VertexEndcap", "VertexEndcapCollection", [0, 1, 2, 3, 4, 5, 6, 7]),
    ("ITBarrel", "InnerTrackerBarrel", "InnerTrackerBarrelCollection", [0, 1, 2]),
    ("ITEndcap", "InnerTrackerEndcap", "InnerTrackerEndcapCollection", [0, 1, 2, 3, 4, 5, 6]),
    ("OTBarrel", "OuterTrackerBarrel", "OuterTrackerBarrelCollection", [0, 1, 2]),
    ("OTEndcap", "OuterTrackerEndcap", "OuterTrackerEndcapCollection", [0, 1, 2, 3]),
]

# Smearing switched off: the remaining randomness is the energy-loss fluctuations
NOSMEAR = {"PoissonSmearing": False, "ElectronicEffects": False, "ThresholdSmearSigma": 0,
           "TimeSmearingSigma": 0, "TimeSmearingModel": 0}

# Algorithm properties of each test configuration. Gaudi rejects properties an algorithm does not
# have: MuonCVXDRealDigitiser only runs configurations without MuonCVXDDigitiser-specific ones.
CONFIGS = {
    "default": {},
    "nosmear": NOSMEAR,
    # Charge digitization: 3-bit table, uniform binning
    "bits3": {**NOSMEAR, "ChargeDigitizeNumBits": 3},
    "uniform": {**NOSMEAR, "ChargeDigitizeBinning": 0},
    # Time smearing: realistic model, derived from the sensor thickness
    "timing_realistic": {"TimeSmearingModel": 2},
    # Realistic model with every term overridden: the timewalk alone is left, 0.1 * TRise = 0.1 ns,
    # which must give the same hits as the constant model with a 0.1 ns sigma
    "timing_overrides": {"TimeSmearingModel": 2, "TRise": 1.0, "SigmaLandau": 0, "SigmaJitter": 0,
                         "SigmaTDC": 0, "SigmaClock": 0},
    "timing_sigma100": {"TimeSmearingModel": 1, "TimeSmearingSigma": 0.1},
    # Ionisation trail re-simulated from the local direction and EnergyLoss, not from Geant4
    "resim": {**NOSMEAR, "ResimulateIonisation": True},
    "ms": {**NOSMEAR, "DoMultipleScattering": True},
    "ms_resim": {**NOSMEAR, "DoMultipleScattering": True, "ResimulateIonisation": True},
}

HIT_COLUMNS = ["subdet", "event", "cellID", "sim_x", "sim_y", "sim_z", "sim_t", "sim_edep", "px", "py", "pz",
               "pathLength", "x", "y", "z", "t", "edep", "npix", "pixels"]


def layer_of(cellid):
    """Layer field of the MAIA tracker cell ID (system:5,side:-2,layer:6,module:11,sensor:8)."""
    return (cellid >> 7) & 0x3F


def pixel_string(pixels):
    """Fired pixels as 'x:y:charge:time;...': x and y in pixel units, charge in electrons, time in ns."""
    return ";".join(f"{x:.1f}:{y:.1f}:{q:.0f}:{t:.9g}" for x, y, q, t in pixels)


def pixel_times(pixels):
    """Times of the fired pixels of a pixel_string."""
    return [float(p.split(":")[3]) for p in pixels.split(";")] if pixels else []
