#!/usr/bin/env python3
"""Re-introduce in a copy of the sources the upstream bugs fixed by the port of MuonCVXDDigitiser.

A copy built with all the fixes reverted must agree with the Marlin processor: this shows that the
port is faithful and that the remaining differences come from the fixes only.

Usage: revert_fixes.py <source directory> <fix> [<fix> ...]
Fixes: cutondeltarays, ladderlength, threshold
"""
import sys

FIXES = {
    # The fluctuation model modifies the delta-ray cut it is given; upstream passed the processor
    # member, so the cut changed from segment to segment, hit to hit and event to event.
    "cutondeltarays": [
        ("MuonCVXDDigitiser/components/MuonCVXDDigitiser.h",
         "    MyG4UniversalFluctuationForSi m_fluctuate;",
         "    MyG4UniversalFluctuationForSi m_fluctuate;\n    mutable double m_cutOnDeltaRaysUpstream{-1.};"),
        ("MuonCVXDDigitiser/components/MuonCVXDDigitiser.cpp",
         "        double tmax = m_cutOnDeltaRays;",
         "        if (m_cutOnDeltaRaysUpstream < 0) m_cutOnDeltaRaysUpstream = m_cutOnDeltaRays;\n"
         "        double& tmax = m_cutOnDeltaRaysUpstream;"),
    ],
    # The ladder length of tracker barrels without sensor length was taken from zHalfSensitive in cm.
    "ladderlength": [
        ("MuonCVXDDigitiser/src/LayerGeometry.cc",
         "2 * z_layout.zHalfSensitive * dd4hep::cm / dd4hep::mm;",
         "2 * z_layout.zHalfSensitive;"),
    ],
    # The smeared threshold was accumulated from pixel to pixel.
    "threshold": [
        ("MuonCVXDDigitiser/components/MuonCVXDDigitiser.cpp",
         "    for (auto& [pixelID, hit] : pixels)\n    {\n        double smear = 0;\n        float origCharge = hit.charge;\n"
         "        // Each pixel has its own threshold, smeared around the nominal one\n"
         "        if (m_thresholdSmearSigma > 0) smear = random.gauss(0., m_thresholdSmearSigma);\n"
         "        float actualThreshold = m_threshold + smear;",
         "    float actualThreshold = m_threshold;\n"
         "    for (auto& [pixelID, hit] : pixels)\n    {\n        double smear = 0;\n        float origCharge = hit.charge;\n"
         "        if (m_thresholdSmearSigma > 0) smear = random.gauss(0., m_thresholdSmearSigma);\n"
         "        actualThreshold = actualThreshold + smear;"),
    ],
}

source = sys.argv[1]
for fix in sys.argv[2:]:
    if fix not in FIXES:
        sys.exit(f"Unknown fix '{fix}', choose among {', '.join(FIXES)}")
    for path, old, new in FIXES[fix]:
        full_path = f"{source}/{path}"
        text = open(full_path).read()
        if text.count(old) != 1:
            sys.exit(f"Cannot revert '{fix}': the code to patch was not found exactly once in {path}. "
                     "Update revert_fixes.py to the current sources.")
        open(full_path, "w").write(text.replace(old, new))
    print(f"Reverted fix '{fix}'")
