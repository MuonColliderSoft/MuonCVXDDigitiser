// Checks of the framework-independent helpers in MuonCVXDDigitiserCore:
// geometry loading, cell ID coding and reproducibility of the fluctuation models.
//
// Usage: test_shared_infrastructure <compact.xml>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "CLHEP/Random/MixMaxRng.h"
#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"

#include "G4UniversalFluctuation.h"
#include "LayerGeometry.h"
#include "MyG4UniversalFluctuationForSi.h"
#include "TrackerCellID.h"

namespace {

int failures = 0;

void check(bool condition, const std::string& what)
{
    if (!condition) {
        std::cerr << "FAILED: " << what << std::endl;
        ++failures;
    }
}

void checkGeometry(const dd4hep::Detector& detector, const TrackerCellID& cellIDCoder,
                   const std::string& name, int zSegmented, bool expectZSegmented)
{
    LayerGeometry geo = loadLayerGeometry(detector, name, zSegmented, {});
    geo.print(std::cout);

    check(geo.numberOfLayers > 0, name + ": has layers");
    check(geo.zSegmented == expectZSegmented, name + ": z-segmentation resolved");
    check(geo.surfaceMap && !geo.surfaceMap->empty(), name + ": has surfaces");

    for (int i = 0; i < geo.numberOfLayers; ++i) {
        const std::string layer = name + " layer " + std::to_string(i);
        check(geo.layerThickness[i] > 0, layer + ": thickness > 0");
        check(std::abs(geo.layerHalfThickness[i] - 0.5 * geo.layerThickness[i]) < 1e-6, layer + ": half thickness");
        if (geo.type.isBarrel) {
            check(geo.laddersInLayer[i] > 0, layer + ": ladders > 0");
            check(geo.layerRadius[i] > 0, layer + ": radius > 0");
            check(geo.layerLadderLength[i] > 0, layer + ": ladder length > 0");
            check(geo.layerLadderWidth[i] > 0, layer + ": ladder width > 0");
            check(geo.pixelsInRow(i, 0.025) > 0 && geo.pixelsInColumn(i, 0.025) > 0, layer + ": pixel counts");
        } else {
            check(geo.layerPetalLength[i] > 0, layer + ": petal length > 0");
            check(geo.layerPetalOuterWidth[i] > 0, layer + ": petal outer width > 0");
        }
    }

    // Every surface of the sub-detector must decode to its system ID and to a layer ID
    // that the default MAIA LayerIDs configuration can map
    std::vector<int> seenLayerIDs;
    for (const auto& [id, surface] : *geo.surfaceMap) {
        check(cellIDCoder.system(id) == geo.detectorID, name + ": surface system ID");
        int layerID = cellIDCoder.layer(id);
        if (std::find(seenLayerIDs.begin(), seenLayerIDs.end(), layerID) == seenLayerIDs.end())
            seenLayerIDs.push_back(layerID);
    }
    std::sort(seenLayerIDs.begin(), seenLayerIDs.end());
    std::cout << name << " surface layer IDs:";
    for (int id : seenLayerIDs) std::cout << " " << id;
    std::cout << "\n\n";
    check((int)seenLayerIDs.size() == geo.numberOfLayers, name + ": one layer ID per geometry layer");

    geo.layerIDs = seenLayerIDs;
    for (int i = 0; i < geo.numberOfLayers; ++i)
        check(geo.layerIndex(seenLayerIDs[i]) == i, name + ": layerIndex mapping");
    check(geo.layerIndex(-999) == -1, name + ": layerIndex of unknown ID");
}

void checkCellID(const TrackerCellID& coder)
{
    for (int side : {-1, 0, 1}) {
        auto cellID = coder.encode(3, side, 5, 1234, 17);
        check(coder.system(cellID) == 3, "cellID system round trip");
        check(coder.side(cellID) == side, "cellID side round trip (side " + std::to_string(side) + ")");
        check(coder.layer(cellID) == 5, "cellID layer round trip");
        check(coder.module(cellID) == 1234, "cellID module round trip");
        check(coder.sensor(cellID) == 17, "cellID sensor round trip");

        auto other = coder.withSensor(cellID, 42);
        check(coder.sensor(other) == 42 && coder.module(other) == 1234 && coder.side(other) == side,
              "cellID withSensor");
    }

    bool threw = false;
    try {
        TrackerCellID bad("system:5,layer:6");
    } catch (const std::exception&) {
        threw = true;
    }
    check(threw, "cellID coder rejects encodings without side/module/sensor");
}

template <typename Model, typename Tmax>
std::vector<double> sample(const Model& model, long seed, Tmax tmaxInit)
{
    CLHEP::MixMaxRng engine(seed);
    std::vector<double> losses;
    // Scan momenta and lengths so that both the Gaussian and the non-Gaussian branches are used
    for (double momentum : {0.5, 5., 100., 10000.}) {
        for (double length : {0.001, 0.005, 0.05}) {
            Tmax tmax = tmaxInit;
            losses.push_back(model.SampleFluctuations(engine, momentum, 105.66, tmax, length, 0.39 * length));
        }
    }
    return losses;
}

void checkFluctuations()
{
    const MyG4UniversalFluctuationForSi forSi;
    const G4UniversalFluctuation universal;

    check(sample(forSi, 1234, 0.03) == sample(forSi, 1234, 0.03), "MyG4UniversalFluctuationForSi reproducible");
    check(sample(forSi, 1234, 0.03) != sample(forSi, 4321, 0.03), "MyG4UniversalFluctuationForSi seed dependent");
    check(sample<G4UniversalFluctuation, const double>(universal, 1234, 0.03) ==
              sample<G4UniversalFluctuation, const double>(universal, 1234, 0.03),
          "G4UniversalFluctuation reproducible");
    check(sample<G4UniversalFluctuation, const double>(universal, 1234, 0.03) !=
              sample<G4UniversalFluctuation, const double>(universal, 4321, 0.03),
          "G4UniversalFluctuation seed dependent");

    for (double loss : sample(forSi, 99, 0.03)) check(std::isfinite(loss) && loss >= 0, "finite non-negative loss");
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <compact.xml>" << std::endl;
        return 2;
    }

    dd4hep::Detector& detector = dd4hep::Detector::getInstance();
    detector.fromCompact(argv[1]);

    const TrackerCellID cellIDCoder(detector.constantAsString("GlobalTrackerReadoutID"));
    std::cout << "Tracker cell ID encoding: " << cellIDCoder.encodingString() << "\n\n";

    checkCellID(cellIDCoder);

    checkGeometry(detector, cellIDCoder, "VertexBarrel", -1, true);
    checkGeometry(detector, cellIDCoder, "VertexEndcap", -1, false);
    checkGeometry(detector, cellIDCoder, "InnerTrackerBarrel", -1, false);
    checkGeometry(detector, cellIDCoder, "InnerTrackerEndcap", -1, false);
    checkGeometry(detector, cellIDCoder, "OuterTrackerBarrel", -1, false);
    checkGeometry(detector, cellIDCoder, "OuterTrackerEndcap", -1, false);
    // Explicit overrides of the tri-state
    check(loadLayerGeometry(detector, "VertexBarrel", 0, {}).zSegmented == false, "ZSegmented=0 on vertex barrel");
    check(loadLayerGeometry(detector, "InnerTrackerBarrel", 1, {}).zSegmented == true, "ZSegmented=1 on tracker barrel");
    check(loadLayerGeometry(detector, "VertexEndcap", 1, {}).zSegmented == false, "ZSegmented=1 ignored on endcap");

    // Tracker barrels have no lengthSensor and fall back to zHalfSensitive, which must be converted to mm
    for (const std::string name : {"InnerTrackerBarrel", "OuterTrackerBarrel"}) {
        LayerGeometry geo = loadLayerGeometry(detector, name, -1, {});
        auto* zPlanarData = detector.detector(name).extension<dd4hep::rec::ZPlanarData>();
        for (int i = 0; i < geo.numberOfLayers; ++i) {
            const double expected = 2 * zPlanarData->layers[i].zHalfSensitive / dd4hep::mm;
            check(std::abs(geo.layerLadderLength[i] - expected) < 1e-3, name + ": ladder length in mm from zHalfSensitive");
        }
    }

    bool threw = false;
    try {
        classifySubDetector("Calorimeter");
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    check(threw, "classifySubDetector rejects names without Barrel/Endcap");

    checkFluctuations();

    if (failures) {
        std::cerr << failures << " check(s) failed" << std::endl;
        return 1;
    }
    std::cout << "All checks passed" << std::endl;
    return 0;
}
