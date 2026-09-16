#include "LayerGeometry.h"

#include <algorithm>
#include <cmath>
#include <ostream>
#include <stdexcept>

#include "DD4hep/DD4hepUnits.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"

using dd4hep::DetElement;
using dd4hep::rec::SurfaceManager;
using dd4hep::rec::ZDiskPetalsData;
using dd4hep::rec::ZPlanarData;

SubDetectorType classifySubDetector(const std::string& subDetName)
{
    SubDetectorType type;

    if (subDetName.find("Barrel") != std::string::npos) {
        type.isBarrel = true;
    } else if (subDetName.find("Endcap") != std::string::npos) {
        type.isBarrel = false;
    } else {
        throw std::invalid_argument("Could not determine sub-detector type for: " + subDetName);
    }

    if (subDetName.find("Vertex") != std::string::npos) {
        type.isVertex = true;
    } else if (subDetName.find("InnerTracker") != std::string::npos) {
        type.isInnerTracker = true;
    } else if (subDetName.find("OuterTracker") != std::string::npos) {
        type.isOuterTracker = true;
    }

    return type;
}

SensorExtent surfaceExtent(dd4hep::rec::ISurface& surface)
{
    SensorExtent extent;
    extent.halfLengthU = 0.5 * surface.length_along_u() / dd4hep::mm;
    extent.halfLengthV = 0.5 * surface.length_along_v() / dd4hep::mm;

    auto* volSurface = dynamic_cast<dd4hep::rec::Surface*>(&surface);
    if (!volSurface) return extent;
    const auto lines = volSurface->getLines();
    if (lines.empty()) return extent;

    const dd4hep::rec::Vector3D origin = surface.origin();
    const dd4hep::rec::Vector3D u = surface.u();
    const dd4hep::rec::Vector3D v = surface.v();
    double halfU = 0;
    double halfV = 0;
    for (const auto& [start, end] : lines) {
        for (const dd4hep::rec::Vector3D& point : {start, end}) {
            const dd4hep::rec::Vector3D d = point - origin;
            halfU = std::max(halfU, std::abs(d.dot(u)));
            halfV = std::max(halfV, std::abs(d.dot(v)));
        }
    }
    // Rounded to 1 nm, so that rounding errors of the outline do not shift the pixel grid
    extent.halfLengthU = std::round(halfU / dd4hep::mm * 1e6) / 1e6;
    extent.halfLengthV = std::round(halfV / dd4hep::mm * 1e6) / 1e6;
    return extent;
}

bool resolveZSegmented(int zSegmented, const SubDetectorType& type)
{
    return (zSegmented < 0) ? (type.isVertex && type.isBarrel) : (zSegmented > 0 && type.isBarrel);
}

LayerGeometry loadLayerGeometry(const dd4hep::Detector& detector,
                                const std::string& subDetName,
                                int zSegmented,
                                const std::vector<int>& layerIDs)
{
    LayerGeometry geo;
    geo.subDetName = subDetName;
    geo.type = classifySubDetector(subDetName);
    geo.zSegmented = resolveZSegmented(zSegmented, geo.type);
    geo.layerIDs = layerIDs;

    DetElement subDetector = detector.detector(subDetName);
    if (!subDetector.isValid()) {
        throw std::runtime_error("Could not find subdetector: " + subDetName);
    }
    geo.detectorID = subDetector.id();

    std::vector<ZPlanarData::LayerLayout> barrelLayers;
    std::vector<ZDiskPetalsData::LayerLayout> endcapLayers;
    if (geo.type.isBarrel) {
        // Barrel-like geometry
        ZPlanarData* zPlanarData = subDetector.extension<ZPlanarData>(false);
        if (!zPlanarData) {
            throw std::runtime_error("Could not find surface of type ZPlanarData for subdetector: " + subDetName);
        }
        barrelLayers = zPlanarData->layers;
        geo.numberOfLayers = barrelLayers.size();
    } else {
        // Endcap-like geometry
        ZDiskPetalsData* zDiskPetalData = subDetector.extension<ZDiskPetalsData>(false);
        if (!zDiskPetalData) {
            throw std::runtime_error("Could not find surface of type ZDiskPetalsData for subdetector: " + subDetName);
        }
        endcapLayers = zDiskPetalData->layers;
        geo.numberOfLayers = endcapLayers.size();
    }

    SurfaceManager* surfMan = detector.extension<SurfaceManager>(false);
    geo.surfaceMap = surfMan ? surfMan->map(subDetector.name()) : nullptr;
    if (!geo.surfaceMap) {
        throw std::runtime_error("Could not find surface map for detector: " + subDetName + " in SurfaceManager");
    }
    geo.sensorExtents.reserve(geo.surfaceMap->size());
    for (const auto& [cellID, surface] : *geo.surfaceMap) {
        geo.sensorExtents.emplace(cellID, surfaceExtent(*surface));
    }

    const int nLayers = geo.numberOfLayers;
    geo.laddersInLayer.resize(nLayers);
    geo.sensorsPerLadder.resize(nLayers);
    geo.layerHalfPhi.resize(nLayers);
    geo.layerHalfThickness.resize(nLayers);
    geo.layerThickness.resize(nLayers);
    geo.layerRadius.resize(nLayers);
    geo.layerLadderLength.resize(nLayers);
    geo.layerLadderWidth.resize(nLayers);
    geo.layerLadderHalfWidth.resize(nLayers);
    geo.layerActiveSiOffset.resize(nLayers);
    geo.layerPhiOffset.resize(nLayers);
    geo.petalsInLayer.resize(nLayers);
    geo.layerPetalLength.resize(nLayers);
    geo.layerPetalInnerWidth.resize(nLayers);
    geo.layerPetalOuterWidth.resize(nLayers);

    int curr_layer = 0;
    if (geo.type.isBarrel) {
        for (const ZPlanarData::LayerLayout& z_layout : barrelLayers)
        {
            // ALE: Geometry is in cm, convert all length to mm
            geo.laddersInLayer[curr_layer] = z_layout.ladderNumber;
            geo.layerHalfPhi[curr_layer] = M_PI / ((double)geo.laddersInLayer[curr_layer]);
            geo.layerThickness[curr_layer] = z_layout.thicknessSensitive * dd4hep::cm / dd4hep::mm;
            geo.layerHalfThickness[curr_layer] = 0.5 * geo.layerThickness[curr_layer];
            geo.layerRadius[curr_layer] = z_layout.distanceSensitive * dd4hep::cm / dd4hep::mm + geo.layerHalfThickness[curr_layer];
            geo.sensorsPerLadder[curr_layer] = z_layout.sensorsPerLadder;
            if (geo.zSegmented) {
                geo.layerLadderLength[curr_layer] = z_layout.lengthSensor * z_layout.sensorsPerLadder * dd4hep::cm / dd4hep::mm;
            } else {
                geo.layerLadderLength[curr_layer] = z_layout.lengthSensor * dd4hep::cm / dd4hep::mm;
            }
            if ((geo.type.isInnerTracker || geo.type.isOuterTracker) && z_layout.lengthSensor == 0) {
                geo.layerLadderLength[curr_layer] = 2 * z_layout.zHalfSensitive * dd4hep::cm / dd4hep::mm;
            }
            geo.layerLadderWidth[curr_layer] = z_layout.widthSensitive * dd4hep::cm / dd4hep::mm;
            geo.layerLadderHalfWidth[curr_layer] = geo.layerLadderWidth[curr_layer] / 2.;
            geo.layerActiveSiOffset[curr_layer] = - z_layout.offsetSensitive * dd4hep::cm / dd4hep::mm;
            geo.layerPhiOffset[curr_layer] = z_layout.phi0;

            curr_layer++;
        }
    } else {
        for (const ZDiskPetalsData::LayerLayout& z_layout : endcapLayers)
        {
            // see k4geo detector/tracker/VertexEndcap_o1_v06_geo.cpp for how the endcap (VXD) geometry is built.
            // Note: petal-like structure but current geometry only defines a single sensitive element for the whole disk.
            geo.layerThickness[curr_layer] = z_layout.thicknessSensitive * dd4hep::cm / dd4hep::mm;
            geo.layerHalfThickness[curr_layer] = 0.5 * geo.layerThickness[curr_layer];
            geo.layerPetalLength[curr_layer] = z_layout.lengthSensitive * dd4hep::cm / dd4hep::mm;
            geo.layerPetalInnerWidth[curr_layer] = z_layout.widthInnerSensitive * dd4hep::cm / dd4hep::mm;
            geo.layerPetalOuterWidth[curr_layer] = z_layout.widthOuterSensitive * dd4hep::cm / dd4hep::mm;
            geo.petalsInLayer[curr_layer] = z_layout.petalNumber;
            if (geo.layerPetalOuterWidth[curr_layer] == 0) {
                float outerEndcapRadius = 112.0 * dd4hep::cm / dd4hep::mm; // FIX find source
                geo.layerPetalOuterWidth[curr_layer] = 2 * outerEndcapRadius * std::tan(M_PI / geo.petalsInLayer[curr_layer]);
            }

            curr_layer++;
        }
    }

    return geo;
}

int LayerGeometry::layerIndex(int layerID) const
{
    auto it = std::find(layerIDs.begin(), layerIDs.end(), layerID);
    if (it == layerIDs.end()) return -1;
    return std::distance(layerIDs.begin(), it);
}

const SensorExtent* LayerGeometry::sensorExtent(unsigned long cellID) const
{
    auto it = sensorExtents.find(cellID);
    return it == sensorExtents.end() ? nullptr : &it->second;
}

int LayerGeometry::pixelsInColumn(int layer, double pixelSizeX) const
{
    if (type.isBarrel) {
        return std::ceil(layerLadderWidth[layer] / pixelSizeX);
    } else {
        return std::ceil(layerPetalOuterWidth[layer] / pixelSizeX);
    }
}

int LayerGeometry::pixelsInRow(int layer, double pixelSizeY) const
{
    if (type.isBarrel) {
        return std::ceil(layerLadderLength[layer] / pixelSizeY);
    } else {
        return std::ceil(layerPetalLength[layer] / pixelSizeY);
    }
}

void LayerGeometry::print(std::ostream& out) const
{
    out << "Geometry of " << subDetName << " (" << (type.isBarrel ? "barrel" : "endcap")
        << ", z-segmentation " << (zSegmented ? "enabled" : "disabled") << ")\n"
        << "Number of layers: " << numberOfLayers << "\n";
    for (int i = 0; i < numberOfLayers; ++i)
    {
        out << "Layer " << i;
        if (i < (int)layerIDs.size()) out << " (layer ID " << layerIDs[i] << ")";
        out << "\n";
        if (type.isBarrel) {
            out << "  Number of ladders: " << laddersInLayer[i] << "\n"
                << "  Sensors per ladder: " << sensorsPerLadder[i] << "\n"
                << "  Radius: " << layerRadius[i] << "\n"
                << "  Ladder length: " << layerLadderLength[i] << "\n"
                << "  Ladder width: " << layerLadderWidth[i] << "\n"
                << "  Ladder half width: " << layerLadderHalfWidth[i] << "\n"
                << "  Phi offset: " << layerPhiOffset[i] << "\n"
                << "  Active Si offset: " << layerActiveSiOffset[i] << "\n"
                << "  Half phi: " << layerHalfPhi[i] << "\n";
        } else {
            out << "  Number of petals: " << petalsInLayer[i] << "\n"
                << "  Petal length: " << layerPetalLength[i] << "\n"
                << "  Petal inner width: " << layerPetalInnerWidth[i] << "\n"
                << "  Petal outer width: " << layerPetalOuterWidth[i] << "\n";
        }
        out << "  Thickness: " << layerThickness[i] << "\n"
            << "  Half thickness: " << layerHalfThickness[i] << "\n";
    }
}
