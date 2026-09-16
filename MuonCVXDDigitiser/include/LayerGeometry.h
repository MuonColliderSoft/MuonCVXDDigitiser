#ifndef LayerGeometry_h
#define LayerGeometry_h 1

#include <iosfwd>
#include <string>
#include <vector>

#include "DD4hep/Detector.h"
#include "DDRec/SurfaceManager.h"

/** Sub-detector type, derived from the SubDetectorName
 *  (e.g. "VertexBarrel", "InnerTrackerEndcap").
 */
struct SubDetectorType
{
    bool isBarrel{false};
    bool isVertex{false};
    bool isInnerTracker{false};
    bool isOuterTracker{false};
};

/** Throws std::invalid_argument if the name contains neither "Barrel" nor "Endcap".
 *  The vertex/inner/outer flags are left false when not recognised: callers that
 *  require one of them must check.
 */
SubDetectorType classifySubDetector(const std::string& subDetName);

/** Resolve the ZSegmented tri-state: < 0 = auto (vertex barrel only), 0 = off, > 0 = on.
 *  Segmentation along z is only meaningful for barrel layers.
 */
bool resolveZSegmented(int zSegmented, const SubDetectorType& type);

/** Per-layer geometry of a planar tracker sub-detector, read from the
 *  DDRec ZPlanarData (barrel) or ZDiskPetalsData (endcap) extensions.
 *  All lengths in mm, angles in rad.
 */
struct LayerGeometry
{
    std::string subDetName;
    SubDetectorType type;
    bool zSegmented{false};
    int detectorID{0};
    const dd4hep::rec::SurfaceMap* surfaceMap{nullptr};

    int numberOfLayers{0};

    // barrel
    std::vector<int>   laddersInLayer{};
    std::vector<int>   sensorsPerLadder{};
    std::vector<float> layerRadius{};
    std::vector<float> layerThickness{};
    std::vector<float> layerHalfThickness{};
    std::vector<float> layerLadderLength{};
    std::vector<float> layerLadderWidth{};
    std::vector<float> layerLadderHalfWidth{};
    std::vector<float> layerPhiOffset{};
    std::vector<float> layerActiveSiOffset{};
    std::vector<float> layerHalfPhi{};

    // endcap
    std::vector<float> layerPetalLength{};
    std::vector<float> petalsInLayer{};
    std::vector<float> layerPetalInnerWidth{};
    std::vector<float> layerPetalOuterWidth{};

    /// Layer IDs as found in the cell ID, in the order of the geometry layers
    std::vector<int> layerIDs{};

    /** Index of the geometry layer for a layer ID from the cell ID,
     *  or -1 if the ID is not in layerIDs.
     */
    int layerIndex(int layerID) const;

    /// Number of pixels across the ladder (barrel) or petal (endcap)
    int pixelsInColumn(int layer, double pixelSizeX) const;
    /// Number of pixels along the ladder (barrel) or petal (endcap)
    int pixelsInRow(int layer, double pixelSizeY) const;

    void print(std::ostream& out) const;
};

/** Load the layer geometry for subDetName.
 *  Throws std::invalid_argument or std::runtime_error if the sub-detector,
 *  its layer extension or its surface map cannot be found.
 */
LayerGeometry loadLayerGeometry(const dd4hep::Detector& detector,
                                const std::string& subDetName,
                                int zSegmented,
                                const std::vector<int>& layerIDs);

#endif //LayerGeometry_h
