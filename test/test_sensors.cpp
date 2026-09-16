// Checks of the pixel clustering in the sensor models used by MuonCVXDRealDigitiser,
// on synthetic pixel patterns.

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "GaudiKernel/MsgStream.h"

#include "HKBaseSensor.h"
#include "TrackerCellID.h"
#include "TrivialSensor.h"

namespace {

int failures = 0;

void check(bool condition, const std::string& what)
{
    if (!condition) {
        std::cerr << "FAILED: " << what << std::endl;
        ++failures;
    }
}

const TrackerCellID cellIDCoder("system:5,side:-2,layer:6,module:11,sensor:8");
MsgStream logStream(nullptr);

// Ladder of 4 pixel rows (x) and 8 pixel columns (y), 25 um pixels, split into ySegments sensors along y
const float pixelSize = 0.025;
const float ladderWidth = 4 * pixelSize;
const float ladderLength = 8 * pixelSize;
const double threshold = 200.;

struct Pixel { int row; int col; };

/// Cluster sizes found by a TrivialSensor for the given fired pixels, sorted
std::vector<int> clusterSizes(int ySegments, const std::vector<Pixel>& fired, std::vector<int>* sensorIDs = nullptr)
{
    TrivialSensor sensor(3, 7, 1, ySegments, ladderLength, ladderWidth, 0.05, pixelSize, pixelSize,
                         cellIDCoder, 1, threshold, 0., 25., logStream);
    check(sensor.GetStatus() == MatrixStatus::ok, "sensor geometry");

    sensor.BeginClockStep();
    for (const auto& p : fired) sensor.UpdatePixel(p.row, p.col, 1000.);
    sensor.EndClockStep();

    SegmentDigiHitList hits;
    sensor.buildHits(hits);

    std::vector<int> sizes;
    for (const auto& hit : hits) {
        sizes.push_back(hit.size);
        check(cellIDCoder.layer(hit.cellID) == 3 && cellIDCoder.module(hit.cellID) == 7 && cellIDCoder.system(hit.cellID) == 1,
              "cluster cell ID");
        if (sensorIDs) sensorIDs->push_back(cellIDCoder.sensor(hit.cellID));
    }
    std::sort(sizes.begin(), sizes.end());
    if (sensorIDs) std::sort(sensorIDs->begin(), sensorIDs->end());
    return sizes;
}

} // namespace

int main()
{
    // A compact L-shaped cluster
    check(clusterSizes(1, {{2, 1}, {2, 2}, {3, 2}}) == std::vector<int>{3}, "L-shaped cluster");

    // Diagonal neighbours are merged (8-connectivity)
    check(clusterSizes(1, {{1, 1}, {2, 2}}) == std::vector<int>{2}, "diagonal neighbours");

    // The first pixel of a row is not a neighbour of the last pixel of the previous row
    check(clusterSizes(1, {{0, 7}, {1, 0}}) == std::vector<int>({1, 1}), "no wrap-around between rows");
    check(clusterSizes(1, {{1, 0}, {1, 7}}) == std::vector<int>({1, 1}), "no wrap-around on the NE neighbour");

    // Pixels on the two sides of a sensor boundary belong to different sensors. The neighbour of the
    // first column of sensor 1 is in sensor 0 and must not be confused with a pixel of sensor 1.
    std::vector<int> sensorIDs;
    check(clusterSizes(2, {{1, 3}, {1, 4}, {0, 7}}, &sensorIDs) == std::vector<int>({1, 1, 1}),
          "clusters split at sensor boundary");
    check(sensorIDs == std::vector<int>({0, 1, 1}), "sensor IDs of the split clusters");

    // Pixels at the far edge of the ladder
    check(clusterSizes(2, {{3, 7}, {3, 6}, {2, 7}}) == std::vector<int>{3}, "cluster at the ladder corner");

    // A pixel below threshold does not fire
    {
        TrivialSensor sensor(0, 0, 1, 1, ladderLength, ladderWidth, 0.05, pixelSize, pixelSize,
                             cellIDCoder, 1, threshold, 0., 25., logStream);
        sensor.BeginClockStep();
        sensor.UpdatePixel(1, 1, 150.);
        sensor.EndClockStep();
        SegmentDigiHitList hits;
        sensor.buildHits(hits);
        check(hits.empty(), "no cluster below threshold");
    }

    // Out-of-range segment configurations are rejected instead of dividing by zero
    {
        TrivialSensor sensor(0, 0, 1, 0, ladderLength, ladderWidth, 0.05, pixelSize, pixelSize,
                             cellIDCoder, 1, threshold, 0., 25., logStream);
        check(sensor.GetStatus() == MatrixStatus::segment_number_error, "zero sensors per ladder rejected");
    }

    if (failures) {
        std::cerr << failures << " check(s) failed" << std::endl;
        return 1;
    }
    std::cout << "All checks passed" << std::endl;
    return 0;
}
