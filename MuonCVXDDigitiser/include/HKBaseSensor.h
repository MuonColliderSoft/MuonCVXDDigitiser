#ifndef HKBaseSensor_h
#define HKBaseSensor_h 1

#include "PixelDigiMatrix.h"
#include "FindUnionAlgorithm.h"
#include <tuple>
#include <unordered_map>

using std::vector;
using std::tuple;
using std::tie;
using std::unordered_map;

tuple<int, int, int, int> GetBound(const ClusterOfPixel& cluster, GridPosition locate);

/* ****************************************************************************

    Cluster Heap

   ************************************************************************* */
struct ChargePoint
{
    int row;
    int col;
    float charge;
};

struct BufferedCluster
{
    vector<ChargePoint> pixels;
    float time;
};

struct ClusterItem
{
    BufferedCluster buffer;
    int size;
};

using ClusterTable = unordered_map<int, ClusterItem>;

using ReferenceTable = unordered_map<LinearPosition, int>;

class ClusterHeap
{
public:
    ClusterHeap(int rows, int cols);
    virtual ~ClusterHeap();
    void AddCluster(ClusterOfPixel& cluster);
    void SetupPixel(int pos_x, int pos_y, PixelData pix);
    vector<BufferedCluster> PopClusters();
    void SetLabel(string dlabel) { debug_label = dlabel; }

private:
    int hash_cnt;
    GridPosition locate;
    string debug_label;
    ClusterTable cluster_table;
    ReferenceTable ref_table;
    vector<int>  ready_to_pop;
};

/* ****************************************************************************

    Hoshen-Kopelman sensor

   ************************************************************************* */

class HKBaseSensor : public PixelDigiMatrix
{
public:
    HKBaseSensor(int layer,
                          int ladder,
                          int xsegmentNumber,
                          int ysegmentNumber,
                          float ladderLength,
                          float ladderWidth,
                          float thickness,
                          double pixelSizeX,
                          double pixelSizeY,
                          string enc_str,
                          int barrel_id,
                          double thr,
                          float fe_slope,
                          float starttime,
                          float t_step,
                          bool hk8_on = true);
    virtual ~HKBaseSensor() {}

    void buildHits(SegmentDigiHitList& output) override;

protected:
    vector<ClusterHeap> heap_table;
    bool HK8_enabled;
};

#endif //HKBaseSensor_h

