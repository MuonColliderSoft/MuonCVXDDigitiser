#include "HKBaseSensor.h"

#include "streamlog/streamlog.h"

#include <math.h>
#include <algorithm>
#include <limits>
#include <sstream>

using std::sort;
using std::stringstream;
using int_limits = std::numeric_limits<int>;

tuple<int, int, int, int> GetBound(const ClusterOfPixel& cluster, GridPosition locate)
{
    int row_min = int_limits::max();
    int row_max = -1;
    int col_min = int_limits::max();
    int col_max = -1;
    for (LinearPosition pos : cluster)
    {
        GridCoordinate item = locate(pos);
        if (item.row < row_min) row_min = item.row;
        if (item.row > row_max) row_max = item.row;
        if (item.col < col_min) col_min = item.col;
        if (item.col > col_max) col_max = item.col;
    }
    return std::make_tuple(row_min, row_max, col_min, col_max);
}

/* ****************************************************************************

    Cluster Heap

   ************************************************************************* */

ClusterHeap::ClusterHeap(int rows, int cols) :
    hash_cnt(0),
    locate(rows, cols),
    debug_label("Undefined"),
    cluster_table(),
    ref_table(),
    ready_to_pop()
{}

ClusterHeap::~ClusterHeap()
{}

void ClusterHeap::AddCluster(ClusterOfPixel& cluster)
{
    for (LinearPosition curr_pos : cluster)
    {
        auto ref_item = ref_table.find(curr_pos);
        if (ref_item == ref_table.end())
        {
            ref_table.emplace(curr_pos, hash_cnt);
        }
        else if (streamlog::out.write<streamlog::ERROR>())
#pragma omp critical
        {
            GridCoordinate gcoord = locate(curr_pos);
            streamlog::out() << "Cluster heap " << debug_label << ": conflict for pixel "
                << gcoord.row << ":" << gcoord.col << std::endl;
        }
    }

    ClusterItem cl_item { {}, cluster.size() };
    cluster_table.emplace(hash_cnt, cl_item);

    hash_cnt++;
    if (hash_cnt == int_limits::max())
    {
        hash_cnt = 0;
    }
}

void ClusterHeap::SetupPixel(int pos_x, int pos_y, PixelData pix)
{
    LinearPosition pos = locate(pos_x, pos_y);
    auto r_item = ref_table.find(pos);
    if (r_item != ref_table.end())
    {
        int cluster_id = r_item->second;
        ClusterItem& c_item = cluster_table[cluster_id];

        ChargePoint c_pix { pos_x, pos_y, pix.charge };
        c_item.buffer.pixels.push_back(c_pix);
        c_item.buffer.time = pix.time;

        if (c_item.buffer.pixels.size() == c_item.size)
        {
            ready_to_pop.push_back(cluster_id);
        }

        ref_table.erase(pos);
    }
    else if (streamlog::out.write<streamlog::ERROR>())
#pragma omp critical
    {
        streamlog::out() << "Cluster heap " << debug_label << ": undefined pixel "
            << pos_x << ":" << pos_y << std::endl;
    }
}

vector<BufferedCluster> ClusterHeap::PopClusters()
{
    vector<BufferedCluster> result;
    for (auto cluster_id : ready_to_pop)
    {
        result.push_back(cluster_table[cluster_id].buffer);
        cluster_table.erase(cluster_id);
    }

    ready_to_pop.clear();
    return result;
}

/* ****************************************************************************

    Hoshen-Kopelman sensor

   ************************************************************************* */

HKBaseSensor::HKBaseSensor(int layer,
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
                            bool hk8_on) :
    PixelDigiMatrix(layer,
                    ladder,
                    xsegmentNumber,
                    ysegmentNumber,
                    ladderLength,
                    ladderWidth,
                    thickness,
                    pixelSizeX,
                    pixelSizeY,
                    enc_str,
                    barrel_id,
                    thr,
                    fe_slope,
                    starttime,
                    t_step),
    heap_table(0, { 0, 0 }),
    HK8_enabled(hk8_on)
{
    if (GetStatus() == MatrixStatus::ok)
    {
        heap_table.resize(GetSegNumX() * GetSegNumY(), { GetSensorRows(), GetSensorCols() });

        for (int h = 0; h < this->GetSegNumX(); h++)
        {
            for (int k = 0; k < this->GetSegNumY(); k++)
            {
                stringstream d_label;
                d_label << "[" << GetLayer() << ":" << GetLadder() << ":" << h << ":" << k << "]";
                heap_table[s_locate(h, k)].SetLabel(d_label.str());
            }
        }
    }

    reset_simtable_at_once = false;
}

void HKBaseSensor::buildHits(SegmentDigiHitList& output)
{
    FindUnionAlgorithm  fu_algo { s_rows, s_colums };
    BitField64 bf_encoder = getBFEncoder();

    if (!IsActive()) return;

    for (int h = 0; h < this->GetSegNumX(); h++)
    {
        for (int k = 0; k < this->GetSegNumY(); k++)
        {

            //Sensor segments ordered row first
            LinearPosition sens_id = s_locate(h, k);
            bf_encoder[LCTrackerCellID::sensor()] = sens_id;

            ClusterHeap& c_heap = heap_table[sens_id];

            if (CheckStatusOnSensor(h, k, PixelStatus::start))
            {
                /* ****************************************************************
                   Hoshen-Kopelman Algorithm:
                   https://www.ocf.berkeley.edu/~fricke/projects/hoshenkopelman/hoshenkopelman.html
                   ************************************************************** */

                fu_algo.init();

                for (int i = 0; i < this->GetSensorRows(); i++)
                {
                    for (int j = 0; j < this->GetSensorCols(); j++)
                    {
                        if (!checkStatus(h, k, i, j, PixelStatus::start))
                        {
                            fu_algo.invalidate(i, j);
                            continue;
                        }

                        bool N_is_on = checkStatus(h, k, i - 1, j, PixelStatus::start);
                        bool W_is_on = checkStatus(h, k, i, j - 1, PixelStatus::start);
                        bool NW_is_on = HK8_enabled ? checkStatus(h, k, i - 1, j - 1, PixelStatus::start) : false;
                        bool NE_is_on = HK8_enabled ? checkStatus(h, k, i - 1, j + 1, PixelStatus::start) : false;

                        if (N_is_on and W_is_on)
                        {
                            fu_algo.merge(i - 1, j, i, j - 1);
                            fu_algo.merge(i, j - 1, i, j);
                        }
                        else if (W_is_on and NE_is_on)
                        {
                            fu_algo.merge(i, j - 1, i - 1, j + 1);
                            fu_algo.merge(i - 1, j + 1, i, j);
                        }
                        else if (NW_is_on and NE_is_on)
                        {
                            fu_algo.merge(i - 1, j - 1, i - 1, j + 1);
                            fu_algo.merge(i - 1, j + 1, i, j);
                        }
                        else if (W_is_on)
                        {
                            fu_algo.merge(i, j - 1, i, j);
                        }
                        else if (N_is_on)
                        {
                            fu_algo.merge(i - 1, j, i, j);
                        }
                        else if (NW_is_on)
                        {
                            fu_algo.merge(i - 1, j - 1, i, j);
                        }
                        else if (NE_is_on)
                        {
                            fu_algo.merge(i - 1, j + 1, i, j);
                        }
                    }
                }

                fu_algo.close();

                /* ****************************************************************
                   Cluster buffering
                   ************************************************************** */

                for (ClusterOfPixel c_item : fu_algo.get_clusters())
                {
                    c_heap.AddCluster(c_item);
                }
            }

            for (auto p_item : GetPixelsFromSensor(h, k, PixelStatus::ready))
            {
                c_heap.SetupPixel(p_item.row, p_item.col, p_item.data);
            }

            for (BufferedCluster c_item : c_heap.PopClusters())
            {
                // Very simple implementation: geometric mean
                SegmentDigiHit digiHit = {
                    0., 0., 0.,
                    c_item.time,
                    bf_encoder.lowWord(),
                    {}
                };

                for (ChargePoint c_point : c_item.pixels)
                {
                    int global_row = SensorRowToLadderRow(h, c_point.row);
                    int global_col = SensorColToLadderCol(k, c_point.col);

                    digiHit.x += PixelRowToX(global_row);
                    digiHit.y += PixelColToY(global_col);

                    digiHit.charge += c_point.charge;

                    fillInHitRelation(digiHit.sim_hits, l_locate(global_row, global_col));
                }

                digiHit.x /= c_item.pixels.size();
                digiHit.y /= c_item.pixels.size();

                output.push_back(std::move(digiHit));
            }
        }
    }
}
