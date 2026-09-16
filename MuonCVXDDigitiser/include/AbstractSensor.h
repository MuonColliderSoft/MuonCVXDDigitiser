#ifndef AbstractSensor_h
#define AbstractSensor_h 1

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <UTIL/BitField64.h>
#include <UTIL/LCTrackerConf.h>
#include <EVENT/SimTrackerHit.h>

using std::string;
using std::vector;
using UTIL::BitField64;
using lcio::LCTrackerCellID;
using EVENT::SimTrackerHit;

enum class PixelStatus : char {
    on,
    off,
    ready,
    start,
    out_of_bounds,
    geometry_error
};

enum class MatrixStatus : char {
    ok,
    pixel_number_error,
    segment_number_error
};

struct PixelData
{
    float charge;
    float time;
    PixelStatus status;
};

using SimHitSet = std::unordered_set<SimTrackerHit*>;

struct SegmentDigiHit
{
    float x;
    float y;
    float charge;
    float time;
    int cellID0;
    int size;
    SimHitSet sim_hits;
};

using SegmentDigiHitList = vector<SegmentDigiHit>;

struct GridCoordinate
{
    int row;
    int col;
};

static inline bool operator==(GridCoordinate a, GridCoordinate b)
{
    return a.row == b.row && a.col == b.col;
}

using LinearPosition = int;

class GridPosition
{
public:
    GridPosition(int rows, int cols) : b_size(cols) {}
    virtual ~GridPosition() {}
    LinearPosition operator()(int row, int col) { return row * b_size + col; }
    LinearPosition operator()(GridCoordinate gc) { return gc.row * b_size + gc.col;}
    GridCoordinate operator()(LinearPosition pos)
    {
        return { pos / b_size, pos % b_size };
    }
private:
    int b_size;
};

using SimHitTable = std::unordered_multimap<LinearPosition, SimTrackerHit*>;

/**
 * @class AbstractSensor
 * @brief Abstract class for segmented sensor
 *
 * Abstract class for segmented sensor.
 */
class AbstractSensor
{
public:

    AbstractSensor( int layer,
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
                    float starttime,
                    float t_step);

    virtual ~AbstractSensor();

    virtual inline int GetLayer() { return _layer; }
    virtual inline int GetLadder() { return _ladder; }
    virtual inline float GetThickness() { return _thickness; }
    virtual inline float GetHalfThickness() { return _thickness / 2; }
    virtual inline float GetLength() { return _ladderLength; }
    virtual inline float GetHalfLength() { return _ladderLength / 2; }
    virtual inline float GetWidth() { return _ladderWidth; }
    virtual inline float GetHalfWidth() { return _ladderWidth / 2; }
    virtual inline double GetPixelSizeX() { return _pixelSizeX; }
    virtual inline double GetPixelSizeY() { return _pixelSizeY; }
    virtual inline int GetLadderRows() { return l_rows; }
    virtual inline int GetLadderCols() { return l_columns; }
    virtual inline int GetSensorRows() { return s_rows; }
    virtual inline int GetSensorCols() { return s_colums; }
    virtual inline int GetSegNumX() { return x_segnum; }
    virtual inline int GetSegNumY() { return y_segnum; }

    virtual int XToPixelRow(double x);
    virtual int YToPixelCol(double y);

    virtual double PixelRowToX(int ix);
    virtual double PixelColToY(int iy);

    virtual inline string GetCellIDFormatStr() { return cellFmtStr; }

    virtual inline MatrixStatus GetStatus() { return status; }

    virtual void InitHitRegister();

    virtual void RegisterHit(int x, int y, SimTrackerHit* hit);

    virtual void Reset() = 0;

    virtual void BeginClockStep() = 0;

    virtual void UpdatePixel(int x, int y, float chrg) = 0;

    virtual void EndClockStep() = 0;

    virtual PixelData GetPixel(int x, int y) = 0;

    virtual bool IsActive() = 0;

    virtual bool CheckStatus(int x, int y, PixelStatus pstat) = 0;

    virtual void buildHits(SegmentDigiHitList& output) = 0;

protected:

    virtual inline int SensorRowToLadderRow(int seg_x, int pos_x) { return seg_x * s_rows + pos_x; }
    virtual inline int SensorColToLadderCol(int seg_y, int pos_y) { return seg_y * s_colums + pos_y; }
    virtual inline int LadderRowToSensorRow(int pos_x, int seg_x) { return pos_x - seg_x * s_rows; }
    virtual inline int LadderColToSensorCol(int pos_y, int seg_y) { return pos_y - seg_y * s_colums; }

    virtual PixelData getPixel(int seg_x, int seg_y, int pos_x, int pos_y);
    virtual bool checkStatus(int seg_x, int seg_y, int pos_x, int pos_y, PixelStatus pstat);
    virtual BitField64 getBFEncoder();
    virtual void fillInHitRelation(SimHitSet& sset, LinearPosition pos);

    virtual bool check(int x, int y);

    int _barrel_id;
    int _layer;
    int _ladder;
    float _thickness;
    double _pixelSizeX;
    double _pixelSizeY;
    float _ladderLength;
    float _ladderWidth;
    int l_rows;
    int l_columns;
    int s_rows;
    int s_colums;
    int x_segnum;
    int y_segnum;
    string cellFmtStr;
    double _thr_level;
    float init_time;
    float clock_step;
    GridPosition l_locate;
    GridPosition s_locate;
    MatrixStatus status;
    bool reset_simtable_at_once;

private:
    SimHitTable simhit_table;
};


#endif //AbstractSensor_h
