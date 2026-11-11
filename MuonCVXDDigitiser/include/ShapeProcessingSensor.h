#ifndef ShapeProcessingSensor_h
#define ShapeProcessingSensor_h 1

#include "HKBaseSensor.h"
#include <vector>

using std::vector;

class ShapeProcessingSensor : public HKBaseSensor
{
public:
    ShapeProcessingSensor(int layer,
                          int ladder,
                          int xsegmentNumber,
                          int ysegmentNumber,
                          float ladderLength,
                          float ladderWidth,
                          float thickness,
                          double pixelSizeX,
                          double pixelSizeY,
                          const TrackerCellID& cellIDCoder,
                          int system_id,
                          double thr,
                          float fe_slope,
                          float starttime,
                          float t_step,
                          MsgStream& log);
    virtual ~ShapeProcessingSensor() {}

protected:
    vector<GridCoordinate> GetContour(const ClusterOfPixel& spot);

private:
    GridCoordinate GetNextPoint(GridCoordinate c, GridCoordinate p);
    GridPosition p_locate;
};

#endif //ShapeProcessingSensor_h