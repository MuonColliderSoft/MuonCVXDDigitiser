#include "AbstractSensor.h"
#include <cmath>

// Value of the side field for barrel detectors (ILDDetID::barrel in LCIO)
static constexpr int BARREL_SIDE = 0;

AbstractSensor::AbstractSensor( int layer,
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
                                float starttime,
                                float t_step,
                                MsgStream& log):
    _system_id(system_id),
    _layer(layer),
    _ladder(ladder),
    _thickness(thickness),
    _pixelSizeX(fabs(pixelSizeX)),
    _pixelSizeY(fabs(pixelSizeY)),
    _ladderLength(ladderLength > 0 ? ladderLength : 0),
    _ladderWidth(ladderWidth > 0 ? ladderWidth : 0),
    _cellIDCoder(cellIDCoder),
    _thr_level(thr),
    init_time(starttime),
    clock_step(t_step),
    l_locate({ 0, 0 }),
    s_locate({ 0, 0 }),
    status(MatrixStatus::ok),
    reset_simtable_at_once(true),
    _log(log),
    simhit_table()
{
    int lwid = floor(ladderWidth * 1e4);
    int psx = floor(pixelSizeX * 1e4);
    int llen = floor(ladderLength * 1e4);
    int psy = floor(pixelSizeY * 1e4);

    l_rows = lwid / psx;
    l_columns = llen / psy;

    x_segnum = xsegmentNumber;
    y_segnum = ysegmentNumber;

    s_rows = xsegmentNumber > 0 ? l_rows / xsegmentNumber : 1;
    s_colums = ysegmentNumber > 0 ? l_columns / ysegmentNumber : 1;
    
    if (lwid % psx > 0 || llen % psy > 0)
    {
        status = MatrixStatus::pixel_number_error;
    }
    else if (xsegmentNumber <= 0 || ysegmentNumber <= 0)
    {
        status = MatrixStatus::segment_number_error;
    }
    else if (l_rows % xsegmentNumber > 0 || l_columns % ysegmentNumber > 0)
    {
        status = MatrixStatus::segment_number_error;
    }

    l_locate = GridPosition(l_rows, l_columns);
    s_locate = GridPosition(x_segnum, y_segnum);
}

AbstractSensor::~AbstractSensor()
{}

int AbstractSensor::XToPixelRow(double x)
{
    return int((x + _ladderWidth / 2) / _pixelSizeX);
}

int AbstractSensor::YToPixelCol(double y)
{
    return int((y + _ladderLength / 2) / _pixelSizeY);
}

double AbstractSensor::PixelRowToX(int ix)
{
    return ((0.5 + double(ix)) * _pixelSizeX) - _ladderWidth / 2;
}

double AbstractSensor::PixelColToY(int iy)
{
    return ((0.5 + double(iy)) * _pixelSizeY) - _ladderLength / 2;
}


PixelData AbstractSensor::getPixel(int seg_x, int seg_y, int pos_x, int pos_y)
{
    return GetPixel(SensorRowToLadderRow(seg_x, pos_x), SensorColToLadderCol(seg_y, pos_y));
}

bool AbstractSensor::checkStatus(int seg_x, int seg_y, int pos_x, int pos_y, PixelStatus pstat)
{
    // A position outside the sensor segment must not be confused with a pixel of the neighbouring segment
    if (pos_x < 0 || pos_x >= s_rows || pos_y < 0 || pos_y >= s_colums)
        return pstat == PixelStatus::out_of_bounds;

    return CheckStatus(SensorRowToLadderRow(seg_x, pos_x),
                        SensorColToLadderCol(seg_y, pos_y),
                        pstat);
}

std::uint64_t AbstractSensor::getLadderCellID()
{
    return _cellIDCoder.encode(_system_id, BARREL_SIDE, _layer, _ladder, 0);
}

bool AbstractSensor::check(int x, int y)
{
    return (0 <= x and x < l_rows) and (0 <= y and y < l_columns);
}

void AbstractSensor::InitHitRegister()
{
    if (reset_simtable_at_once) simhit_table.clear();
}

void AbstractSensor::RegisterHit(int x, int y, const edm4hep::SimTrackerHit& hit)
{
    simhit_table.emplace(l_locate(x, y), hit);
}

void AbstractSensor::fillInHitRelation(SimHitSet& sset, LinearPosition pos)
{
    auto hit_range = simhit_table.equal_range(pos);
    for (auto it = hit_range.first; it != hit_range.second; ++it)
    {
        sset.emplace(it->second.getObjectID().index, it->second);
    }

    if (!reset_simtable_at_once) simhit_table.erase(pos);
}
