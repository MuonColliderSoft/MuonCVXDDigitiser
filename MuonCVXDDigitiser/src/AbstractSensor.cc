#include "AbstractSensor.h"
#include <cmath>

#include <UTIL/ILDConf.h>

using lcio::ILDDetID;

AbstractSensor::AbstractSensor( int layer,
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
                                float t_step):
    _barrel_id(barrel_id),
    _layer(layer),
    _ladder(ladder),
    _thickness(thickness),
    _pixelSizeX(fabs(pixelSizeX)),
    _pixelSizeY(fabs(pixelSizeY)),
    _ladderLength(ladderLength > 0 ? ladderLength : 0),
    _ladderWidth(ladderWidth > 0 ? ladderWidth : 0),
    cellFmtStr(enc_str),
    _thr_level(thr),
    init_time(starttime),
    clock_step(t_step),
    l_locate({ 0, 0 }),
    s_locate({ 0, 0 }),
    status(MatrixStatus::ok),
    reset_simtable_at_once(true),
    simhit_table()
{
    // The ladder must be a whole number of pixels, within 1 um: the lengths are floats, e.g. a
    // 2528.4 mm ladder is stored as 2528.39990 mm, which truncating to 0.1 um units rejected
    const double tolerance = 1e-3; // mm
    l_rows = std::lround(ladderWidth / pixelSizeX);
    l_columns = std::lround(ladderLength / pixelSizeY);
    bool whole_pixels = std::abs(l_rows * pixelSizeX - ladderWidth) < tolerance
                        and std::abs(l_columns * pixelSizeY - ladderLength) < tolerance;

    x_segnum = xsegmentNumber;
    y_segnum = ysegmentNumber;

    s_rows = xsegmentNumber > 0 ? l_rows / xsegmentNumber : 1;
    s_colums = ysegmentNumber > 0 ? l_columns / ysegmentNumber : 1;
    
    if (!whole_pixels)
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
    // A position outside the sensor segment must not be read as a pixel of the next row or of the
    // neighbouring segment
    if (pos_x < 0 || pos_x >= s_rows || pos_y < 0 || pos_y >= s_colums)
        return pstat == PixelStatus::out_of_bounds;

    return CheckStatus(SensorRowToLadderRow(seg_x, pos_x),
                        SensorColToLadderCol(seg_y, pos_y),
                        pstat);
}

BitField64 AbstractSensor::getBFEncoder()
{
    BitField64 bf_encoder { cellFmtStr };
    bf_encoder.reset();
    bf_encoder[LCTrackerCellID::subdet()] = _barrel_id;
    bf_encoder[LCTrackerCellID::side()] = ILDDetID::barrel;
    bf_encoder[LCTrackerCellID::layer()] = _layer;
    bf_encoder[LCTrackerCellID::module()] = _ladder;
    return bf_encoder;
}

bool AbstractSensor::check(int x, int y)
{
    return (0 <= x and x < l_rows) and (0 <= y and y < l_columns);
}

void AbstractSensor::InitHitRegister()
{
    if (reset_simtable_at_once) simhit_table.clear();
}

void AbstractSensor::RegisterHit(int x, int y, SimTrackerHit* hit)
{
    simhit_table.emplace(l_locate(x, y), hit);
}

void AbstractSensor::fillInHitRelation(SimHitSet& sset, LinearPosition pos)
{
    auto hit_range = simhit_table.equal_range(pos);
    for (auto it = hit_range.first; it != hit_range.second; ++it)
    {
        sset.emplace(it->second);
    }

    if (!reset_simtable_at_once) simhit_table.erase(pos);
}
