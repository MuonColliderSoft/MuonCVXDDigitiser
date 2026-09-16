#include "PixelDigiMatrix.h"
#include <cmath>
#include "streamlog/streamlog.h"

PixelDigiMatrix::PixelDigiMatrix(int layer,
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
                                 float t_step):
    AbstractSensor( layer,
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
                    starttime,
                    t_step),
    delta_c(t_step * fe_slope),
    clock_cnt(0),
    pixels(),
    expir_table(),
    charge_buffer(),
    start_table()
{}

PixelDigiMatrix::~PixelDigiMatrix()
{}

void PixelDigiMatrix::Reset()
{
    pixels.clear();
}

void PixelDigiMatrix::BeginClockStep()
{
    if (expir_table.find(clock_cnt) != expir_table.end())
    {
        for (auto pixpos : expir_table[clock_cnt])
        {
            pixels.erase(pixpos.second);
        }
        expir_table.erase(clock_cnt);
    }

    charge_buffer.clear();
    start_table.clear();
}

void PixelDigiMatrix::UpdatePixel(int x, int y, float chrg)
{
    if (!check(x, y)) return;

    LinearPosition lpos = l_locate(x, y);
    auto cb_iter = charge_buffer.find(lpos);
    if (cb_iter == charge_buffer.end())
    {
        charge_buffer.emplace(lpos, chrg);
    }
    else
    {
        charge_buffer[lpos] = cb_iter->second + chrg;
    }
}

void PixelDigiMatrix::EndClockStep()
{
    clock_cnt += 1;

    if (charge_buffer.size() == 0) return;

    for (auto chrg_item : charge_buffer)
    {
        LinearPosition lpos = chrg_item.first;
        float chrg = chrg_item.second;

        LinearPosition s_pos = sensor_for_pixel(lpos);
        ClockTicks pix_expir = 0;

        auto pix = pixels.find(lpos);
        if (pix == pixels.end())
        {
            pix_expir = calc_end_clock(chrg);
            if (pix_expir == 0) continue;

            pix_expir += clock_cnt;
            PixelRawData p_data { chrg, clock_cnt, pix_expir };
            pixels.emplace(lpos, p_data);

            start_table.emplace(s_pos, lpos);
        }
        else
        {
            // Pixel pile-up
            float new_chrg = (pix->second).charge + chrg;
            pixels[lpos].charge = new_chrg;

            pix_expir = calc_end_clock(new_chrg);
            if (pix_expir == 0) continue;

            pix_expir += pixels[lpos].t_begin;
            ClockTicks prev_expir = pixels[lpos].t_end;
            auto d_range = expir_table[prev_expir].equal_range(s_pos);
            for (auto it = d_range.first; it != d_range.second; it++)
            {
                if (it->second == lpos)
                {
                    expir_table[prev_expir].erase(it);
                    break;
                }
            }

            pixels[lpos].t_end = pix_expir;

        }

        if (expir_table.find(pix_expir) == expir_table.end())
        {
            SensorBin sBin {};
            expir_table.emplace(pix_expir, sBin);
        }
        expir_table[pix_expir].emplace(s_pos, lpos);
    }
}

PixelData PixelDigiMatrix::GetPixel(int x, int y)
{
    if (status != MatrixStatus::ok)
    {
        return { 0, 0, PixelStatus::geometry_error };
    }
    if (!check(x, y))
    {
        return { 0, 0, PixelStatus::out_of_bounds };
    }

    LinearPosition lpos = l_locate(x, y);
    auto pstat = calc_status(lpos);

    PixelData result { 0, 0, pstat };

    if (pstat != PixelStatus::off)
    {
        result.time = init_time + pixels[lpos].t_begin * clock_step;
    }
    if  (pstat == PixelStatus::ready)
    {
        result.charge = (pixels[lpos].t_end - pixels[lpos].t_begin) * delta_c;
    }
    return result;
}

bool PixelDigiMatrix::IsActive()
{
    return pixels.size() > 0;
}

bool PixelDigiMatrix::CheckStatus(int x, int y, PixelStatus pstat)
{
    if (!check(x, y)) return pstat == PixelStatus::out_of_bounds;

    return pstat == calc_status(l_locate(x, y));
}

bool PixelDigiMatrix::CheckStatusOnSensor(int seg_x, int seg_y, PixelStatus pstat)
{
    if (pstat == PixelStatus::ready and expir_table.find(clock_cnt) != expir_table.end())
    {
        auto d_range = expir_table[clock_cnt].equal_range(s_locate(seg_x, seg_y));
        return d_range.first != d_range.second;
    }
    if (pstat == PixelStatus::start)
    {
        auto d_range = start_table.equal_range(s_locate(seg_x, seg_y));
        return d_range.first != d_range.second;
    }

    return false;
}

vector<LocatedPixel> PixelDigiMatrix::GetPixelsFromSensor(int seg_x, int seg_y, PixelStatus pstat)
{
    vector<LocatedPixel> result;

    if (pstat == PixelStatus::ready and expir_table.find(clock_cnt) != expir_table.end())
    {
        auto d_range = expir_table[clock_cnt].equal_range(s_locate(seg_x, seg_y));
        for (auto it = d_range.first; it != d_range.second; it++)
        {
            GridCoordinate g_pos = l_locate(it->second);
            LocatedPixel l_pix
            {
                LadderRowToSensorRow(g_pos.row, seg_x),
                LadderColToSensorCol(g_pos.col, seg_y),
                GetPixel(g_pos.row, g_pos.col)
            };
            result.push_back(l_pix);
        }
    }

    if (pstat == PixelStatus::start)
    {
        auto d_range = start_table.equal_range(s_locate(seg_x, seg_y));
        for (auto it = d_range.first; it != d_range.second; it++)
        {
            GridCoordinate g_pos = l_locate(it->second);
            LocatedPixel l_pix
            {
                LadderRowToSensorRow(g_pos.row, seg_x),
                LadderColToSensorCol(g_pos.col, seg_y),
                GetPixel(g_pos.row, g_pos.col)
            };
            result.push_back(l_pix);
        }
    }

    return result;
}

PixelStatus PixelDigiMatrix::calc_status(LinearPosition lpos)
{
    auto pix = pixels.find(lpos);
    if ( pix == pixels.end()) return PixelStatus::off;

    if ((pix->second).t_begin == clock_cnt) return PixelStatus::start;

    if ((pix->second).t_end == clock_cnt) return PixelStatus::ready;

    return PixelStatus::on;
}

ClockTicks PixelDigiMatrix::calc_end_clock(float charge)
{
    // TODO missing dispersion
    if (charge <= _thr_level) return 0;

    return std::ceil((charge - _thr_level) / delta_c);
}

LinearPosition PixelDigiMatrix::sensor_for_pixel(LinearPosition pos)
{
    return s_locate((pos / l_columns) / s_rows, (pos % l_columns) / s_colums);
}
