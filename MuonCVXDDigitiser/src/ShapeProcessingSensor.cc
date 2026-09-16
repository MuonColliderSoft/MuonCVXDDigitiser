#include "ShapeProcessingSensor.h"

#include "streamlog/streamlog.h"

//TODO remove
#include <sstream>

ShapeProcessingSensor::ShapeProcessingSensor(int layer,
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
                                             float t_step) :
    HKBaseSensor(layer,
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
    p_locate({ GetSensorRows(), GetSensorCols() })
{}

/* ***************************************************************************************
 * 
 * Radial sweep algorithm
 * TODO missing holes detection in cluster
 * 
 * ************************************************************************************ */
GridCoordinate ShapeProcessingSensor::GetNextPoint(GridCoordinate c, GridCoordinate p)
{
    int idx = (c.row - p.row) + (c.col - p.col) * 3;
    switch(idx)
    {
    case 1:
        return { c.row - 1, c.col + 1 };
    case -2:
        return { c.row,     c.col + 1 };
    case -3:
        return { c.row + 1, c.col + 1 };
    case -4:
        return { c.row + 1, c.col     };
    case -1:
        return { c.row + 1, c.col - 1 };
    case 2:
        return { c.row    , c.col - 1 };
    case 3:
        return { c.row - 1, c.col - 1 };
    case 4:
        return { c.row - 1, c.col     };
    }
    
    //unreachable statement
    return { 0, 0 };
}

vector<GridCoordinate> ShapeProcessingSensor::GetContour(const ClusterOfPixel& spot)
{
    vector<GridCoordinate> result;
    std::stringstream logstr;

    int row_min;
    int row_max;
    int col_min;
    int col_max;
    tie(row_min, row_max, col_min, col_max) = GetBound(spot, p_locate);
    logstr << "Box for cluster: " << row_min << ":" << col_min  << " " ;
    logstr << row_max << ":" << col_max << " Size: " << spot.size() << std::endl;

    // It's better to put an empty frame around the spot
    int b_rows = row_max - row_min + 3;
    int b_cols = col_max - col_min + 3;

    vector<char> buffer;
    buffer.assign(b_rows * b_cols, 0);

    GridCoordinate curr_p { 0, 0 };
    GridCoordinate prev_p { 0, 0 };
    for (LinearPosition item : spot)
    {
        GridCoordinate a_coord = p_locate(item);
        int l_row = a_coord.row - row_min + 1;
        int l_col = a_coord.col - col_min + 1;
        buffer[l_row * b_cols + l_col] = 1;

        // Starting points
        if (curr_p.row == 0 && a_coord.row == row_min)
        {
            curr_p = { l_row, l_col };
            prev_p = { l_row - 1, l_col };
        }
    }

    int v_cnt = 0;
    do
    {
        logstr << "    " << curr_p.row << ":" << curr_p.col;
        if (v_cnt % 10 == 9) logstr << std::endl ;
        v_cnt++;

        result.push_back(curr_p);

        GridCoordinate next_p = prev_p;
        do
        {
            next_p = GetNextPoint(curr_p, next_p);
        }
        while (!buffer[next_p.row * b_cols + next_p.col]);

        prev_p = curr_p;
        curr_p = next_p;
    }
    while (!(result.size() > 1 && result[0] == prev_p && result[1] == curr_p));

    result.pop_back();
    logstr << std::endl << "Contour size: " << result.size() << std::endl;

    if(streamlog::out.write<streamlog::DEBUG7>())
#pragma omp critical
    {
        streamlog::out() << logstr.str() << std::endl;
    }

    return result;
}
