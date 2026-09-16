#ifndef PixelDigiMatrix_h
#define PixelDigiMatrix_h 1

#include <string>
#include <vector>
#include <unordered_map>

#include "AbstractSensor.h"

using std::string;
using std::unordered_map;
using std::unordered_multimap;
using std::vector;

struct LocatedPixel
{
    int row;
    int col;
    PixelData data;
};

using ClockTicks = int;

/**
 * @class PixelDigiMatrix
 * @brief Simulation of the chip RD53A
 *
 * This class implements the basic behaviour of the chip [RD53A](https://cds.cern.ch/record/2113263).
 * It simulates a matrix of pixels. Each pixel can collect a charge and perform a linear charge
 * depletion mechanism in order to measure the charge itself. A single threshold detection, with smearing,
 * is available for each pixel.
 * The matrix of pixels corresponds to the entire ladder; the ladder is divided into a grid of sensors.
 * Each matrix is identify by a couple of ID: the layer ID and the ladder ID.
 * This class must be operated by an agent which feeds it with the charge and synchronize the actions
 * through a clock.
 */
class PixelDigiMatrix : public AbstractSensor
{
public:
    /**
     * @brief The main constructor
     *
     * This constructor creates a matrix of pixels for a given ladder
     * within a layer of the vertex barrel.
     * @param layer The ID of the layer containing the pixel matrix
     * @param ladder The ID of the ladder matching this matrix of pixels
     * @param xsegmentNumber The number of sensors per ladder width
     * @param ysegmentNumber The number of sensors per ladder length
     * @param ladderLength The length of the ladder
     * @param ladderWidth The width of the ladder
     * @param thickness The thickness the ladder
     * @param pixelSizeX The width of a pixel
     * @param pixelSizeY The length of a pixel
     * @param enc_str The format string used to encode the CellID for any sensor of a ladder
     * @param barrel_id The ID of the vertex barrel inside the detector
     * @param thr The threshold for any pixel of the ladder
     * @param fe_slope The charge depletion slope of the FE
     * @param starttime The start time for the matrix evolution
     * @param t_step The clock period of the chip
     */
    PixelDigiMatrix(int layer,
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
                    float t_step);
    virtual ~PixelDigiMatrix();

    void Reset() override;

    void BeginClockStep() override;

    /**
     * @brief The charge aggregation call.
     *
     * This method must be called by the agent when a quantity of charge must be gathered
     * for a given pixel of the ladder
     * @param x The row number of the pixel
     * @param y The column number of the pixel
     * @param chrg The charge to be aggregated in the pixel
     */
    void UpdatePixel(int x, int y, float chrg) override;

    void EndClockStep() override;

    PixelData GetPixel(int x, int y) override;

    bool IsActive() override;

    bool CheckStatus(int x, int y, PixelStatus pstat) override;

protected:

    bool CheckStatusOnSensor(int seg_x, int seg_y, PixelStatus pstat);

    vector<LocatedPixel> GetPixelsFromSensor(int seg_x, int seg_y, PixelStatus pstat);

    float delta_c;
    int clock_cnt;

private:

    struct PixelRawData
    {
        float charge;
        ClockTicks t_begin;
        ClockTicks t_end;
    };

    using SensorBin = unordered_multimap<LinearPosition, LinearPosition>;

    PixelStatus calc_status(LinearPosition lpos);
    ClockTicks calc_end_clock(float charge);
    LinearPosition sensor_for_pixel(LinearPosition pos);

    unordered_map<LinearPosition, PixelRawData> pixels;
    unordered_map<ClockTicks, SensorBin> expir_table;
    unordered_map<LinearPosition, float> charge_buffer;
    SensorBin start_table; 
};

#endif //PixelDigiMatrix_h
