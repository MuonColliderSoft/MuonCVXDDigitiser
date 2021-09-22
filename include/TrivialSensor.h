#ifndef TrivialSensor_h
#define TrivialSensor_h 1

#include "AbstractSensor.h"

class TrivialSensor : public AbstractSensor
{
public:
    TrivialSensor(int layer,
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

    ~TrivialSensor();

    void Reset() override;

    void BeginClockStep() override;

    void UpdatePixel(int x, int y, float chrg) override;

    void EndClockStep() override;

    PixelData GetPixel(int x, int y) override;

    bool IsActive() override;

    bool CheckStatus(int x, int y, PixelStatus pstat) override;

    void buildHits(SegmentDigiHitList& output) override;

private:
    vector<float> pixels;
    int charged_pix;
    int clock_cnt;
};

#endif //TrivialSensor_h