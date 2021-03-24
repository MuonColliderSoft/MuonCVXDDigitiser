#ifndef PixelDigiMatrix_h
#define PixelDigiMatrix_h 1

#include <vector>
#include <functional>

using std::string;

enum class PixelStatus : char {
    ok,
    undefined,
    out_of_bounds,
    geometry_error
};

struct PixelData
{
    float charge;
    float time;
    PixelStatus status;
};

enum class MatrixStatus : char {
    ok,
    pixel_number_error,
    segment_number_error
};

struct SegmentDigiHit
{
    float x;
    float y;
    float charge;
    float time;
    int cellID0;
    int segment_x;  // redundant information
    int segment_y;  // redundant information
};

typedef std::vector<PixelData> EnergyMatrix;
typedef std::function<PixelData(PixelData pIn)> PixelTransformation;
typedef std::vector<SegmentDigiHit> SegmentDigiHitList;

class PixelDigiMatrix
{
public:
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
                    float s_level);
    virtual ~PixelDigiMatrix();

    virtual void buildHits(SegmentDigiHitList& output) = 0;

    inline int GetLayer() { return _layer; }
    inline int GetLadder() { return _ladder; }
    inline float GetThickness() { return _thickness; }
    inline float GetHalfThickness() { return _thickness / 2; }
    inline float GetLength() { return _ladderLength; }
    inline float GetHalfLength() { return _ladderLength / 2; }
    inline float GetWidth() { return _ladderWidth; }
    inline float GetHalfWidth() { return _ladderWidth / 2; }
    inline double GetPixelSizeX() { return _pixelSizeX; }
    inline double GetPixelSizeY() { return _pixelSizeY; }
    inline int GetLadderRows() { return l_rows; }
    inline int GetLadderCols() { return l_columns; }
    inline int GetSensorRows() { return s_rows; }
    inline int GetSensorCols() { return s_colums; }
    inline int GetSegNumX() { return x_segnum; }
    inline int GetSegNumY() { return y_segnum; }
    inline MatrixStatus GetStatus() { return status; }

    inline string GetCellIDFormatStr() { return cellFmtStr; }

    void Reset();
    void SetTime(float time);
    void UpdatePixel(int x, int y, float chrg);
    void Apply(PixelTransformation l_expr);
    PixelData GetPixel(int x, int y);
    float GetMaxCharge();

    inline int XToPixelRow(double x) { return int((x + _ladderWidth / 2) / _pixelSizeX); }
    inline int YToPixelCol(double y) { return int((y + _ladderLength / 2) / _pixelSizeY); }

    inline double PixelRowToX(int ix) { return ((0.5 + double(ix)) * _pixelSizeX) - _ladderWidth / 2; }
    inline double PixelColToY(int iy) { return ((0.5 + double(iy)) * _pixelSizeY) - _ladderLength / 2; }

protected:

    inline int SensorRowToLadderRow(int seg_x, int pos_x) { return seg_x * s_rows + pos_x; }
    inline int SensorColToLadderCol(int seg_y, int pos_y) { return seg_y * s_colums + pos_y; }
    PixelData GetPixel(int seg_x, int seg_y, int pos_x, int pos_y);

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
    float _satur_level;

private:
    inline int index(int x, int y) { return x * l_columns + y; }
    bool check(int x, int y);
    
    EnergyMatrix pixels;
    MatrixStatus status;
    float max_charge;
    bool  charge_valid;
};

#endif //PixelDigiMatrix_h
