#ifndef DetElemSlidingWindow_h
#define DetElemSlidingWindow_h 1

#include <list>

#include "HitTemporalIndexes.h"
#include "AbstractSensor.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"
#include "G4UniversalFluctuation.h"

#include <UTIL/CellIDDecoder.h>

using dd4hep::rec::SurfaceMap;
using UTIL::CellIDDecoder;

struct TimedSignalPoint
{
    double x;
    double y;
    double sigmaX;
    double sigmaY;
    double charge;
    SimTrackerHit* sim_hit;
};

typedef std::list<TimedSignalPoint> TimedSignalPointList;

class DetElemSlidingWindow
{
public:
    DetElemSlidingWindow(HitTemporalIndexes& htable,
                         AbstractSensor& sensor,
                         float wsize,
                         float starttime,
                         double tanLorentzAngleX,
                         double tanLorentzAngleY,
                         double cutOnDeltaRays,
                         double diffusionCoefficient,
                         double electronsPerKeV,
                         double segmentLength,
                         double energyLoss,
                         double widthOfCluster,
                         double electronicNoise,
                         double maxTrkLen,
                         double maxEnergyDelta,
                         const SurfaceMap* s_map,
                         bool zSegmented);
    virtual ~DetElemSlidingWindow();
    bool active();
    int process();
    float get_time();

private:
    void StoreSignalPoints(SimTrackerHit* hit);
    void UpdatePixels();
    double randomTail( const double qmin, const double qmax );

    float curr_time;
    float time_click;

    HitTemporalIndexes& _htable;
    AbstractSensor& _sensor;
    double _tanLorentzAngleX;
    double _tanLorentzAngleY;
    double _cutOnDeltaRays;
    double _diffusionCoefficient;
    double _electronsPerKeV;
    double _segmentLength;
    double _energyLoss;
    double _widthOfCluster;
    double _electronicNoise;
    double _maxTrkLen;
    double _deltaEne;
    TimedSignalPointList signals;
    const SurfaceMap* surf_map;
    CellIDDecoder<SimTrackerHit> cell_decoder;
    G4UniversalFluctuation* _fluctuate;
    bool _zSegmented;
};






#endif //DetElemSlidingWindow_h
