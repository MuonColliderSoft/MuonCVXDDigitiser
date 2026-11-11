#ifndef DetElemSlidingWindow_h
#define DetElemSlidingWindow_h 1

#include <list>

#include "GaudiKernel/MsgStream.h"

#include "HitTemporalIndexes.h"
#include "AbstractSensor.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"
#include "G4UniversalFluctuation.h"

class EventRandom;

using dd4hep::rec::SurfaceMap;

struct TimedSignalPoint
{
    double x;
    double y;
    double sigmaX;
    double sigmaY;
    double charge;
    edm4hep::SimTrackerHit sim_hit;
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
                         double maxTrkLen,
                         double maxEnergyDelta,
                         const SurfaceMap* s_map,
                         bool zSegmented,
                         const G4UniversalFluctuation& fluctuate,
                         EventRandom& random,
                         MsgStream& log);
    virtual ~DetElemSlidingWindow();
    bool active();
    int process();
    float get_time();
    /// Number of SimTrackerHits skipped because they are not on their surface
    int get_off_surface() const { return off_surface; }

private:
    void StoreSignalPoints(const edm4hep::SimTrackerHit& hit);
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
    double _maxTrkLen;
    double _deltaEne;
    TimedSignalPointList signals;
    const SurfaceMap* surf_map;
    const G4UniversalFluctuation& _fluctuate;
    bool _zSegmented;
    EventRandom& _random;
    MsgStream& _log;
    int off_surface;
};

#endif //DetElemSlidingWindow_h
