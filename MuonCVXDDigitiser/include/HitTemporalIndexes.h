#ifndef HitTemporalIndexes_h
#define HitTemporalIndexes_h 1

#include <limits>
#include <optional>
#include <queue>
#include <unordered_map>
#include <vector>

#include "edm4hep/SimTrackerHit.h"
#include "edm4hep/SimTrackerHitCollection.h"

#include "TrackerCellID.h"

using std::priority_queue;
using std::unordered_map;

class CmpTrackTime
{
public:
    CmpTrackTime(){}
    bool operator()(const edm4hep::SimTrackerHit& alfa, const edm4hep::SimTrackerHit& beta) const
    {
        return alfa.getTime() > beta.getTime();
    }
};

typedef priority_queue<edm4hep::SimTrackerHit, std::vector<edm4hep::SimTrackerHit>, CmpTrackTime> hit_queue;

/**
 * SimTrackerHits of a collection, grouped by layer and ladder and ordered by time.
 * Layers are identified by the layer ID of the cell ID.
 */
class HitTemporalIndexes
{
public:
    HitTemporalIndexes(const edm4hep::SimTrackerHitCollection& STHcol, const TrackerCellID& cellIDCoder);
    virtual ~HitTemporalIndexes();
    std::optional<edm4hep::SimTrackerHit> CurrentHit(int layer, int ladder);
    void DisposeHit(int layer, int ladder);
    int GetHitNumber(int layer, int ladder);
    float GetMinTime();
    float GetMinTime(int layer, int ladder);

    static constexpr float MAXTIME { std::numeric_limits<float>::max() };

private:
    inline int GetKey(int layer, int ladder);

    unordered_map<int, hit_queue> htable;
};

#endif //HitTemporalIndexes_h
