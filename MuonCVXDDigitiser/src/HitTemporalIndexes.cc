#include "HitTemporalIndexes.h"

#include <algorithm>

using std::min;

HitTemporalIndexes::HitTemporalIndexes(const edm4hep::SimTrackerHitCollection& STHcol,
                                       const TrackerCellID& cellIDCoder):
    htable()
{
    for (const edm4hep::SimTrackerHit simTrkHit : STHcol)
    {
        int layer = cellIDCoder.layer(simTrkHit.getCellID());
        int ladder = cellIDCoder.module(simTrkHit.getCellID());
        htable[GetKey(layer, ladder)].push(simTrkHit);
    }
}

HitTemporalIndexes::~HitTemporalIndexes()
{}

std::optional<edm4hep::SimTrackerHit> HitTemporalIndexes::CurrentHit(int layer, int ladder)
{
    int tkey = GetKey(layer, ladder);
    auto item = htable.find(tkey);
    if (item != htable.end() && !item->second.empty())
    {
        return item->second.top();
    }
    return std::nullopt;
}

void HitTemporalIndexes::DisposeHit(int layer, int ladder)
{
    int tkey = GetKey(layer, ladder);
    auto item = htable.find(tkey);
    if (item != htable.end() && !item->second.empty())
    {
        item->second.pop();
    }
}

int HitTemporalIndexes::GetHitNumber(int layer, int ladder)
{
    int tkey = GetKey(layer, ladder);
    auto item = htable.find(tkey);
    return item != htable.end() ? item->second.size() : -1;
}

float HitTemporalIndexes::GetMinTime()
{
    float min_time { MAXTIME };
    for (auto& item : htable)
    {
        if (!item.second.empty()) min_time = min(min_time, item.second.top().getTime());
    }
    return min_time;
}

float HitTemporalIndexes::GetMinTime(int layer, int ladder)
{
    int tkey = GetKey(layer, ladder);
    auto item = htable.find(tkey);
    if (item != htable.end() && !item->second.empty()) return item->second.top().getTime();
    return MAXTIME;
}

int HitTemporalIndexes::GetKey(int layer, int ladder)
{
    // TODO use the full cellID
    // The module field of tracker cell IDs can exceed 1000 (11 bits in MAIA)
    return (layer << 16) + ladder;
}
