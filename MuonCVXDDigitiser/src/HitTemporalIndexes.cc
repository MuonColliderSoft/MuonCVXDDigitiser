#include "HitTemporalIndexes.h"

using std::min;

HitTemporalIndexes::HitTemporalIndexes(const LCCollection* STHcol):
    cellid_decoder(STHcol),
    htable()
{
    for (int i = 0; i < STHcol->getNumberOfElements(); ++i)
    {
        SimTrackerHit* simTrkHit = dynamic_cast<SimTrackerHit*>(STHcol->getElementAt(i));
        int layer = cellid_decoder(simTrkHit)["layer"];
        int ladder = cellid_decoder(simTrkHit)["module"];
        int tkey = GetKey(layer, ladder);

        auto item = htable.find(tkey);
        if (item == htable.end())
        {
            hit_queue* n_queue = new hit_queue();
            n_queue->push(simTrkHit);
            htable.emplace(tkey, n_queue);
        }
        else
        {
            item->second->push(simTrkHit);
        }
    }
}

HitTemporalIndexes::~HitTemporalIndexes()
{
    for (auto item : htable)
    {
        delete(item.second);
    }
}

SimTrackerHit* HitTemporalIndexes::CurrentHit(int layer, int ladder)
{
    int tkey = GetKey(layer, ladder);
    auto item = htable.find(tkey);
    if (item != htable.end() && !item->second->empty())
    {
        return item->second->top();
    }
    return nullptr;
}

void HitTemporalIndexes::DisposeHit(int layer, int ladder)
{
    int tkey = GetKey(layer, ladder);
    auto item = htable.find(tkey);
    if (item != htable.end() && !item->second->empty())
    {
        item->second->pop();
    }
}

int HitTemporalIndexes::GetHitNumber(int layer, int ladder)
{
    int tkey = GetKey(layer, ladder);
    auto item = htable.find(tkey);
    return item != htable.end() ? item->second->size() : -1;
}

float HitTemporalIndexes::GetMinTime()
{
    float min_time { MAXTIME };
    for (auto item : htable)
    {
        min_time = min(min_time, item.second->top()->getTime());
    }
    return min_time;
}

float HitTemporalIndexes::GetMinTime(int layer, int ladder)
{
    int tkey = GetKey(layer, ladder);
    auto item = htable.find(tkey);
    if (item != htable.end()) return item->second->top()->getTime();
    return MAXTIME;
}

int HitTemporalIndexes::GetKey(int layer, int ladder)
{
    // TODO use cellID0 and/or cellID1
    return layer * 1000 + ladder;
}

float HitTemporalIndexes::MAXTIME { std::numeric_limits<float>::max() };
