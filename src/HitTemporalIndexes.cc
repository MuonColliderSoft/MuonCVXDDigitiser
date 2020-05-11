#include "HitTemporalIndexes.h"

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
        return item->second->pop();
    }
}

int HitTemporalIndexes::GetKey(int layer, int ladder)
{
    return layer * 1000 + ladder;
}

