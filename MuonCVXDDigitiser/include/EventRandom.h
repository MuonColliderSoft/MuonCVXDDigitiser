#ifndef EventRandom_h
#define EventRandom_h 1

#include <cstddef>

#include "CLHEP/Random/MixMaxRng.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandPoisson.h"

/** Random number source for the digitisation of one event.
 *
 *  Create one per event, seeded from the event and run number (IUniqueIDGenSvc),
 *  and pass it down to every helper. Nothing is shared between events, so the
 *  output does not depend on the scheduling of events across threads.
 *
 *  The Gaussian generator is an instance rather than the static RandGauss::shoot(),
 *  whose cached second deviate is static thread-local state that would otherwise
 *  leak from one event into the next one processed on the same thread.
 */
class EventRandom
{
public:
    explicit EventRandom(std::size_t seed) :
        m_engine(static_cast<long>(seed)),
        m_gauss(m_engine)
    {}

    EventRandom(const EventRandom&) = delete;
    EventRandom& operator=(const EventRandom&) = delete;

    double gauss(double mean, double stdDev) { return m_gauss.fire(mean, stdDev); }
    double flat() { return CLHEP::RandFlat::shoot(&m_engine); }
    double flat(double low, double high) { return CLHEP::RandFlat::shoot(&m_engine, low, high); }
    long poisson(double mean) { return CLHEP::RandPoisson::shoot(&m_engine, mean); }

    /// For code taking a CLHEP engine directly (e.g. the fluctuation models)
    CLHEP::HepRandomEngine& engine() { return m_engine; }

private:
    CLHEP::MixMaxRng m_engine;
    CLHEP::RandGauss m_gauss;
};

#endif //EventRandom_h
