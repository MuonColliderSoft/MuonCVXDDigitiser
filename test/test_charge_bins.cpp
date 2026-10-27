// Unit test of chargeDigitizationBins(): choice of the charge digitization table by number of
// bits, sub-detector and sensor thickness. The MAIA_v0 vertex sensors are all 50 um, whose table
// equals the generic one, so the thickness-dependent tables are only exercised here.

#include <iostream>
#include <string>
#include <vector>

#include "ChargeDigitizationBins.h"

namespace {

int failures = 0;

void check(bool condition, const std::string& what)
{
    if (!condition) {
        std::cerr << "FAILED: " << what << std::endl;
        ++failures;
    }
}

bool increasing(const std::vector<double>& bins)
{
    for (size_t i = 1; i < bins.size(); ++i)
        if (!(bins[i] > bins[i - 1])) return false;
    return true;
}

}  // namespace

int main()
{
    const std::vector<double> generic4 = chargeDigitizationBins(4, true, 0.050);

    // Number of bin edges and shape of every table. The 8-bit table has 266 edges instead of 256
    // (since it was introduced in 2021): kept as it is until the table is regenerated, at which
    // point this expectation must become 256.
    const std::vector<std::pair<int, size_t>> sizes = {{3, 8}, {4, 16}, {5, 32}, {6, 64}, {8, 266}};
    for (auto [bits, size] : sizes) {
        for (bool vertex : {true, false}) {
            const auto bins = chargeDigitizationBins(bits, vertex, 0.050);
            const std::string what = std::to_string(bits) + " bits, " + (vertex ? "vertex" : "tracker");
            check(bins.size() == size, what + ": " + std::to_string(size) + " bin edges");
            check(!bins.empty() && bins.front() == 500., what + ": first edge at 500 electrons");
            check(increasing(bins), what + ": increasing edges");
        }
    }

    // No table for an unsupported number of bits
    check(chargeDigitizationBins(7, true, 0.050).empty(), "7 bits: no table");
    check(chargeDigitizationBins(2, false, 0.050).empty(), "2 bits: no table");

    // 4-bit vertex tables by sensor thickness (last edge identifies the table)
    const std::vector<std::pair<double, double>> lastEdge = {
        {0.050, 29756}, {0.075, 44634}, {0.100, 59512}, {0.200, 119024}, {0.400, 238048}};
    for (auto [thickness, last] : lastEdge) {
        const auto bins = chargeDigitizationBins(4, true, thickness);
        check(bins.size() == 16 && bins.back() == last,
              "4 bits, vertex, " + std::to_string(thickness) + " mm: last edge " + std::to_string(last));
    }
    // The thickness comes from a float: matched on the value rounded to the micron
    check(chargeDigitizationBins(4, true, double(0.075f)).back() == 44634, "4 bits, vertex, 75 um as float");
    check(chargeDigitizationBins(4, true, 0.0754).back() == 44634, "4 bits, vertex, 75.4 um rounds to 75 um");
    // Thicknesses without a dedicated table keep the generic one
    check(chargeDigitizationBins(4, true, 0.060) == generic4, "4 bits, vertex, 60 um: generic table");
    // Only 4 bits has thickness-dependent tables
    check(chargeDigitizationBins(3, true, 0.100) == chargeDigitizationBins(3, true, 0.050),
          "3 bits, vertex: no thickness dependence");

    // Trackers: generic table, twice as wide around the first edge, for any thickness
    const auto tracker4 = chargeDigitizationBins(4, false, 0.100);
    bool doubled = tracker4.size() == generic4.size();
    for (size_t i = 0; doubled && i < generic4.size(); ++i)
        doubled = tracker4[i] == 2. * (generic4[i] - 500.) + 500.;
    check(doubled, "4 bits, tracker: generic table twice as wide");
    check(chargeDigitizationBins(4, false, 0.100) == chargeDigitizationBins(4, false, 0.050),
          "4 bits, tracker: no thickness dependence");

    if (failures == 0) std::cout << "All charge-bin checks passed" << std::endl;
    return failures == 0 ? 0 : 1;
}
