#ifndef ChargeDigitizationBins_h
#define ChargeDigitizationBins_h 1

#include <vector>

/** Bin edges (in electrons) of the variable-width charge digitization of one layer.
 *
 * The generic table for numBits (3, 4, 5, 6 or 8 bits). With 4 bits, vertex-detector layers use a
 * dedicated table for their sensor thickness (50, 75, 100, 200 or 400 um, matched on the thickness
 * rounded to the micron). For the inner and outer trackers the bins are twice as wide, around the
 * first edge. Empty when there is no table for numBits.
 *
 * @param numBits number of bits of the charge digitization
 * @param isVertex whether the layer belongs to the vertex detector
 * @param thicknessMM sensor thickness of the layer, in mm
 */
std::vector<double> chargeDigitizationBins(int numBits, bool isVertex, double thicknessMM);

#endif
