#ifndef DigitiserTypes_h
#define DigitiserTypes_h 1

#include <vector>

/** Point along the track path inside the sensor where ionisation charge is released.
 *  Local sensor coordinates in mm, eloss in DD4hep energy units (divide by dd4hep::keV for keV).
 */
struct IonisationPoint
{
    double x;
    double y;
    double z;
    double eloss;
};

/** Charge cloud arriving at the readout plane after drift and diffusion.
 *  Local sensor coordinates and widths in mm, charge in electrons.
 */
struct SignalPoint
{
    double x;
    double y;
    double sigmaX;
    double sigmaY;
    double charge;
};

typedef std::vector<IonisationPoint> IonisationPointVec;
typedef std::vector<SignalPoint> SignalPointVec;

#endif //DigitiserTypes_h
