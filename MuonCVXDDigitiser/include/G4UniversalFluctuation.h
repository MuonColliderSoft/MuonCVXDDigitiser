//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4UniversalFluctuation
//
// Author:        Vladimir Ivanchenko
// 
// Creation date: 03.01.2002
//
// Modifications:
//
// 13-05-20 thread-safety (P.Andreetto)
// 09-12-02 remove warnings (V.Ivanchenko)
// 28-12-02 add method Dispersion (V.Ivanchenko)
// 07-02-03 change signature (V.Ivanchenko)
// 13-02-03 Add name (V.Ivanchenko)
// Modified for standalone use in ORCA. d.k. 6/04
//
// Implementation of energy loss fluctuations
// -------------------------------------------------------------------
//

#ifndef G4UniversalFluctuation_h
#define G4UniversalFluctuation_h 

class G4UniversalFluctuation {
public:

    G4UniversalFluctuation();
    ~G4UniversalFluctuation();

    // momentum in MeV/c, mass in MeV, tmax (delta cut) in MeV, 
    // length in mm, meanLoss eloss in MeV.
    double SampleFluctuations(const double momentum,
                              const double mass,
                              const double tmax,
                              const double length,
                              const double meanLoss);

private:

    double chargeSquare;

    // data members to speed up the fluctuation calculation
    double ipotFluct;
    double electronDensity;
    //  G4double zeff;

    double f1Fluct;
    double f2Fluct;
    double e1Fluct;
    double e2Fluct;
    double rateFluct;
    double e1LogFluct;
    double e2LogFluct;
    double ipotLogFluct;
    double e0;

    double minNumberInteractionsBohr;
    double theBohrBeta2;
    double minLoss;
    double problim;
    double sumalim;
    double alim;
    int    nmaxCont1;
    int    nmaxCont2;
};

#endif //G4UniversalFluctuation_h

