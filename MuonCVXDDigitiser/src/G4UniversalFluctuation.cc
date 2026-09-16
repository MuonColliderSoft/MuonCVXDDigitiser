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
// -------------------------------------------------------------------
//
// GEANT4 Class file
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
// 28-12-02 add method Dispersion (V.Ivanchenko)
// 07-02-03 change signature (V.Ivanchenko)
// 13-02-03 Add name (V.Ivanchenko)
// Modified for standalone use in ORCA. d.k. 6/04
//
// -------------------------------------------------------------------
 

#include "G4UniversalFluctuation.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandFlat.h"
#include <cmath>

using std::max;
using CLHEP::RandGaussQ;
using CLHEP::RandPoisson;
using CLHEP::RandFlat;
using CLHEP::twopi_mc2_rcl2;
using CLHEP::eV;
using CLHEP::keV;
using CLHEP::electron_mass_c2;
using CLHEP::proton_mass_c2;

// The constructor setups various constants pluc eloss parameters
// for silicon.   
G4UniversalFluctuation::G4UniversalFluctuation():
    chargeSquare(1.),           //Assume all particles have charge 1
    ipotFluct(0.0001736),       //GEANT4 (for Silicon): material->GetIonisation()->GetMeanExcitationEnergy();
    electronDensity(6.797E+20), //GEANT4 (for Silicon): material->GetElectronDensity();
    f1Fluct(0.8571),            //GEANT4 (for Silicon): material->GetIonisation()->GetF1fluct();
    f2Fluct(0.1429),            //GEANT4 (for Silicon): material->GetIonisation()->GetF2fluct();
    e1Fluct(0.000116),          //GEANT4 (for Silicon): material->GetIonisation()->GetEnergy1fluct();
    e2Fluct(0.00196),           //GEANT4 (for Silicon): material->GetIonisation()->GetEnergy2fluct();
    rateFluct(0.4),             //GEANT4 (for Silicon): material->GetIonisation()->GetRateionexcfluct();
    e1LogFluct(-9.063),         //GEANT4 (for Silicon): material->GetIonisation()->GetLogEnergy1fluct();
    e2LogFluct(-6.235),         //GEANT4 (for Silicon): material->GetIonisation()->GetLogEnergy2fluct();
    ipotLogFluct(-8.659),       //GEANT4 (for Silicon): material->GetIonisation()->GetLogMeanExcEnergy();
    e0(1.E-5),                  //GEANT4 (for Silicon): material->GetIonisation()->GetEnergy0fluct();
    minNumberInteractionsBohr(10.0),
    theBohrBeta2(50.0 * keV / proton_mass_c2),
    minLoss(0.000001 * eV),
    problim(0.01),
    sumalim(-log(0.01)),
    alim(10.),
    nmaxCont1(4),
    nmaxCont2(16)
{}

G4UniversalFluctuation::~G4UniversalFluctuation()
{}

// The main dedx fluctuation routine.
// Arguments: momentum in MeV/c, mass in MeV, delta ray cut (tmax) in
// MeV, silicon thickness in mm, mean eloss in MeV. 
double G4UniversalFluctuation::SampleFluctuations(const double momentum,
                                                  const double mass,
                                                  const double tmax,
                                                  const double length,
                                                  const double meanLoss)
{
    //  calculate actual loss from the mean loss
    //  The model used to get the fluctuation is essentially the same
    // as in Glandz in Geant3.

    // shortcut for very very small loss 
    if(meanLoss < minLoss) return meanLoss;

    double gam2 = pow(momentum, 2) / pow(mass, 2) + 1.0;
    double beta2 = 1.0 - 1.0 / gam2;

    // Validity range for delta electron cross section
    double loss, siga;

    // Gaussian fluctuation 
    if (meanLoss >= minNumberInteractionsBohr * tmax || tmax <= ipotFluct * minNumberInteractionsBohr)
    {
        siga = (1.0 / beta2 - 0.5) * twopi_mc2_rcl2 * tmax * length * electronDensity * chargeSquare ;
        siga = sqrt(siga);
        do
        {
            loss = RandGaussQ::shoot(meanLoss, siga);
        }
        while (loss < 0. || loss > 2. * meanLoss);

        return loss;
    }

    // Non Gaussian fluctuation 
    double w1 = tmax / ipotFluct;
    double w2 = log(2. * electron_mass_c2 * (gam2 - 1.0));

    double C = meanLoss * (1. - rateFluct) / (w2 - ipotLogFluct - beta2);

    double a1 = C * f1Fluct * (w2 - e1LogFluct - beta2) / e1Fluct;
    double a2 = C * f2Fluct * (w2 - e2LogFluct - beta2) / e2Fluct;
    double a3 = rateFluct * meanLoss * (tmax - ipotFluct) / (ipotFluct * tmax * log(w1));
    if (a1 < 0.) a1 = 0.;
    if (a2 < 0.) a2 = 0.;
    if (a3 < 0.) a3 = 0.;

    double suma = a1 + a2 + a3;

    loss = 0. ;
    int p3 = 0;

    if (suma < sumalim)             // very small Step
    {
        if(tmax == ipotFluct)
        {
            a3 = meanLoss / e0;

            if (a3 > alim)
            {
                siga = sqrt(a3) ;
                p3 = max(0, int(RandGaussQ::shoot(a3, siga) + 0.5));
            }
            else
            {
                p3 = RandPoisson::shoot(a3);
            }
            loss = p3 * e0;

            if (p3 > 0)
            {
                loss += (1. - 2. * RandFlat::shoot()) * e0;
            }
        }
        else
        {
            double adacut = tmax - ipotFluct + e0;
            a3 = meanLoss * (adacut - e0) / (adacut * e0 * log(adacut / e0));

            if (a3 > alim)
            {
                siga = sqrt(a3);
                p3 = max(0, int(RandGaussQ::shoot(a3, siga) + 0.5));
            }
            else
            {
                p3 = RandPoisson::shoot(a3);
            }

            if (p3 > 0)
            {
                double w = (adacut - e0) / adacut;
                double corrfac = 1.;
                if (p3 > nmaxCont2)
                {
                    corrfac = double(p3) / double(nmaxCont2);
                    p3 = nmaxCont2;
                }

                for (int i = 0; i < p3; i++)
                {
                    loss += 1. / (1. - w * RandFlat::shoot());
                }
                loss *= e0 * corrfac;  
            }        
        }
    }
    else                              // not so small Step
    {
        // excitation type 1
        int p1 = 0;
        if (a1 > alim)
        {
            siga = sqrt(a1) ;
            p1 = max(0, int(RandGaussQ::shoot(a1, siga) + 0.5));
        }
        else
        {
            p1 = RandPoisson::shoot(a1);
        }

        // excitation type 2
        int p2 = 0;
        if (a2 > alim)
        {
            siga = sqrt(a2) ;
            p2 = max(0, int(RandGaussQ::shoot(a2, siga) + 0.5));
        }
        else
        {
            p2 = RandPoisson::shoot(a2);
        }
        loss = p1 * e1Fluct + p2 * e2Fluct;

        // smearing to avoid unphysical peaks
        if (p2 > 0)
        {
            loss += (1. - 2. * RandFlat::shoot()) * e2Fluct;   
        }
        else if (loss > 0.)
        {
            loss += (1. - 2. * RandFlat::shoot()) * e1Fluct;
        }   

        // ionisation .......................................
        if(a3 > 0.)
        {
            if (a3 > alim)
            {
                siga = sqrt(a3) ;
                p3 = max(0, int(RandGaussQ::shoot(a3, siga) + 0.5));
            }
            else
            {
                p3 = RandPoisson::shoot(a3);
            }

            if (p3 > 0)
            {
                double na = 0.; 
                double alfa = 1.;
                double d_p3 = double(p3);
                if (p3 > nmaxCont2)
                {
                    double rfac = d_p3 / (double(nmaxCont2 + p3));
                    na = RandGaussQ::shoot(d_p3 * rfac, double(nmaxCont1) * rfac);
                    if (na > 0.)
                    {
                        alfa = w1 * double(nmaxCont2 + p3) / (w1 * double(nmaxCont2) + d_p3);
                        double alfa1 = alfa * log(alfa ) / (alfa - 1.);
                        double ea = na * ipotFluct * alfa1;
                        double sea = ipotFluct * sqrt(na * (alfa - pow(alfa1, 2)));
                        loss += RandGaussQ::shoot(ea,sea);
                    }
                }

                int nb = int(d_p3 - na);
                if (nb > 0)
                {
                    w2 = alfa * ipotFluct;
                    double w  = (tmax - w2) / tmax;      
                    for (int k = 0; k < nb; k++)
                    {
                        loss +=  w2 / (1. - w * RandFlat::shoot());
                    }
                }
            }
        }
    } 

    return loss;
}
