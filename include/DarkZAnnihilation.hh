#pragma once

/*
 * DarkZAnnihilation.hh
 *
 *  Created on: Oct 6, 2020
 *      Author: celentan
 *
 *  Modified on: Sep 27, 2023
 *      Author: bbantoob
 *      - Change to inherit from DarkMatterAnnihilation class
 *      - Included missing implementations from abstract methods in new class
 *      namely PreFactor
 */

#include "DarkMatterAnnihilation.hh"


class DarkZAnnihilation: public DarkMatterAnnihilation
{

public:

    DarkZAnnihilation(double MAIn, double EThreshIn, double SigmaNormIn = 1., double ANuclIn = 207., double ZNuclIn = 82., double DensityIn = 11.35, double epsilIn = 0.0001, int IDecayIn = 0,double rIn=1./3,double alphaD=0.5,int IBranchingIn=0, double fIn=0.1);
    virtual ~DarkZAnnihilation();

    virtual double PreFactor(double s);
    virtual double GetSigmaTot(double E0);
    virtual double sMin();
    virtual double q(double s);
    virtual double CrossSectionDSDX(double Xev, double E0);
    virtual double CrossSectionDSDXDU(double Xev, double UThetaEv, double E0);
    virtual double Width();
    virtual void SetMA(double MAIn);

    virtual double AngularDistributionResonant(double eta,double E0);

private:
    double piF2(double q2); //q2 in GeV^2

};

