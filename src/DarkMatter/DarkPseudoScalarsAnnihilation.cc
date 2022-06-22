/*
 * DarkPseudoScalarsAnnihilation.cc
 *
 *  Created on: Oct 29, 2020
 *      Author: celentan
 */

#include "DarkMatter.hh"
#include "DarkPseudoScalarsAnnihilation.hh"
#include "Utils.hh"

#include <iostream>
#include <cmath>

DarkPseudoScalarsAnnihilation::DarkPseudoScalarsAnnihilation(double MAIn, double EThreshIn, double SigmaNormIn, double ANuclIn, double ZNuclIn,
        double DensityIn, double epsilIn, int IDecayIn, int IBranchingIn, double rIn, double fIn, double alphaDIn) :
        DarkMatter(MAIn, EThreshIn, SigmaNormIn, ANuclIn, ZNuclIn, DensityIn, epsilIn, IDecayIn), iBranchingType(IBranchingIn), r(rIn), f(fIn), alphaD(alphaDIn) {
    DMType = 4; //A.C.
    ParentPDGID = -11;
    DaughterPDGID = 11;

    //default values
    r = rIn;
    f = fIn;
    mChi = MAIn / 3;
    mChi1 = mChi;
    mChi2 = mChi;

    if (iBranchingType != 0) {
        std::cerr << " DarkPseudoScalarAnnihilation iBranchingType!=0 not yet supported" << std::endl;
        exit(1);
    }

    if (iBranchingType == 2) {

        mChi1 = MA * r;
        mChi2 = (1. + f) * mChi1;
    } else {

        mChi = MA * r;
        mChi2 = mChi;
        mChi1 = mChi;
    }

    deltaMchi = mChi2 - mChi1;
    std::cout << "Initialized DarkPseudoScalarsAnnihilation (e+ e- -> A' -> DM DM) for material density = " << DensityIn << std::endl;
    std::cout << std::endl;
}

DarkPseudoScalarsAnnihilation::~DarkPseudoScalarsAnnihilation() {
    ;
}

//Input: E0, positron energy in GeV
//output: total annihilation cross-section in pbarn.
//Since the framework assumes this method is returning the total cross section per nucleous, for the moment I scale this by Z.
double DarkPseudoScalarsAnnihilation::TotalCrossSectionCalc(double E0) {
    double ss = 2. * Mel * E0;
    if (sqrt(ss) < 2. * mChi)
        return 0.;   // A.C. e+e- -> S -> chi chi can happen also for an S and chi with large mass,
                     // i.e. through the off-shell tail of the resonance, but this still needs to be kinematically allowed
    double qq = sqrt(ss) / 2. * sqrt(1 - 4 * mChi * mChi / (ss));
    double gg = this->Width();

    double sigma = 4 * M_PI * alphaEW * epsil * epsil * alphaD;
    sigma = sigma * qq / sqrt(ss);
    sigma = sigma / ((ss - MA * MA) * (ss - MA * MA) + MA * MA * gg * gg);

    sigma = sigma * (ss / 2);   // A.C. this is for final state fermions (default)

    //here sigma is in G4 internal units, 1 /Energy^2. Move to pBarn;
    sigma = sigma * GeVtoPb;

    //A.C. correct here for atomic effects
    sigma = sigma * ZNucl;
    return sigma;
}

double DarkPseudoScalarsAnnihilation::GetSigmaTot(double E0) {
    return TotalCrossSectionCalc(E0);
}

bool DarkPseudoScalarsAnnihilation::EmissionAllowed(double E0, double DensityMat) // Different kinematic limit here
        {

    if (sqrt(2. * Mel * E0) < 2. * mChi)
        return false;
    if (E0 < EThresh)
        return false;
    if (NEmissions)
        return false; // For G4 DM classes
    if (fabs(DensityMat - Density) > 0.1)
        return false;
    return true;
}

double DarkPseudoScalarsAnnihilation::CrossSectionDSDX(double XEv, double E0) {
    if (XEv > 0.9999)
        return 1.;
    return 0.;
}

double DarkPseudoScalarsAnnihilation::CrossSectionDSDXDU(double XEv, double UThetaEv, double E0) {
    if (XEv > 0.9999)
        return 1.;
    return 0.;
}

double DarkPseudoScalarsAnnihilation::Width() {
    double ret;
    ret = MA * epsil * epsil * alphaEW * 1. / 2;
    if (MA / 2. > mChi) {
        ret += MA * alphaD * sqrt(1 - 4 * mChi * mChi / (MA * MA)) / 2.;
    }
    return ret;
}

void DarkPseudoScalarsAnnihilation::SetMA(double MAIn) {
    std::cout << "DarkPseudoScalarsAnnihilation::SetMA was called with MAIn = " << MAIn << std::endl;
    if (iBranchingType == 2) {
        mChi1 = MA * r;
        mChi2 = (1. + f) * mChi1;
    } else {
        mChi = MA * r;
        mChi1 = mChi;
        mChi2 = mChi;
    }
}
