/*
 * DarkAxialsAnnihilation.cc
 *
 *  Created on: Oct 6, 2020
 *      Author: celentan
 *  Fixed: Nov 2, 2020
 */

#include "DarkMatter.hh"
#include "DarkMatterAnnihilation.hh"
#include "DarkAxialsAnnihilation.hh"
#include "Utils.hh"

#include <iostream>
#include <cmath>

DarkAxialsAnnihilation::DarkAxialsAnnihilation(double MAIn, double EThreshIn, double SigmaNormIn, double ANuclIn, double ZNuclIn, double DensityIn,
                                               double epsilIn, int IDecayIn, double rIn, double alphaDIn, int IBranchingIn, double fIn) :
                                               DarkMatterAnnihilation(MAIn, EThreshIn, SigmaNormIn, ANuclIn, ZNuclIn, DensityIn, epsilIn, IDecayIn, rIn,alphaDIn, IBranchingIn,fIn)
{
    DMType = 3; //A.C.

    std::cout << "Initialized DarkAxialAnnihilation (e+ e- -> A' -> DM DM) for material density = " << DensityIn << std::endl;
    std::cout << std::endl;
}

DarkAxialsAnnihilation::~DarkAxialsAnnihilation() {
    ;
}

//Convenience private method to be shared among TotalCrossSectionCalc and GetSigmaMax.
//This is the total cross section without the BW denominator
double DarkAxialsAnnihilation::PreFactor(double E0) {
    double ss = 2. * Mel * E0;

    if (sqrt(ss) < 2. * mChi)
        return 0.;   // A.C. e+e- -> A' -> chi chi can happen also for an A' and chi with large mass,
                     // i.e. through the off-shell tail of the resonance, but this still needs to be kinematically allowed
    double qq = sqrt(ss) / 2. * sqrt(1 - 4 * mChi * mChi / (ss));

    double sigma = 4 * M_PI * alphaEW * epsilBench * epsilBench * alphaD;
    sigma = sigma * qq / sqrt(ss);

    sigma = sigma * (8. / 3. * qq * qq); // A.C. this is for final state fermions (default)

    //here sigma is 1 /Energy^2. Move to pBarn;
    sigma = sigma * GeVtoPb;

    //A.C. correct here for atomic effects
    sigma = sigma * ZNucl;
    return sigma;
}

double DarkAxialsAnnihilation::GetSigmaTot(double E0) {
    return TotalCrossSectionCalc(E0);
}

bool DarkAxialsAnnihilation::EmissionAllowed(double E0, double DensityMat) // Different kinematic limit here
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

double DarkAxialsAnnihilation::CrossSectionDSDX(double XEv, double E0) {
    (void)(E0);
    if (XEv > 0.9999)
        return 1.;
    return 0.;
}

double DarkAxialsAnnihilation::CrossSectionDSDXDU(double XEv, double UThetaEv, double E0) {
    (void)(UThetaEv);
    (void)(E0);
    if (XEv > 0.9999)
        return 1.;
    return 0.;
}

double DarkAxialsAnnihilation::Width() {
    double ret;
    ret = MA * epsil * epsil * alphaEW * 1. / 3;
    if (MA / 2. > mChi) {
        ret += MA * alphaD * (1 - 4 * mChi * mChi / (MA * MA)) * sqrt(1 - 4 * mChi * mChi / (MA * MA)) / 3.;
    }
    return ret;
}

void DarkAxialsAnnihilation::SetMA(double MAIn) {
    std::cout << "DarkAxialsAnnihilation::SetMA was called with MAIn = " << MAIn << std::endl;
    if (iBranchingType == 2) {
        mChi1 = MA * r;
        mChi2 = (1. + f) * mChi1;
    } else {
        mChi = MA * r;
        mChi1 = mChi;
        mChi2 = mChi;
    }
}

double DarkAxialsAnnihilation::AngularDistributionResonant(double eta,double E0){
  double val=1;
  double ss = 2. * Mel * E0;

  switch (iBranchingType){
      case 0:
          //Fermionic LDM. Angular distribution f(eta) ~ 1+eta*eta. Max for eta=+-1
          if (sqrt(ss) < 2.*mChi){
              printf("DarkScalarsAnnihilation::AngularDistribution error with threshold, E0=%f, m=%f\n",E0,mChi);
              exit(1);
          }
          val=(1+eta*eta)/2; //Must be maximum == 1
          break;
      case 1:
          //Scalar LDM TODO
          val=1;
          break;
      default:
          break;
  }
  return val;
}

