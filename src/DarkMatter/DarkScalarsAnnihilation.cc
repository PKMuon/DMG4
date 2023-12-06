/*
 * DarkScalarsAnnihilation.cc
 *
 *  Created on: Oct 29, 2020
 *      Author: celentan
 */

#include "DarkMatterAnnihilation.hh"
#include "DarkScalarsAnnihilation.hh"
#include "Utils.hh"

#include <iostream>
#include <cmath>

DarkScalarsAnnihilation::DarkScalarsAnnihilation(double MAIn, double EThreshIn, double SigmaNormIn, double ANuclIn, double ZNuclIn, double DensityIn,
    double epsilIn, int IDecayIn, double rIn, double alphaDIn, int IBranchingIn, double fIn) :
    DarkMatterAnnihilation(MAIn, EThreshIn, SigmaNormIn, ANuclIn, ZNuclIn, DensityIn, epsilIn, IDecayIn, rIn, alphaDIn, IBranchingIn, fIn) {
  DMType = 2; //A.C.

  std::cout << "Initialized DarkScalarsAnnihilation (e+ e- -> A' -> DM DM) for material density = " << DensityIn << std::endl;
  std::cout << std::endl;

}

DarkScalarsAnnihilation::~DarkScalarsAnnihilation() {
  ;
}

//Convenience private method to be shared among TotalCrossSectionCalc and GetSigmaMax.
//This is the total cross section without the BW denominator
//s: e+e- invariant mass squared
double DarkScalarsAnnihilation::PreFactor(double ss) {

  if (ss < this->sMin())
    return 0.; // A.C. e+e- -> S -> chi chi can happen also for an S and chi with large mass,
               // i.e. through the off-shell tail of the resonance, but this still needs to be kinematically allowed
  double qq = this->q(ss);

  double sigma = 4 * M_PI * alphaEW * epsilBench * epsilBench * alphaD;
  sigma = sigma * qq / sqrt(ss);

  switch (iBranchingType) {
  case 0:
    sigma = sigma * (2 * qq * qq); // A.C. this is for final state fermions (default)
    break;
  case 1:
    sigma = sigma * (MA * MA) / 4;   // A.C. this is for final state scalars
    break;
  }

  //here sigma is in  1 /Energy^2. Move to pBarn;
  sigma = sigma * GeVtoPb;

  return sigma;
}

//E0: positron TOTAL energy in lab frame
double DarkScalarsAnnihilation::GetSigmaTot(double E0) {
  return TotalCrossSectionCalc(E0);
}



double DarkScalarsAnnihilation::CrossSectionDSDX(double XEv, double E0) {
  (void) (E0);
  if (XEv > 0.9999)
    return 1.;
  return 0.;
}

double DarkScalarsAnnihilation::CrossSectionDSDXDU(double XEv, double UThetaEv, double E0) {
  (void) (UThetaEv);
  (void) (E0);
  if (XEv > 0.9999)
    return 1.;
  return 0.;
}

double DarkScalarsAnnihilation::Width() {
  double ret;
  ret = MA * epsil * epsil * alphaEW * 1. / 2;
  if (MA / 2. > mChi) {
    switch (iBranchingType) {
    case 0:
      ret += MA * alphaD * (1 - 4 * mChi * mChi / (MA * MA)) * sqrt(1 - 4 * mChi * mChi / (MA * MA)) / 2.; // A.C. this is for final state fermions (default)
      break;
    case 1:
      ret += MA * alphaD * sqrt(1 - 4 * mChi * mChi / (MA * MA)) / 4.; // A.C. this is for final state scalars
      break;
    }
  }
  return ret;
}

void DarkScalarsAnnihilation::SetMA(double MAIn) {
  std::cout << "DarkScalarsAnnihilation::SetMA was called with MAIn = " << MAIn << std::endl;

  if (iBranchingType == 2) {
    mChi1 = MA * r;
    mChi2 = (1. + f) * mChi1;
  } else {
    mChi = MA * r;
    mChi1 = mChi;
    mChi2 = mChi;
  }
}

//e+e- --> S --> LDM LDM is uniform in cos(theta) CM, for both scalar and fermionic LDM
//E0: positron TOTAL energy in lab frame
double DarkScalarsAnnihilation::AngularDistributionResonant(double eta,double E0){
  double ret=1;
  return ret;
}
