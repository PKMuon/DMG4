/*
 * DarkAxialsAnnihilation.cc
 *
 *  Created on: Oct 6, 2020
 *      Author: celentan
 *  Fixed: Nov 2, 2020
 */

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
//s: e+e- invariant mass GeV^2
double DarkAxialsAnnihilation::PreFactor(double ss) {

    if (ss<this->sMin())
        return 0.;   // A.C. e+e- -> A' -> chi chi can happen also for an A' and chi with large mass,
                     // i.e. through the off-shell tail of the resonance, but this still needs to be kinematically allowed
    double qq = this->q(ss);

    double sigma = 4 * M_PI * alphaEW * epsilBench * epsilBench * alphaD;
    sigma = sigma * qq / sqrt(ss);

    if (iBranchingType==0){
      sigma = sigma * (8. / 3. * qq * qq); // A.C. this is for final state fermions (default)
    }
    else if (iBranchingType==1){
      sigma = sigma * (2. / 3. * qq * qq); // A.C. this is for final state fermions (default)
    }
    //here sigma is 1 /Energy^2. Move to pBarn;
    sigma = sigma * GeVtoPb;

    return sigma;
}

//E0: positron TOTAL energy in lab frame
double DarkAxialsAnnihilation::GetSigmaTot(double E0) {
    return TotalCrossSectionCalc(E0);
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
      if (iBranchingType==0){
        ret += MA * alphaD * (1 - 4 * mChi * mChi / (MA * MA)) * sqrt(1 - 4 * mChi * mChi / (MA * MA)) / 3.;
      }
      else if (iBranchingType==1){
        ret += MA * alphaD * (1 - 4 * mChi * mChi / (MA * MA)) * sqrt(1 - 4 * mChi * mChi / (MA * MA)) / 12.;
      }
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

//E0: positron TOTAL energy in lab frame
//TODO: atomic effects are probably breaking the "ss"
double DarkAxialsAnnihilation::AngularDistributionResonant(double eta,double E0){
  double val=1;
  double ss = 2. * Mel * E0+2*Mel*Mel;


  if (ss < this->sMin()){
    printf("DarkScalarsAnnihilation::AngularDistribution error with threshold, E0=%f, m=%f\n",E0,mChi);
    exit(1);
  }

  switch (iBranchingType){
      case 0:
          //Fermionic LDM. Angular distribution f(eta) ~ 1+eta*eta. Max for eta=+-1
          val=(1+eta*eta)/2; //Must be maximum == 1
          break;
      case 1:
          //Scalar LDM. Angular distribution f(eta) ~ 1-eta*eta. Max for eta=0
          val=(1-eta*eta);
          break;
      default:
          break;
  }
  return val;
}

