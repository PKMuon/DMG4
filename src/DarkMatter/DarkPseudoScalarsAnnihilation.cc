/*
 * DarkPseudoScalarsAnnihilation.cc
 *
 *  Created on: Oct 29, 2020
 *      Author: celentan
 */

#include "DarkMatterAnnihilation.hh"
#include "DarkPseudoScalarsAnnihilation.hh"
#include "Utils.hh"

#include <iostream>
#include <cmath>

DarkPseudoScalarsAnnihilation::DarkPseudoScalarsAnnihilation(double MAIn, double EThreshIn, double SigmaNormIn, double ANuclIn, double ZNuclIn, double DensityIn,
                                                             double epsilIn, int IDecayIn, double rIn, double alphaDIn, int IBranchingIn, double fIn) :
                                                             DarkMatterAnnihilation(MAIn, EThreshIn, SigmaNormIn, ANuclIn, ZNuclIn, DensityIn, epsilIn, IDecayIn, rIn,alphaDIn, IBranchingIn,fIn)
{
    DMType = 4; //A.C.

    std::cout << "Initialized DarkPseudoScalarsAnnihilation (e+ e- -> A' -> DM DM) for material density = " << DensityIn << std::endl;
    std::cout << std::endl;
}

DarkPseudoScalarsAnnihilation::~DarkPseudoScalarsAnnihilation() {
    ;
}


//Convenience private method to be shared among TotalCrossSectionCalc and GetSigmaMax.
//This is the total cross section without the BW denominator
//s: e+e- invariant mass squared
double DarkPseudoScalarsAnnihilation::PreFactor(double ss){

    if (ss < this->sMin())
        return 0.;   // A.C. e+e- -> S -> chi chi can happen also for an S and chi with large mass,
                   // i.e. through the off-shell tail of the resonance, but this still needs to be kinematically allowed
    double qq = this->q(ss);

    double sigma = 4 * M_PI * alphaEW * epsilBench * epsilBench * alphaD;
    sigma = sigma * qq / sqrt(ss);

    switch (iBranchingType) {
      case 0:
        sigma = sigma * (ss / 2);   // A.C. this is for final state fermions (default)
        break;
      case 1:
        sigma = sigma * (MA * MA) / 4;   // A.C. this is for final state scalars
        break;
    }
    //here sigma is in G4 internal units, 1 /Energy^2. Move to pBarn;
    sigma *= GeVtoPb;


    return sigma;
}

//E0: positron TOTAL energy in lab frame
double DarkPseudoScalarsAnnihilation::GetSigmaTot(double E0) {
    return TotalCrossSectionCalc(E0);
}


double DarkPseudoScalarsAnnihilation::CrossSectionDSDX(double XEv, double E0) {
    (void)(E0);
    if (XEv > 0.9999)
        return 1.;
    return 0.;
}

double DarkPseudoScalarsAnnihilation::CrossSectionDSDXDU(double XEv, double UThetaEv, double E0) {
    (void)(UThetaEv);
    (void)(E0);
    if (XEv > 0.9999)
        return 1.;
    return 0.;
}

double DarkPseudoScalarsAnnihilation::Width() {
    double ret;
    ret = MA * epsil * epsil * alphaEW * 1. / 2;
    if (MA / 2. > mChi) {
      if (iBranchingType==0){
        ret += MA * alphaD * sqrt(1 - 4 * mChi * mChi / (MA * MA)) / 2.;
      }
      if (iBranchingType==1){
             ret += MA * alphaD * sqrt(1 - 4 * mChi * mChi / (MA * MA)) / 4.;
      }
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

//E0: positron TOTAL energy in lab frame
double DarkPseudoScalarsAnnihilation::AngularDistributionResonant(double eta,double E0){
   double val=1;
   double ss = 2. * Mel * E0+2*Mel*Mel;

   //TODO annihilation
   if (ss < this->sMin()){
     printf("DarkScalarsAnnihilation::AngularDistribution error with threshold, E0=%f, m=%f\n",E0,mChi);
     exit(1);
   }
   switch (iBranchingType){
       case 0:
       case 1:
           //Fermionic and scalar LDM. Angular distribution f(eta) ~ 1
           val=1; //Must be maximum == 1
           break;
       default:
           break;
   }
   return val;
}


