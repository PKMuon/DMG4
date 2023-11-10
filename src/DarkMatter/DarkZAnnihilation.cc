/*
 * DarkZAnnihilation.cc
 *
 *  Created on: Oct 6, 2023
 *      Author: celentan
 *  B-L implementation added on: Nov 10, 2023 by B. Banto
 *
 *  In this class, the parameter "epsilon" is actually the coupling "g" between Z' and mu-tau current (see 2206.03101 and 2207.09979)
 *
 *  The different BranchingTypes implemented refer to the following models:
 *  0: Lmu-Ltau vanilla (decay to neutrinos)
 *  1: Lmu-Ltau scalar DM
 *  10: B-L vanilla (decay to neutrinos)
 *  11: TODO: B-L scalar DM
 *  12: TODO: B-L fermionic DM
 *
 */

#include "DarkMatter.hh"
#include "DarkMatterAnnihilation.hh"
#include "DarkZAnnihilation.hh"
#include "Utils.hh"

#include <iostream>
#include <cmath>



DarkZAnnihilation::DarkZAnnihilation(double MAIn, double EThreshIn, double SigmaNormIn, double ANuclIn, double ZNuclIn, double DensityIn,
                                                 double epsilIn, int IDecayIn, double rIn, double alphaDIn, int IBranchingIn, double fIn) :
                                                 DarkMatterAnnihilation(MAIn, EThreshIn, SigmaNormIn, ANuclIn, ZNuclIn, DensityIn, epsilIn, IDecayIn, rIn,alphaDIn, IBranchingIn,fIn){
  DMType = 11; //A.C.
  ParentPDGID = -11;

  std::cout << "Initialized DarkZAnnihilation (e+ e- -> Z' -> DM DM) for material density = " << DensityIn << std::endl;
  std::cout << "mA: " << MA * 1E3 << " MeV " << std::endl;
  std::cout << "IBranchingType: " << iBranchingType << std::endl;
  if (iBranchingType == 0 || iBranchingType == 10) {
    std::cout << "decay to neutrinos, mass is negligible m_nu << mA" << std::endl;
  }
  else if (iBranchingType == 1 || iBranchingType == 11) {
    std::cout << "mChi1: " << mChi1 * 1E3 << " MeV " << std::endl;
    std::cout << "mChi2: " << mChi2 * 1E3 << " MeV " << std::endl;
  } else {
    std::cout << "mChi: " << mChi * 1E3 << " MeV " << std::endl;
  }
  std::cout << "Width: " << this->Width() * 1E3 << " MeV " << std::endl;

  std::cout << std::endl;
}

DarkZAnnihilation::~DarkZAnnihilation() {
  ;
}

//Convenience private method to be shared among TotalCrossSectionCalc and GetSigmaMax.
//This is the total cross section without the BW denominator
//E0: positron total energy in GeV
double DarkZAnnihilation::PreFactor(double E0) {

  double ss = 2. * Mel * E0 + 2*Mel*Mel;
  double E1 = 0., E2 = 0.;

  double sMin = 0;
  if (iBranchingType == 1 || iBranchingType == 11) //scalar DM
    sMin = 4 * mChi * mChi;

  if (ss < sMin)
    return 0.;   // A.C. e+e- -> Z' -> ff can happen also for an A' and chi with large mass,
                 // i.e. through the off-shell tail of the resonance, but this still needs to be kinematically allowed

  double sigma = 1;
  switch (iBranchingType) {
  case 0: //Lmu-Ltau, vanilla, neutrinos (nMu+nTau contribution). See 2206.03101 Eq. A.4.
    sigma = (4 * M_PI * alphaEW) / 3.;
    sigma = sigma * epsil * epsil / (4 * M_PI); //alphaZ' for Z' -> nunu
    sigma = sigma * epsil * epsil * piF2(ss); //|Pi(S)|^2 for e+e- --> Z'
    sigma = sigma * ss;
    break;
  case 1: //Lmu-Ltau, scalar DM.  See 2206.03101 Eq. A.3.
    sigma = (M_PI * alphaEW) / 3.;
    sigma = sigma * alphaD; //alphaD for Z' -> DM-DM
    sigma = sigma * epsil * epsil * piF2(ss); //|Pi(S)|^2 for e+e- --> Z'
    sigma = sigma * ss;
    sigma = sigma * pow(1-sMin/ss,3./2);
    break;
  case 10: //B-L, vanilla, neutrinos
    sigma = sigma * epsil * epsil / 3.; // for e+e- --> Z'
    sigma = sigma * epsil * epsil / (4 * M_PI); //alphaZ' for Z' -> nunu
    sigma = sigma * ss;
    break;
  case 11: //B-L, scalar DM
    break;
  default:
    break;
  }

  //here sigma is in G4 internal units, 1 /Energy^2. Move to pBarn;
  sigma *= GeVtoPb;

  //A.C. correct here for atomic effects
  sigma = sigma * ZNucl;
  return sigma;
}

double DarkZAnnihilation::GetSigmaTot(double E0) {
  return TotalCrossSectionCalc(E0);
}


//E0: positron total energy in GeV
bool DarkZAnnihilation::EmissionAllowed(double E0, double DensityMat) // Different kinematic limit here
{
  switch (iBranchingType){
  case 0: //Lmu-Ltau, vanilla, neutrino
  case 1:
    if (sqrt(2. * Mel * E0 + 2*Mel*Mel) < 2. * mChi)
      return false;
    break;
  case 10: //B-L, vanilla, neutrino
  case 11:
    if (sqrt(2. * Mel * E0 + 2*Mel*Mel) < 2. * mChi)
      return false;
    break;
  }

  if (E0 < EThresh)
    return false;
  if (NEmissions)
    return false; // For G4 DM classes
  if (fabs(DensityMat - Density) > 0.1)
    return false;
  return true;
}

double DarkZAnnihilation::CrossSectionDSDX(double XEv, double E0) {
  (void) (E0);
  if (XEv > 0.9999)
    return 1.;
  return 0.;
}

double DarkZAnnihilation::CrossSectionDSDXDU(double XEv, double UThetaEv, double E0) {
  (void) (UThetaEv);
  (void) (E0);
  if (XEv > 0.9999)
    return 1.;
  return 0.;
}

double DarkZAnnihilation::Width() {
  double ret=0;
  switch (iBranchingType) {
  case 0: //Lmu-Ltau, vanilla
      ret=MA/3*epsil*epsil/(4*M_PI); //Z->nu nu (mu+tau)
      if (MA > 2*Mmu)
        ret += MA/3 * epsil*epsil/(4*M_PI)*(1+2*Mmu*Mmu/(MA*MA))*sqrt(1-4*Mmu*Mmu/(MA*MA)); //Z->mu mu
    break;
  case 1: //Lmu-Ltau, DM
    ret=MA/3*epsil*epsil/(4*M_PI);//Z->nu nu
    if (MA > 2*Mmu)
      ret += MA/3 * epsil*epsil/(4*M_PI)*(1+2*Mmu*Mmu/(MA*MA))*sqrt(1-4*Mmu*Mmu/(MA*MA)); //Z->mu mu
    ret+= MA/12*alphaD*pow((1-4*mChi*mChi/(MA*MA)),3./2); //Z->DM DM
    break;
  case 10: //B-L, vanilla
      ret=MA*epsil*epsil/(4*M_PI); //Z->nu nu
      if (MA > 2*Mel)
        ret += MA/3 * epsil*epsil/(4*M_PI)*(1+2*Mel*Mel/(MA*MA))*sqrt(1-4*Mel*Mel/(MA*MA)); //Z->el el
      if (MA > 2*Mmu)
        ret += MA/3 * epsil*epsil/(4*M_PI)*(1+2*Mmu*Mmu/(MA*MA))*sqrt(1-4*Mmu*Mmu/(MA*MA)); //Z->mu mu
    break;
  case 11: //B-L, DM
    break;
  }

  return ret;
}

void DarkZAnnihilation::SetMA(double MAIn) {
  std::cout << "DarkZAnnihilation::SetMA was called with MAIn = " << MAIn << std::endl;

  // B.B: Not used at the moment!
  if (iBranchingType == 100) {
    mChi1 = MA * r;
    mChi2 = (1. + f) * mChi1;
  } else {
    mChi = MA * r;
    mChi1 = mChi;
    mChi2 = mChi;
  }

}

double DarkZAnnihilation::AngularDistributionResonant(double eta, double E0) {

  double ss = 2. * Mel * E0;
  double qq;
  double val = 0;
  switch (iBranchingType) {
  case 0: //Lmu-Ltau, vanilla, neutrino --> 1+eta*eta
    val=(1+eta*eta)/2; //Max: eta=1
    break;
  case 1: //Lmu-Ltau, DM, scalar
    //Scalar LDM, Angular distribution f(eta) ~ 1-eta*eta. Max: eta=0;
    val = 1-eta*eta;
    break;
  case 10: //B-L, vanilla, neutrino
    val = 1;
    break;
  case 11: //B-L, DM, scalar
    //Scalar LDM, Angular distribution f(eta) ~ 1-eta*eta. Max: eta=0;
    val = 1-eta*eta;
    break;
  }
  return val;
}

//This function is the result of the loop-factor |Pi(q^2)|^2 for the Lmu-Ltau model, computed for g=1.
//For a given value of g, the returned value must be multiplied by g*g
//See: https://gitlab.cern.ch/P348/DMG4/-/issues/17
double DarkZAnnihilation::piF2(double m_q2) { //q2 in GeV^

  double ee = sqrt(4 * M_PI * alphaEW);

  std::complex<double> one(1., 0.);
  std::complex<double> q2(m_q2, 0.);
  std::complex<double> rmu(Mmu * Mmu / m_q2, 0.);
  std::complex<double> rtau(Mtau * Mtau / m_q2, 0.);

  std::complex<double> ret = 0.5 * log(Mtau * Mtau / (Mmu * Mmu));
  ret = ret + 2. * (rmu - rtau);
  ret = ret - (one + 2. * rmu) * sqrt(one - 4. * rmu) * acoth(sqrt(one - 4. * rmu));
  ret = ret + (one + 2. * rtau) * sqrt(one - 4. * rtau) * acoth(sqrt(one - 4. * rtau));

  ret = ret * ee / (2 * M_PI * M_PI);

  ret = ret / 3.;

  return norm(ret);
}

