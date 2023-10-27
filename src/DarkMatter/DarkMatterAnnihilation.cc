// This is a base class to simulate resonant annihilation production of dark photons and dark scalars A 
// To be used in a Geant4 application.
//
// Concrete implementations: DarkPhotonsAnnihilation, DarkScalarsAnnihilation, ...
//
#include "DarkMatter.hh"
#include "DarkMatterAnnihilation.hh"
#include "Randomize.hh"
#include <iostream>
#include <map>
#include <vector>

DarkMatterAnnihilation::DarkMatterAnnihilation(double MAIn, double EThreshIn, double SigmaNormIn, double ANuclIn, double ZNuclIn, double DensityIn,
    double epsilIn, int IDecayIn, double rIn, double alphaDIn, int IBranchingIn, double fIn) :
    DarkMatter(MAIn, EThreshIn, SigmaNormIn, ANuclIn, ZNuclIn, DensityIn, epsilIn, IDecayIn), iBranchingType(IBranchingIn), r(rIn), f(fIn), alphaD(alphaDIn) {
  // Default values, to be redefined in concrete implementations if necessary
  DMType = 1; //A.C.
  ParentPDGID = -11;
  DaughterPDGID = 11;
  mChi = MAIn / 3;
  mChi1 = mChi;
  mChi2 = mChi;

  if (iBranchingType == 2) {
    mChi1 = MA * r;
    mChi2 = (1. + f) * mChi1;
  } else {
    mChi = MA * r;
    mChi1 = mChi;
    mChi2 = mChi;
  }

  deltaMchi = mChi2 - mChi1;

  //std::cout << "Initialized DarkMatterAnnihilation (e+ e- -> A' -> DM DM) for material density = " << DensityIn << std::endl;
}

DarkMatterAnnihilation::~DarkMatterAnnihilation() {
  ;
}

/*This method returns a random cosine for e+e- --> A' --> ff in the CM frame
 * E0: positron total energy in LAB frame (GeV units)
 */
double DarkMatterAnnihilation::SimulateEmissionResonant(double E0) {

  int maxiter = 25000000;
  double eta;
  double fcomp, frand;
  for (int iii = 1; iii < maxiter; iii++) {
    eta = G4UniformRand() * 2 - 1; //between -1 and 1;
    fcomp = this->AngularDistributionResonant(eta, E0);
    frand = G4UniformRand();
    if (frand < fcomp)
      return eta;
  }
  return 1;
}

//Input: E0, positron total energy in GeV
double DarkMatterAnnihilation::BreitWignerDenominator(double E0) {
  double ss = 2. * Mel * E0 + 2 * Mel * Mel;
  double gg = this->Width();
  double BWden = ((ss - MA * MA) * (ss - MA * MA) + MA * MA * gg * gg);
  return BWden;
}

//Input: E0, positron total energy in GeV
//output: total annihilation cross-section in pbarn.
//Since the framework assumes this method is returning the total cross section per nucleous, for the moment I scale this by Z.
double DarkMatterAnnihilation::TotalCrossSectionCalc(double E0) {
  double sigma = this->PreFactor(E0);
  sigma = sigma / this->BreitWignerDenominator(E0);
  //A.C. correct here for atomic effects
  sigma = sigma * ZNucl;
  return sigma;
}

//A.C. useful function to directly return the maximum cross section value (for numerical precision)
double DarkMatterAnnihilation::GetTotalCrossSectionMax() {
  double Eres = (MA * MA - 2 * Mel * Mel) / (2. * Mel);
  double sigma = this->PreFactor(Eres);
  double gg = this->Width();
  sigma = sigma / (MA * MA * gg * gg);
  //A.C. correct here for atomic effects
  sigma = sigma * ZNucl;
  return sigma;
}

/*A.C. this function returns the total cross section for the annihilation process, accounting for atomic effects.
 *The electrons in a given atomic shell are assumed to be moving around the nuclei with spherical symmetry.
 *The electrons in a given atomic shell are assumed to have energies distributed according to a certain distribution, handled in DMProcessAnnihilation.cc
 *See: https://gitlab.cern.ch/P348/DMG4/-/issues/14?work_item_iid=21

 INPUT:

 E0:  positron total energy in GeV
 Z: for each shell (map index), a vector with the number of electrons in that shell
 ene: for each shell (map index), a vector with the kinetic energies of the electrons in that shell, that were computed by MC method, in GeV. Use this to integrate the energy dependency
 */

double DarkMatterAnnihilation::GetSigmaTotAtomicEffects(double E0, std::map<int, int> &Z, std::map<int, std::vector<double>> &ene) {
  double sigma = 0;
  for (int is = 0; is < Z.size(); is++) {
    int ZeleShell = Z[is];
    std::vector<double> eneShell = ene[is];
    double sigmaShell = GetSigmaTotAtomicEffectsOneShell(E0, ZeleShell, eneShell);
    sigma = sigma + sigmaShell;
  }
  return sigma;
}

double DarkMatterAnnihilation::GetSigmaTotAtomicEffectsOneShell(double E0, int &Zshell, std::vector<double> &eneShell) {
  double Ep = E0;
  double Pp = sqrt(Ep * Ep - Mel * Mel);
  double W = this->Width();
  double pre = this->PreFactor(E0);

  //MC-integration of the energy dependency
  double sigmaShell = 0;
  for (double Tm : eneShell) {
    double Em = Tm + Mel;
    double Pm = sqrt(Em * Em - Mel * Mel);

    double arg1 = (2 * Pp * Pm - 2 * Mel * Mel - 2 * Ep * Em + MA * MA) / (MA * W);
    double arg2 = (-2 * Pp * Pm - 2 * Mel * Mel - 2 * Ep * Em + MA * MA) / (MA * W);

    sigmaShell = sigmaShell + 1. / (4 * MA * W * Pp * Pm) * (atan(arg1) - atan(arg2));
  }
  sigmaShell /= eneShell.size();

  //Now add the common factors
  sigmaShell = sigmaShell * Zshell * pre;

  return sigmaShell;

}

double DarkMatterAnnihilation::GetTotalCrossSectionMaxAtomicEffects(std::map<int, int> &Z, std::map<int, std::vector<double>> &ene) {
  double Eres = (MA * MA - 2 * Mel * Mel) / (2. * Mel);
  return GetSigmaTotAtomicEffects(Eres, Z, ene);
}
