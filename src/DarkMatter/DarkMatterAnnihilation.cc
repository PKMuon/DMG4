// This is a base class to simulate resonant annihilation production of dark photons and dark scalars A 
// To be used in a Geant4 application.
//
// Concrete implementations: DarkPhotonsAnnihilation, DarkScalarsAnnihilation, ...
//
#include "DarkMatterAnnihilation.hh"
#include "Randomize.hh"
#include <iostream>
#include <map>
#include <vector>
#include <cmath>

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

bool DarkMatterAnnihilation::EmissionAllowed(double E0, double DensityMat){
  if (E0 < EThresh)
     return false;
   if (NEmissions)
     return false; // For G4 DM classes
   if (fabs(DensityMat - Density) > 0.1)
     return false;
   return true;

}

double DarkMatterAnnihilation::q(double s){
  double qq=0;
  if (iBranchingType == 2){
      qq=s/4+(pow(mChi1,4)+pow(mChi2,4)-2*pow(mChi1,2)*pow(mChi1,2))/s-(pow(mChi1,2)+pow(mChi2,2))/2;
      return sqrt(qq);
   }
   else{
     qq=sqrt(s/4.-mChi*mChi);
   }
   return qq;
}

double DarkMatterAnnihilation::sMin(){
  double smin=0;
  if (iBranchingType == 2){
    smin=(mChi1+mChi2)*(mChi1+mChi2);
  }
  else{
    smin=(4*mChi*mChi);
  }
  return smin;
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

//Input: s, e+e- invariant mass squared GeV^2
double DarkMatterAnnihilation::BreitWignerDenominator(double ss) {
  double gg = this->Width();
  double BWden = ((ss - MA * MA) * (ss - MA * MA) + MA * MA * gg * gg);
  return BWden;
}

//Input: E0, positron total energy in GeV
//output: total annihilation cross-section in pbarn.
//Since the framework assumes this method is returning the total cross section per nucleous, for the moment I scale this by Z.
double DarkMatterAnnihilation::TotalCrossSectionCalc(double E0) {

  double ss = 2. * Mel * E0 + 2 * Mel * Mel;


  double sigma = this->PreFactor(ss);
  sigma = sigma / this->BreitWignerDenominator(ss);
  //A.C. correct here for atomic effects
  sigma = sigma * ZNucl;
  return sigma;
}

//A.C. useful function to directly return the maximum cross section value (for numerical precision)
double DarkMatterAnnihilation::GetTotalCrossSectionMax() {
  double sigma = this->PreFactor(MA*MA);
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

double DarkMatterAnnihilation::GetSigmaTotAtomicEffects(double E0,const std::map<int, int>& Z,const std::map<int, std::vector<double>>& ene) {
  double sigma = 0;
  for (int is = 0; is < Z.size(); is++) {
    int ZeleShell = Z.at(is);
    const std::vector<double>& eneShell = ene.at(is);

    //double sigmaShell = GetSigmaTotAtomicEffectsOneShell(E0, ZeleShell, eneShell);
    double sigmaShell = GetSigmaTotAtomicEffectsOneShellFull(E0, ZeleShell, eneShell);

    sigma = sigma + sigmaShell;
  }

  return sigma;
}




//This implementation considers the full dependency
double DarkMatterAnnihilation::GetSigmaTotAtomicEffectsOneShellFull(double E0,int Zshell,const std::vector<double> &eneShell) {
  double Ep = E0;
  double Pp = sqrt(Ep * Ep - Mel * Mel);
  double W = this->Width();


  //MC-integration of the energy dependency
  double sigmaShell = 0;
  for (double Tm : eneShell) {
    double Em = Tm + Mel;
    double Pm = sqrt(Em * Em - Mel * Mel);

    double smin=this->sMin();
    double zmin=-1;
    double zmax=(2*Mel*Mel+2*Ep*Em-smin)/(2*Ep*Pm);
    if (zmax>1){
      zmax=1;
    }
    if (zmax<-1){ //the requirement regarding smin is not satisfied for any value of z. Therefore, this electron cannot contribute at all (but has to be accounted for in the average)
      sigmaShell = sigmaShell+0;
      continue;
    }

    //integration limits
    long double s_zmax = 2 * Mel * Mel + 2 * Ep * (Em - zmax * Pm);
    long double s_zmin = 2 * Mel * Mel + 2 * Ep * (Em - zmin * Pm);
    long double x_zmax = s_zmax - 4 * mChi * mChi;
    if (x_zmax < 0)
      x_zmax = 0; //numerical precision
    long double x_zmin = s_zmin - 4 * mChi * mChi;
    if (x_zmin < 0)
      x_zmin = 0;  //numerical precision

    x_zmax=sqrt(x_zmax);
    x_zmin=sqrt(x_zmin);




    //middle value of z in the allowed z-integral region (if all values are allowed, z=0)
    //    double zmid=(zmax+zmin)/2;
    //     double s_zmid=2*Mel*Mel+2*Ep*(Em-zmid*Pm);
    const int nz=20;
    const double dz=(zmax-zmin)/nz;
    double pre=0;
    for (int iz=0;iz<nz;iz++){
      double zThis=zmin+dz*iz+dz/2;
      double sThis=2*Mel*Mel+2*Ep*(Em-zThis*Pm);
      double preThis=this->PreFactor(sThis)/this->q(sThis);
      if ((preThis<0)||(!isfinite(preThis))){
        preThis=0; //numerical precision
      }
      pre=pre+preThis;
    }
    pre=pre/(nz);
    pre=pre/(8*Ep*Pm); //jacobian



    long double delta2=MA*MA-4*mChi*mChi;
    long double eta2=sqrt(delta2*delta2+MA*MA*W*W);

    double arg1=this->funIntegral(x_zmin,delta2,eta2);
    double arg2=this->funIntegral(x_zmax,delta2,eta2);

    double thisSigmaShell = pre*(arg1-arg2);

    if ((thisSigmaShell<0)||(!isfinite(thisSigmaShell))){
      thisSigmaShell=0; //work-around for numerical precision
    }


    sigmaShell = sigmaShell + thisSigmaShell;

  }
  sigmaShell /= eneShell.size();

  //Now add the common factors
  sigmaShell = sigmaShell * Zshell;

  return sigmaShell;



}

//This implementation considers only the dependency in the BW denominator
double DarkMatterAnnihilation::GetSigmaTotAtomicEffectsOneShell(double E0,int Zshell,const std::vector<double> &eneShell) {
  double Ep = E0;
  double Pp = sqrt(Ep * Ep - Mel * Mel);
  double W = this->Width();
  double ss = 2*Mel*Mel+2*Ep*Mel;
  double pre = this->PreFactor(ss);

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

double DarkMatterAnnihilation::GetTotalCrossSectionMaxAtomicEffects(const std::map<int, int> &Z,const std::map<int, std::vector<double>> &ene) {
  double Eres = (MA * MA - 2 * Mel * Mel) / (2. * Mel);
  return GetSigmaTotAtomicEffects(Eres, Z, ene);
}


//See: https://gitlab.cern.ch/P348/DMG4/-/issues/27 -> this is "P(x)" without B(x)/(8*Ep*Pm)
double DarkMatterAnnihilation::funIntegral(long double x, long double delta2, long double eta2){


  long double W = this->Width();
  long double p1=1/(2*sqrt(2*eta2+2*delta2)) * log((x*x-sqrt(2*eta2+2*delta2)*x+eta2)/(x*x+sqrt(2*eta2+2*delta2)*x+eta2));
  long double p2=sqrt(eta2+delta2)/(sqrt(2)*MA*W)*(atan((2*x*(sqrt(2*eta2+2*delta2))-2*eta2-2*delta2)/(2*MA*W))+atan((2*x*(sqrt(2*eta2+2*delta2))+2*eta2+2*delta2)/(2*MA*W)));

  double ret=(double)(p1+p2);


  return ret;

}
