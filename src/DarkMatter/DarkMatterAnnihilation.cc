// This is a base class to simulate resonant annihilation production of dark photons and dark scalars A 
// To be used in a Geant4 application.
//
// Concrete implementations: DarkPhotonsAnnihilation, DarkScalarsAnnihilation, ...
//
#include "DarkMatter.hh"
#include "DarkMatterAnnihilation.hh"
#include "Randomize.hh"
#include <iostream>

DarkMatterAnnihilation::DarkMatterAnnihilation(double MAIn, double EThreshIn, double SigmaNormIn, double ANuclIn, double ZNuclIn, double DensityIn,
                                               double epsilIn, int IDecayIn, double rIn, double alphaDIn, int IBranchingIn, double fIn)
                                                : DarkMatter(MAIn, EThreshIn, SigmaNormIn, ANuclIn, ZNuclIn, DensityIn, epsilIn, IDecayIn),
                                                iBranchingType(IBranchingIn), r(rIn), f(fIn), alphaD(alphaDIn)
{
  // Default values, to be redefined in concrete implementations if necessary
  DMType = 1; //A.C.
  ParentPDGID = -11;
  DaughterPDGID = 11;
  mChi = MAIn/3;
  mChi1=mChi;
  mChi2=mChi;

  if (iBranchingType==2){
      mChi1 = MA * r;
      mChi2 = (1. + f) * mChi1;
  }else{
      mChi = MA * r;
      mChi1=mChi;
      mChi2=mChi;
  }

  deltaMchi=mChi2-mChi1;

  //std::cout << "Initialized DarkMatterAnnihilation (e+ e- -> A' -> DM DM) for material density = " << DensityIn << std::endl;
}

DarkMatterAnnihilation::~DarkMatterAnnihilation()
{;}


/*This method returns a random cosine for e+e- --> A' --> ff in the CM frame
 * E0: positron energy in LAB frame (GeV units)
 */
double DarkMatterAnnihilation::SimulateEmissionResonant(double E0){

    int maxiter = 25000000;
    double eta;
    double fcomp,frand;
    for( int iii = 1; iii < maxiter; iii++) {
        eta = G4UniformRand()*2 -1; //between -1 and 1;
        fcomp=this->AngularDistributionResonant(eta,E0);
        frand=G4UniformRand();
        if (frand<fcomp) return eta;
     }
    return 1;
}

/*This method returns the angular distribution for e+e- --> A' --> ff in the CM frame
 * It has to be implemented in the derived classes (the default method is a dummy implementation)
 * eta: cosine of the f in the CM frame
 * E0: positron beam energy in LAB frame (GeV units)
 * IMPORTANT: it has to be normalized so that the maximum is one.
 */
double DarkMatterAnnihilation::AngularDistributionResonant(double eta,double E0){
    return (1+eta*eta)/2.;
}

//Input: E0, positron energy in GeV
//output: total annihilation cross-section in pbarn.
//Since the framework assumes this method is returning the total cross section per nucleous, for the moment I scale this by Z.
double DarkMatterAnnihilation::TotalCrossSectionCalc(double E0) {

    double ss = 2. * Mel * E0;
    double sigma=this->PreFactor(E0);
    double gg = this->Width();
    sigma=sigma/((ss - MA * MA) * (ss - MA * MA) + MA * MA * gg * gg);
    return sigma;
}

//A.C. useful function to directly return the maximum cross section value (for numerical precision)
double DarkMatterAnnihilation::GetTotalCrossSectionMax(){

  double Eres = (MA*MA/(2.*Mel));
  double sigma=this->PreFactor(Eres);
  double gg = this->Width();
  sigma=sigma/(MA * MA * gg * gg);
  return sigma;
}


