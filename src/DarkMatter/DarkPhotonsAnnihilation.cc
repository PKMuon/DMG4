/*
 * DarkPhotonsAnnihilation.cc
 *
 *  Created on: Oct 6, 2020
 *      Author: celentan
 *  Fixed: Nov 2, 2020
 */

#include "DarkMatterAnnihilation.hh"
#include "DarkPhotonsAnnihilation.hh"
#include "Utils.hh"


#include <iostream>
#include <cmath>


DarkPhotonsAnnihilation::DarkPhotonsAnnihilation(double MAIn, double EThreshIn, double SigmaNormIn, double ANuclIn, double ZNuclIn, double DensityIn,
                                                 double epsilIn, int IDecayIn, double rIn, double alphaDIn, int IBranchingIn, double fIn) :
                                                 DarkMatterAnnihilation(MAIn, EThreshIn, SigmaNormIn, ANuclIn, ZNuclIn, DensityIn, epsilIn, IDecayIn, rIn,alphaDIn, IBranchingIn,fIn)
{
  DMType = 1; //A.C.

  std::cout << "Initialized DarkPhotonsAnnihilation (e+ e- -> A' -> DM DM) for material density = " << DensityIn << std::endl;
  std::cout << "mA: "<<MA*1000<<" MeV "<<std::endl;
  std::cout << "IbranchingType: "<<iBranchingType<<std::endl;
  if (iBranchingType==2){
    std::cout<<"mChi1: "<<mChi1*1000<<" MeV "<<std::endl;
    std::cout<<"mChi2: "<<mChi2*1000<<" MeV "<<std::endl;
  }else{
    std::cout<<"mChi: "<<mChi*1000<<" MeV "<<std::endl;
  }
  std::cout<<"Width: "<<this->Width()*1E3<<" MeV "<<std::endl;




  std::cout << std::endl;
}


DarkPhotonsAnnihilation::~DarkPhotonsAnnihilation()
{;}


//Convenience private method to be shared among TotalCrossSectionCalc and GetSigmaMax.
//This is the total cross section without the BW denominator
//s: e+e- invariant mass squared
double DarkPhotonsAnnihilation::PreFactor(double ss){

  double qq=0.,E1=0.,E2=0.;

  if (ss<this->sMin()){
    return 0;
  }

  qq = this->q(ss);

  switch (iBranchingType) {
  case 2:
      E1=(ss-mChi2*mChi2+mChi1*mChi1)/(2*sqrt(ss));
      E2=(ss+mChi2*mChi2-mChi1*mChi1)/(2*sqrt(ss));
      break;
  }


  double sigma = 4 * M_PI * alphaEW * epsilBench * epsilBench * alphaD;
  sigma = sigma * qq / sqrt(ss);



  switch (iBranchingType){
  case 0:
      sigma = sigma * (ss - 4. / 3. * qq * qq); // A.C. this is for final state fermions (default)
      break;
  case 1:
      sigma = sigma * (2./3. * qq * qq);      // A.C. this is for final state scalars
      break;
  case 2:
      sigma = sigma *2*(mChi1 * mChi2 + E1*E2+qq*qq/3);                      //inelastic DM
      break;
  }



  //here sigma is in G4 internal units, 1 /Energy^2. Move to pBarn;
  sigma *= GeVtoPb;

  return sigma;
}

//E0: positron TOTAL energy in lab frame
double DarkPhotonsAnnihilation::GetSigmaTot(double E0) {
  return TotalCrossSectionCalc(E0);
}



double DarkPhotonsAnnihilation::CrossSectionDSDX(double XEv, double E0) {
  (void)(E0);
  if (XEv > 0.9999) return 1.;
  return 0.;
}


double DarkPhotonsAnnihilation::CrossSectionDSDXDU(double XEv, double UThetaEv, double E0) {
  (void)(UThetaEv);
  (void)(E0);
  if (XEv > 0.9999) return 1.;
  return 0.;
}


double DarkPhotonsAnnihilation::Width() {
  double ret;
  ret = MA * epsil * epsil * alphaEW * 1. / 3;

  switch (iBranchingType){
  case 0:
      if (MA/2. > mChi)
      ret += MA * alphaD * (1 + 2 * mChi * mChi / (MA * MA)) * sqrt(1 - 4 * mChi * mChi / (MA * MA)) / 3.;
      break;
  case 1:
      if (MA/2. > mChi){
          ret += MA *alphaD/12 * pow(1-4*mChi*mChi/(MA*MA),3./2);
      }
      break;
  case 2:
      if (MA > (mChi1+mChi2)){
          double E1=(MA*MA-mChi2*mChi2+mChi1*mChi1)/(2*MA);
          double qq=sqrt(E1*E1-mChi1*mChi1);
          ret+= alphaD*2*qq/(3*MA*MA)*(3*mChi1*mChi2+MA*MA-mChi1*mChi1/2-mChi2*mChi2/2-pow(mChi1*mChi1-mChi2*mChi2,2)/(2*MA*MA));
      }
      break;
  }

  return ret;
}


void DarkPhotonsAnnihilation::SetMA(double MAIn) {
  std::cout << "DarkPhotonsAnnihilation::SetMA was called with MAIn = " << MAIn << std::endl;



  if (iBranchingType==2){

       mChi1 = MA * r;
       mChi2 = (1. + f) * mChi1;
   }else{

       mChi = MA * r;
       mChi1=mChi;
       mChi2=mChi;
   }


}

//E0: positron TOTAL energy in lab frame
double DarkPhotonsAnnihilation::AngularDistributionResonant(double eta,double E0){


    double ss = 2. * Mel * E0+2*Mel*Mel;

    //TODO for atomic effects
    if (ss<this->sMin()){
      printf("DarkPhotonsAnnihilation::AngularDistribution error with threshold, E0=%f, m=%f\n",E0,mChi);
      exit(1);
    }

    double qq=this->q(ss);
    double val=0;
    switch (iBranchingType){
     case 0:
         //Fermionic LDM. Angular distribution f(eta) ~ s+4*qq*qq*eta*eta+4*m*m. Max for eta=+-1
         val=ss+4*qq*qq*eta*eta+4*mChi*mChi;
         val/=(ss+4*qq*qq+4*mChi*mChi);
         break;
     case 1:
         //Scalar LDM, Angular distribution f(eta) ~ 1-eta*eta. Max: eta=0;
         val=1-eta*eta;
         break;
     case 2:
         //Asymmetric LDM. Angular distribution f(eta)=(m1*m2+E1*E2+qq*qq*eta*eta). Max for eta=+-1
         double E1=(ss+mChi1*mChi1-mChi2*mChi2)/(2*sqrt(ss));
         double E2=(ss-mChi1*mChi1+mChi2*mChi2)/(2*sqrt(ss));


         val=(mChi1*mChi2+E1*E2+qq*qq*eta*eta);
         val/=(mChi1*mChi2+E1*E2+qq*qq);
         break;
    }

    return val;
}

