#include "DarkMatterAnnihilation.hh"
#include "DarkMassSpin2Annihilation.hh"
#include "Utils.hh"

#include <iostream>
#include <cmath>


DarkMassSpin2Annihilation::DarkMassSpin2Annihilation( double MAIn
                                                    , double EThreshIn
                                                    , double SigmaNormIn
                                                    , double ANuclIn
                                                    , double ZNuclIn
                                                    , double DensityIn
                                                    , double epsilIn
                                                    , int IDecayIn
                                                    , double rIn
                                                    , double alphaDIn
                                                    , int IBranchingIn
                                                    , double fIn ) 
                                                    :DarkMatterAnnihilation(MAIn, EThreshIn, SigmaNormIn, ANuclIn, ZNuclIn, DensityIn, epsilIn, IDecayIn, rIn,alphaDIn, IBranchingIn,fIn)
{
  DMType = 5;

  mChi = MA * r;
  mChi2 = mChi;
  mChi1 = mChi;

  std::cout << "Initialized DarkMassSpin2Annihilation " 
            << "(e+ e- -> A' -> DM DM) for material density = "
            << DensityIn << std::endl << std::endl;
}


DarkMassSpin2Annihilation::~DarkMassSpin2Annihilation(  ){ ; }


// Total cross-section in pBarn
//s: e+e- invariant mass squared
double DarkMassSpin2Annihilation::PreFactor( double ss )
{
  double  rsChi = mChi*mChi / ss, rsEl = Mel*Mel / ss;
  // A.C. e+e- -> G -> chi chi can happen also for an G and chi with large 
  // mass, i.e. through the off-shell tail of the resonance, but this still 
  // needs to be kinematically allowed.
  if ( sqrt(ss) < 2.0 * mChi && sqrt(ss) < 2.0 * Mel ){ return 0.0; }
  double sigma =   ( 1.0/(256.0*M_PI) ) * ss*ss*ss 
  //               * sqrt( 1.0 - 4.0 * rsEl ) * (1.0 + (8.0/3.0) * rsEl )
                 * pow( 1.0 - 4.0*rsChi, 3.0/2.0 ) * (1.0 + (8.0/3.0) * rsChi);
  
  // here sigma is in  1 /Energy^2. Move to pBarn
  sigma = sigma * GeVtoPb * 4*M_PI * alphaD * epsilBench*epsilBench;

  return sigma;
}

//E0: positron TOTAL energy in lab frame
double DarkMassSpin2Annihilation::GetSigmaTot( double E0 )
{
  return TotalCrossSectionCalc( E0 );
}


// Different kinematic limit here
//E0: positron TOTAL energy in lab frame
bool DarkMassSpin2Annihilation::EmissionAllowed( double E0, double DensityMat )
{
  if ( sqrt(2. * Mel * E0+2*Mel*Mel) < 2.0 * mChi ){ return false; }
  if ( E0 < EThresh ){ return false; }
  if ( NEmissions ){ return false; }   // For G4 DM classes
  if ( fabs(DensityMat - Density) > 0.1 ){ return false; }
  return true;
}


double DarkMassSpin2Annihilation::CrossSectionDSDX( double XEv, double E0 )
{
  (void)(E0);
  if ( XEv > 0.9999 ){ return 1.; }   
  return 0.;
}


double DarkMassSpin2Annihilation::CrossSectionDSDXDU( double XEv
                                                    , double UThetaEv
                                                    , double E0 )
{
  (void)(UThetaEv); (void)(E0);
  if ( XEv > 0.9999 ){ return 1.; }   
  return 0.;
}


// decay width of decay into a pair of dark fermions
double DarkMassSpin2Annihilation::Width()
{
  double rMAChi = mChi*mChi / (MA*MA);
  return  ( 1.0/(160.0*M_PI) ) * MA*MA*MA * 4.0*M_PI * alphaD 
        * pow( 1.0 - 4.0 * rMAChi, 3.0/2.0 ) * ( 1.0 + ( 8.0 / 3.0 ) * rMAChi );
}

void DarkMassSpin2Annihilation::SetMA(double MAIn) {
    std::cout << "DarkMassSpin2Annihilation::SetMA was called with MAIn = " << MAIn << std::endl;
    mChi = MA * r;
    mChi1 = mChi;
    mChi2 = mChi;
}

//E0: positron TOTAL energy in lab frame
double DarkMassSpin2Annihilation::AngularDistributionResonant(double E0,double eta){
  return 1.;
}
