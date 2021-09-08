// This is a class to simulate dark scalars production by electrons in matter eN -> eNA
// Description is in the base class DarkMatter code
// To be used in a Geant4 application.
//
//
#include "DarkMatter.hh"
#include "DarkScalars.hh"
#include "Utils.hh"

#include "Randomize.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>

#include <iostream>
#include "G4ios.hh"


#include "KFactorsScalars.code"


DarkScalars::DarkScalars(double MAIn, double EThreshIn, double SigmaNormIn, double ANuclIn, double ZNuclIn, double DensityIn,
                         double epsilIn, int IDecayIn)
: DarkMatter(MAIn, EThreshIn, SigmaNormIn, ANuclIn, ZNuclIn, DensityIn, epsilIn, IDecayIn)
{
  DMType = 2;
  ParentPDGID = 11;
  DaughterPDGID = 11;
  std::cout << "Initialized DarkScalars off electrons and positrons for material density = " << DensityIn << std::endl;
  std::cout << std::endl;
}


DarkScalars::~DarkScalars()
{;}


double DarkScalars::TotalCrossSectionCalc(double E0)
{
  //double ThetaMaxA;
  //double ThetaMaxEl;
  double sigmaTot;

  if(E0 < 2.*MA) return 0.;

  if(MA > 0.001) { // analytical calculation above 1 MeV

    double tmin = MA*MA*MA*MA/(4.*E0*E0);
    double tmax = MA*MA;

    double aa = 111.*pow(ZNucl,-1./3.)/Mel;
    double d = 0.164*pow(ANucl,-2./3.);
    double ta = pow(1./aa,2.);
    double td = d;
    double fluxAnalytical = ZNucl*ZNucl*(-((td*td*(((ta - td)*(ta + td + 2.0*tmax)*(tmax - tmin))/((ta + tmax)*(td + tmax)) + (ta + td + 2.0*tmin)*(log(ta + tmax)
                          - log(td + tmax) - log(ta + tmin) + log(td + tmin))))/((ta-td)*(ta-td)*(ta-td))));

    double beta = sqrt(1. - MA*MA/(E0*E0));
    double cutoff1 = Mel/MA;
    double cutoff2 = MA/E0;
    double cutoff = cutoff2;
    if(cutoff1 > cutoff2) cutoff = cutoff1;

    sigmaTot= GeVtoPb*(2./3.)*alphaEW*alphaEW*alphaEW*epsilBench*epsilBench*fluxAnalytical*beta*log(1./(cutoff*cutoff))/(MA*MA);
    if(sigmaTot < 0.) sigmaTot=0.;

    double KFactor = KfactorScalarsApproximate(MA, E0);

    G4cout << "Total CS scalars calc, E = " << E0 << "  M = " << MA << "  KFactor = " << KFactor << G4endl;

    return sigmaTot / KFactor; // This K-factor decreases the cross section for MA > ~2 MeV

  } else {

    G4cout << "No analytical calculations below 1 MeV, exiting" << G4endl;
    exit(1);
  }
}


double DarkScalars::GetSigmaTot(double E0)
{
  if(MA > 0.001) {
    return GetSigmaTot0(E0);
  } else {
    //return TotCSScalar(MA);
    G4cout << "CS for masses below 1 MeV is not yet implemented" << G4endl;
    exit(1);
  }
}


double DarkScalars::CrossSectionDSDX(double XEv, double E0)
{
  if(MA > 0.001) {
    double momentumOfDP=sqrt(XEv*XEv*E0*E0-MA*MA);
    double umaxtilde = -MA*MA*(1.0-XEv)/XEv - Mel*Mel*XEv;
    double Numerator = Mel*Mel*(2. -XEv)*(2. -XEv) - 2.*umaxtilde*XEv;
    double Denominator = 3.*umaxtilde*umaxtilde;
    double sigma = momentumOfDP*Numerator/Denominator;
    return sigma;
  } else {
    std::cout << "DarkScalars: Error: differential cs for scalars below 1 MeV is not implemented, exiting" << std::endl;
    exit(1);
  }
}


double DarkScalars::CrossSectionDSDXDU(double XEv, double UThetaEv, double E0)
{
  if(MA > 0.001) {
    double Uxtheta = 2.*E0*E0*UThetaEv*XEv + MA*MA*(1. - XEv)/XEv + Mel*Mel*XEv;
    double AA = XEv*XEv / (2*Uxtheta*Uxtheta);
    double BB = (1. - XEv)*(1. - XEv)*(MA*MA - 4.0*Mel*Mel)/(Uxtheta*Uxtheta*Uxtheta*Uxtheta);
    double CC = MA*MA - Uxtheta*XEv/(1. - XEv) + Mel*Mel*XEv*XEv/(1. - XEv);
    double sigma = sqrt(XEv*XEv - MA*MA/(E0*E0)) * (AA + BB*CC);
    return sigma;
  } else {
    std::cout << "DarkScalars: Error: differential cs for scalars below 1 MeV is not implemented, exiting" << std::endl;
    exit(1);
  }
}


double DarkScalars::Width()
{
  return 1./2.*1./137.*MA*epsil*epsil*sqrt(1.-4.*Mel*Mel/(MA*MA))*(1.-4.*Mel*Mel/(MA*MA));
}
