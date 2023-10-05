 // This is a class to simulate flavour-changing dark Long-lived Phi boson production by muons in matter mu N -> tau N phi
// Description is in the base class DarkMatter code
// To be used in a Geant4 application.
//
//
#include "DarkMatter.hh"
#include "DarkLFCScalars.hh"
#include "Utils.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_dilog.h>

#include <iostream>
#include <cmath>

// Auxiliary structures and functions:


// Additional structure holding E0 value to be passed into GSL callbacks.
struct BoundParms {
    DarkLFCScalars * this_;
    double E0;
};


// A callback wrapping function for DarkLFCScalars::CrossSectionDSDXDTheta()
static double _DarkLFCScalarsDsDxDThetaMuon(double x[], size_t dim, void * parms_) {
    (void)dim; // to avoid warning
    BoundParms * parms = reinterpret_cast<BoundParms*>(parms_);
    // Forward invocation to target method
    return parms->this_->CrossSectionDSDXDTheta( x[0], x[1], parms->E0 );
}


// Class methods:  ---------------
DarkLFCScalars::DarkLFCScalars(double MAIn, double EThreshIn, double SigmaNormIn, double ANuclIn, double ZNuclIn, double DensityIn,
                         double epsilIn, int IDecayIn)
: DarkMatter(MAIn, EThreshIn, SigmaNormIn, ANuclIn, ZNuclIn, DensityIn, epsilIn, IDecayIn)
{
  DMType = 32; // scalar
  ParentPDGID = 13;
  DaughterPDGID = 0;
  epsilBench = 0.003; // From https://doi.org/10.1140/epjc/s10052-023-11891-3

  tMax =       10000.;  // tmax initial; Value 10000. means that tmax = E0*E0 will be taken
  ThetaMax =     0.3;   // Max. angle of Z

  std::cout << "Initialized Dark Z boson for material density = " << DensityIn << std::endl;
  std::cout << "Energy cutoff = " << EThresh << " GeV" << std::endl;
  std::cout << std::endl;
}


DarkLFCScalars::~DarkLFCScalars()
{;}


double DarkLFCScalars::TotalCrossSectionCalc(double E0)
{
  return TotalCrossSectionCalc_WW2(E0); // Integral of ds/dxdTheta
}


// Below is the 2 - dimensional integration of WW ds/dxdTheta
//
double DarkLFCScalars::TotalCrossSectionCalc_WW2(double E0)
{
  if(E0 < 2.*MA) return 0.;

  double Xmin1 = MA/E0;
  if(EThresh/E0 > Xmin1) Xmin1 = EThresh/E0;
  double Xmax1 = 1. - Mtau/E0;
  if(Xmax1 < Xmin1) return 0.;

  double PrefactorAlphaEW = alphaEW*alphaEW;

  double xl[2] = { Xmin1, 0.};
  double xu[2] = { Xmax1, ThetaMax};

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_monte_function G;
  BoundParms parms = {this, E0};
  G.f = _DarkLFCScalarsDsDxDThetaMuon;
  G.dim = 2;
  G.params = &parms;

  gsl_rng_env_setup();

  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  double res, err, sigmaTot;

  // monte_miser: adaptive MC integration

  size_t calls = 5000000;
  gsl_monte_miser_state* stat = gsl_monte_miser_alloc(2);
  gsl_monte_miser_integrate(&G, xl, xu, 2, calls, r, stat, &res, &err);
  gsl_monte_miser_free(stat);
  sigmaTot = GeVtoPb*res*PrefactorAlphaEW;

  // Other vegas integration methods:
  // gsl_monte_plain : plain MC
  // gsl_monte_vegas

#if 0
  gsl_monte_vegas_state* stat = gsl_monte_vegas_alloc(2);
  gsl_monte_vegas_integrate(&G, xl, xu, 2, 10000, r, s, &res, &err);
  do {
    gsl_monte_vegas_integrate (&G, xl, xu, 2, calls/5, r, s, &res, &err);
  } while (fabs (gsl_monte_vegas_chisq(stat) - 1.0) > 0.5);
  double sigmaTot = GeVtoPb*res*PrefactorAlphaEW;
  gsl_monte_vegas_free(stat);
#endif

  gsl_rng_free (r);
  std::cout << "Total CS calc, E = " << E0 << "  M = " << MA << " CS = " << sigmaTot << std::endl;
  return sigmaTot;
}



double DarkLFCScalars::GetSigmaTot(double E0)
{
  return GetSigmaTot0(E0);
}


double DarkLFCScalars::CrossSectionDSDX(double XEv, double E0)
{
  std::cout << "DarkLFCScalars: wrong value of IApprox, exiting" << std::endl;
  exit(1);
}

// Below is the IWW formula for the double differential cross-section
// from 1705.01633, see e.g. corresponding second line of Eq.(25)
//
double DarkLFCScalars::CrossSectionDSDXDU(double XEv, double UThetaEv, double E0)
{
  if(XEv*E0 <= MA) return 0.;
  double Uxtheta = E0*E0*2.*UThetaEv*XEv + MA*MA*(1.0-XEv)/XEv + Mmu*Mmu*XEv;
  double AA = (1. - XEv + XEv*XEv/2.) / (Uxtheta*Uxtheta);
  double BB = (1. - XEv)*(1. - XEv)*(MA*MA + 2.0*Mmu*Mmu)/(Uxtheta*Uxtheta*Uxtheta*Uxtheta);
  double CC = MA*MA - Uxtheta*XEv/(1. - XEv);
  double sigma = sqrt(XEv*XEv - MA*MA/(E0*E0)) * (AA + BB*CC);
  return sigma;
}


// WW cross section for the total cross section
//
double DarkLFCScalars::CrossSectionDSDXDTheta(double XEv, double ThetaEv, double E0)
{

  if(E0*XEv < EThresh) return 0.;
  double x2=XEv*XEv;
  double theta2=ThetaEv*ThetaEv;
  double d = 0.164*pow(ANucl,-2./3);
  double MA2= MA*MA;
  double Mmu2= Mmu*Mmu;
  double Mtau2 = Mtau*Mtau;
  double E02= E0*E0;
  double u = -XEv*E02*theta2 - MA2*(1.0-XEv)/XEv + Mmu2*(1.0-XEv);
  double s = Mmu2 - (u-Mtau2)/(1.0-XEv);
  double gV = epsil;
  double gA = 0.;
  double tmax = MA2;
  const double pi = 3.141592654;
  double tmin= pow((u-Mtau2), 2.0)/(4.0*E02*(1.0-XEv)*(1.0-XEv));
  // I've calculated ChiWWAnalytical by using mathematica's "Integrate[...]" function
  // and converted the resulted expression to C-like form
  double ChiWWAnalytical = ZNucl*ZNucl*( (d*tmin/tmax + d*(d+tmin)/(d+tmax) + (d+2*tmin)*log(tmax) - (d+2*tmin)*log(d+tmax))
                                       - (2*d + (d+2*tmin)*log(tmin) - (d+2*tmin)*log(d+tmin)) )/d;
 
  // Trace of the amplitude from Eq. A.36 in https://arxiv.org/abs/2211.00664
  double Factor1= 4*Mmu*Mtau*(gA*gA - gV*gV)/(pow((Mmu2 - s),2.)*pow(Mtau2 - u,2));
  double Factor2= Mmu2*Mmu2*(MA2 + u) + 2*Mmu2*Mmu*(Mtau2*Mtau - Mtau*u) + Mmu2*(Mtau2*Mtau2 - 2*MA2*s - 2*Mtau2*u + u*(u-2*s)) + 2*Mmu*Mtau*s*(u - Mtau2) + s*(MA2*s + Mtau2*Mtau2 - 2*Mtau2*u + u*(s+u));
  double AmplZpr2WWVEGAS = Factor1*Factor2/(8*pi*s*s);
  // one should multiply the prefactor written below by factor
  // alphaEW^2 to get diff_CS_WW_phi in GeV^(-2)
  // (see e.g. calling VEGAS MC  function)
  double PrefactorWithoutAlphaEW=sqrt(x2-MA2/E02)/(1.0-XEv);
  double DsDxDthetaWithoutAlphaEW=sin(ThetaEv)*PrefactorWithoutAlphaEW*AmplZpr2WWVEGAS*ChiWWAnalytical;

  double ResTemporary;
  if (DsDxDthetaWithoutAlphaEW < 0.0 ) {
    ResTemporary = 0.0;
  } else {
    ResTemporary = DsDxDthetaWithoutAlphaEW;
  }
  return ResTemporary;
}


double DarkLFCScalars::Width()
{
  const double muMass = Mmu;
  const double tauMass = Mtau;
  const double massRatio2 = muMass*muMass/(MA*MA);
  double width          = 0.;
  double nuWidth        = 0.;
  double muWidth        = 0.;
  if (MA < 2.*tauMass) {
    nuWidth = epsil*epsil*alphaEW*(1./3.)*MA;
    if (MA > 2.*muMass) {
      double factor = (1.+2.*massRatio2)*sqrt(1.-4.*massRatio2);
      muWidth = nuWidth*factor;
    }
    width = nuWidth+muWidth; // in GeV
  } else {
    std::cout << "DarkLFCScalars: width for the mass above 2 Mtau is not implemented, exiting" << std::endl;
    exit(1);
  }
  return width;
}
