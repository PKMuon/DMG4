// This is a class to simulate dark X scalar boson production by muons in matter muN -> eNX
// Description is in the base class DarkMatter code
// To be used in a Geant4 application.
//
//
#include "DarkMatter.hh"
#include "DarkScalarLFC.hh"
#include "Utils.hh"

#include "G4MuonMinus.hh"
#include "G4TauMinus.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_dilog.h>

#include <stdexcept>
#include <iostream>
#include "G4ios.hh"

#include "KFactorsLFCmu.code"
#include "KFactorsLFCelectron.code"

// Auxiliary structures and functions:


// Additional structures holding E0 value to be passed into GSL callbacks.
struct BoundParms {
  DarkScalarLFC * this_;
  double E0;
};

struct BoundParms2 {
  DarkScalarLFC * this_;
  double E0;
  double X;
};

struct BoundParms3 {
  DarkScalarLFC * this_;
  double E0;
  double psi;
};


// photon flux for IWW, depends only on mass and energy
double DarkScalarLFC::chiIWW(double E0) {

  double tmin = MA*MA*MA*MA/(4.*E0*E0);
  // double tmax = MA*MA+MParent*MParent;
  double tmax = tMax;
  if(fabs(tMax - 10000.) < 0.001) tmax = E0*E0;

  double aa = 111.*pow(ZNucl,-1./3)/Mel;
  double d = 0.164*pow(ANucl,-2./3);
  double ta = pow(1./aa,2.);
  double td = d;
  double fluxAnalytical = ZNucl*ZNucl*(-((td*td*(((ta - td)*(ta + td + 2.0*tmax)*(tmax - tmin))/((ta + tmax)*(td + tmax)) + (ta + td + 2.0*tmin)*(log(ta + tmax) - log(td + tmax) - log(ta + tmin) + log(td + tmin))))/((ta-td)*(ta-td)*(ta-td))));

  return fluxAnalytical;

}

// A callback wrapping function for DarkScalarLFC::CrossSectionDSDX_IWW()
static double _DarkScalarLFCDsDxMuon(double x1, void * parms_) {
  //BoundParms * parms = (BoundParms*) parms_;  // or, equivalently, in C++ style
  BoundParms * parms = reinterpret_cast<BoundParms*>(parms_);
  // Forward invocation to target method
  return parms->this_->CrossSectionDSDX_IWW( x1, parms->E0 );
}


// A callback wrapping function for DarkScalarLFC::CrossSectionDSDX_WW()
static double _DarkScalarLFCDsDxMuon_WW(double x1, void * parms_) {
  //BoundParms * parms = (BoundParms*) parms_;  // or, equivalently, in C++ style
  BoundParms * parms = reinterpret_cast<BoundParms*>(parms_);
  // Forward invocation to target method
  return parms->this_->CrossSectionDSDX_WW( x1, parms->E0 );
}


// A callback wrapping function for DarkScalarLFC::CrossSectionDSDXDTheta()
static double _DarkScalarLFCDsDxDThetaMuon(double x[], size_t dim, void * parms_) {
  BoundParms * parms = reinterpret_cast<BoundParms*>(parms_);
  // Forward invocation to target method
  // x[1] = theta is jacobian of dcos to dtheta
  return parms->this_->CrossSectionDSDXDTheta( x[0], x[1], parms->E0 );
}

// A callback wrapping function for DarkScalarLFC::CrossSectionDSDXDTheta() for use in DSDX_WW integration over theta
static double _DarkScalarLFCDsDx_integrand(double theta1, void * parms_) {
  BoundParms2 * parms = reinterpret_cast<BoundParms2*>(parms_);
  // Forward invocation to target method
  // theta1 is jacobian from cos theta -> theta
  return parms->this_->CrossSectionDSDXDTheta( parms->X, theta1, parms->E0 );
}
static double _DarkScalarLFCDsDpsi_integrand(double y, void * parms_) {
  BoundParms3 * parms = reinterpret_cast<BoundParms3*>(parms_);
  double auxpsi = 0.5*parms->psi*parms->psi;
  return parms->this_->CrossSectionDSDXDPSI_IWW( 1.0-y, auxpsi, parms->E0 );
}

// A callback wrapping function for DarkScalarLFC::CrossSectionDSDXDpsi()
static double _DarkScalarLFCDsDxDPsiMuon(double x[], size_t dim, void * parms_) {
  BoundParms * parms = reinterpret_cast<BoundParms*>(parms_);
  // Forward invocation to target method
  double auxpsi = 0.5*x[1]*x[1];
  return parms->this_->CrossSectionDSDXDPSI_WW( x[0], auxpsi, parms->E0 );
}


// Class methods:  ---------------


DarkScalarLFC::DarkScalarLFC(double MAIn, double EThreshIn, double SigmaNormIn, double ANuclIn, double ZNuclIn, double DensityIn,
    double epsilIn, int IDecayIn, int ParentPDGIDIn)
: DarkMatter(MAIn, EThreshIn, SigmaNormIn, ANuclIn, ZNuclIn, DensityIn, epsilIn, IDecayIn)
{
  DMType = 45;
  // set electron or muon mode
  if(ParentPDGIDIn == 13 || ParentPDGIDIn == 11) ParentPDGID = ParentPDGIDIn;
  else if (ParentPDGIDIn == -9999) ParentPDGID = 11; // set default to electron mode
  else throw std::invalid_argument("ParentPDGID for DarkScalarLFC must be 11 or 13.");
  DaughterPDGID = 0;

  // set initial and final lepton masses
  if(ParentPDGID == 13) {MChild = Mel; MParent = Mmu;}
  if(ParentPDGID == 11) {MChild = Mmu; MParent = Mel;}

  IApprox =        2;   // Approximation: 1 - IWW; 2 - WW (default is 2)
  IMethodTotalCS = 4;   // Method for total CS: 1 - ds/dxdTheta; 2 - ds/dxdPsi; 3 - ds/dx (default is 1); 4 - use Kfactors
  tMax =       10000.;  // tmax initial; Value 10000. means that tmax = E0*E0 will be taken
  ThetaMax =     0.1;   // Max. angle of Z
  PsiMax =       3.0;   // Max. angle of recoil lepton

  std::cout << "Initialized Dark scalar LFC particle with ParentPDGID " << ParentPDGID << " for material density = " << DensityIn << std::endl;
  if(IApprox == 1) std::cout << "Using IWW approximation" << std::endl;
  if(IApprox == 2) std::cout << "Using WW approximation" << std::endl; 
  if(IMethodTotalCS == 1) std::cout << "ds/dxdTheta is used for total CS" << std::endl;
  if(IMethodTotalCS == 2) std::cout << "ds/dxdPsi is used for total CS" << std::endl;
  if(IMethodTotalCS == 3) std::cout << "ds/dx is used for total CS" << std::endl;
  std::cout << "ma = " << MA << " GeV" << std::endl;
  std::cout << "Energy cutoff = " << EThresh << " GeV" << std::endl;
  std::cout << std::endl;
}


DarkScalarLFC::~DarkScalarLFC()
{;}


double DarkScalarLFC::TotalCrossSectionCalc(double E0)
{
  if(IMethodTotalCS == 4) return TotalCrossSectionCalc_ETL(E0);

  if(IApprox == 1) {
    if(IMethodTotalCS == 3) return TotalCrossSectionCalc_IWW(E0); 
  }
  if(IApprox == 2) {
    if(IMethodTotalCS == 1) return TotalCrossSectionCalc_WW2(E0); // Integral of ds/dxdTheta
    if(IMethodTotalCS == 2) return TotalCrossSectionCalc_WW3(E0); // Integral of ds/dxdPsi
    if(IMethodTotalCS == 3) return TotalCrossSectionCalc_WW(E0);  // Integral of ds/dx
    std::cout << "DarkScalarLFC: wrong value of IMethodTotalCS, exiting" << std::endl;
    exit(1);
  } 
  else {
    std::cout << "DarkScalarLFC: wrong value of IApprox, exiting" << std::endl;
  }
  exit(1);
}

// This method computes the ETL for the total cross-section using the K-factors approximation
// no the IWW approach
//
double DarkScalarLFC::TotalCrossSectionCalc_ETL(double E0)
{
  double sigmaTot = TotalCrossSectionCalc_IWW(E0);
  double KFactor = 1;
  if(ParentPDGID == 13) KFactor = KfactorLFCmu2eApproximate(MA, E0);
  if(ParentPDGID == 11) KFactor = KfactorLFCe2muApproximate(MA, E0);
  else std::cout << "No KFactor for ParentPDGID = " << ParentPDGID << ". Using KFactor = 1." << std::endl;
  double result = sigmaTot / KFactor;
  return result;
}


// We integrate below the differential cross-section over x and t (see e.g. second line of Eq.(30) in 1705.01633,
// where x=E_Z'/E_0, t is a transfer momentum squared
//
double DarkScalarLFC::TotalCrossSectionCalc_IWW(double E0)
{
  // A. Ponten eq. (26 & 27) (Master thesis draft)
  if(E0 < 2.*MA) return 0.;

  gsl_integration_workspace* w1 = gsl_integration_workspace_alloc (1000);
  double result1, error1;

  double Xmin1=MA/E0;
  double Xmax1 = 1. - MA*MA*MA*MA/(8.*E0*E0*E0*ANucl) - MChild/E0;

  gsl_function F1;
  BoundParms parms = { this, E0 };
  F1.function = _DarkScalarLFCDsDxMuon;
  F1.params = &parms;

  gsl_integration_qags (&F1, Xmin1, Xmax1, 0, 1e-7, 1000, w1, &result1, &error1);

  double IntDsDx = result1;
  gsl_integration_workspace_free (w1);

  double Prefactor = chiIWW(E0)*epsil*epsil*alphaEW*alphaEW/(4*M_PI);

  double sigmaTot= GeVtoPb*Prefactor*IntDsDx;
  return sigmaTot;
}


// Below is the integration of WW ds/dx
//
double DarkScalarLFC::TotalCrossSectionCalc_WW(double E0)
{
  double sigmaTot;

  //if(E0 < 2.*MA) return 0.;
  //if(E0 < MA*99.) return 0.; // Some approximations probably don't work for smaller energy

  //to switch off default error handler, store old error handler in old_handler:
  gsl_error_handler_t * old_handler=gsl_set_error_handler_off();

  gsl_integration_workspace* w1 = gsl_integration_workspace_alloc (1000);
  double result1, error1;
  double Xmin1 = MA/E0;
  if(EThresh/E0 > Xmin1) Xmin1 = EThresh/E0;
  double Xmax1 = 1. - MA*MA*MA*MA/(8.*E0*E0*E0*ANucl) - MChild/E0;
  if(Xmax1 < Xmin1) return 0.;

  gsl_function F1;
  BoundParms parms = {this, E0};
  F1.function = _DarkScalarLFCDsDxMuon_WW; 
  F1.params = &parms;

  //gsl_integration_qags (&F1, Xmin1, Xmax1, 0, 1e-7, 1000, w1, &result1, &error1);
  double relerr=1.0e-9;   //initial error tolerance (relative error)
  int status=1;
  while(status) {
    status=gsl_integration_qags (&F1, Xmin1, Xmax1, 0, relerr, 1000, w1, &result1, &error1);
    relerr *= 1.2;
    if(status) G4cout << "Increased tolerance=" << relerr << G4endl;
  }
  //if integration routine returns error code, integration is repeated
  //using increased error tolerance, message is printed out
  gsl_set_error_handler(old_handler); //reset error handler (might be unneccessary.)

  double IntDsDx = result1;
  gsl_integration_workspace_free (w1);

  // include prefactor to get total cross section
  double prefactor = E0*epsil*epsil*alphaEW*alphaEW/(2.*M_PI);
  sigmaTot= GeVtoPb*prefactor*IntDsDx;

  return sigmaTot;
}


// Below is the 2 - dimensional integration of WW ds/dxdTheta
//
double DarkScalarLFC::TotalCrossSectionCalc_WW2(double E0)
{
  if(E0 < 2.*MA) return 0.;

  double Xmin1 = MA/E0;
  if(EThresh/E0 > Xmin1) Xmin1 = EThresh/E0;
  double Xmax1 = 1. - MA*MA*MA*MA/(8.*E0*E0*E0*ANucl) - MChild/E0;
  if(Xmax1 < Xmin1) return 0.;


  double xl[2] = { Xmin1, 0.};
  double xu[2] = { Xmax1, ThetaMax};

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_monte_function G;
  BoundParms parms = {this, E0};
  G.f = _DarkScalarLFCDsDxDThetaMuon; 
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

  // insert  prefactor to get correct total cross section
  double prefactor = E0*epsil*epsil*alphaEW*alphaEW/(2.*M_PI);
  sigmaTot = GeVtoPb*res*prefactor;

  // Other vegas integration methods:
  // gsl_monte_plain : plain MC
  // gsl_monte_vegas

#if 0
  gsl_monte_vegas_state* stat = gsl_monte_vegas_alloc(2);
  gsl_monte_vegas_integrate(&G, xl, xu, 2, 10000, r, s, &res, &err);
  do {
    gsl_monte_vegas_integrate (&G, xl, xu, 2, calls/5, r, s, &res, &err);
  } while (fabs (gsl_monte_vegas_chisq(stat) - 1.0) > 0.5);
  double sigmaTot = GeVtoPb*res*PrefactorEpsilonAlphaEWE0;
  gsl_monte_vegas_free(stat);
#endif

  gsl_rng_free (r);
  return sigmaTot;
}

// Below is the 2 - dimensional integration of WW ds/dxdPsi
//
double DarkScalarLFC::TotalCrossSectionCalc_WW3(double E0)
{
  if(E0 < 2.*MA) return 0.;

  double Xmin1 = MA/E0;
  if(EThresh/E0 > Xmin1) Xmin1 = EThresh/E0;
  double Xmax1 = 1. - MA*MA*MA*MA/(8.*E0*E0*E0*ANucl) - MChild/E0;
  if(Xmax1 < Xmin1) return 0.;


  double xl[2] = { Xmin1, 0.};
  double xu[2] = { Xmax1, PsiMax};

  const gsl_rng_type *T;
  gsl_rng *r;

  gsl_monte_function G;
  BoundParms parms = {this, E0};
  G.f = _DarkScalarLFCDsDxDPsiMuon; 
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
  double prefactor =  epsil*epsil*alphaEW*alphaEW/(2.*M_PI)*E0*E0; 
  sigmaTot = GeVtoPb*res*prefactor;

  // Other vegas integration methods:
  // gsl_monte_plain : plain MC
  // gsl_monte_vegas

#if 0
  gsl_monte_vegas_state* stat = gsl_monte_vegas_alloc(2);
  gsl_monte_vegas_integrate(&G, xl, xu, 2, 10000, r, s, &res, &err);
  do {
    gsl_monte_vegas_integrate (&G, xl, xu, 2, calls/5, r, s, &res, &err);
  } while (fabs (gsl_monte_vegas_chisq(stat) - 1.0) > 0.5);
  double sigmaTot = GeVtoPb*res*PrefactorEpsilonAlphaEWE0;
  gsl_monte_vegas_free(stat);
#endif

  gsl_rng_free (r);
  return sigmaTot;
}


double DarkScalarLFC::GetSigmaTot(double E0)
{
  return GetSigmaTot0(E0);
}


double DarkScalarLFC::CrossSectionDSDX(double XEv, double E0)
{
  if(IApprox == 1) return CrossSectionDSDX_IWW(XEv, E0);
  if(IApprox == 2) return CrossSectionDSDX_WW(XEv, E0);
  std::cout << "DarkScalarLFC: wrong value of IApprox, exiting" << std::endl;
  exit(1);
}


// Below is cross section obtained by integration of IWW ds/dxdTheta
// chiIWW and constant prefactors not included
double DarkScalarLFC::CrossSectionDSDX_IWW(double XEv, double E0)
{
  // NOTE: this function returns dsdx_IWW/(alpha^2*h^2/(4pi)*chi_IWW) 
  if(XEv*E0 <= MA) return 0.;
  // A. Ponten eq. (27) (Master thesis draft)
  double umin = -XEv*E0*E0*ThetaMax*ThetaMax-MA*MA*(1.-XEv)/XEv+MParent*MParent*(1.-XEv)-MChild*MChild;
  double umax = -MA*MA*(1.-XEv)/XEv+MParent*MParent*(1.-XEv)-MChild*MChild;
  double deltaM2 = MA*MA-MParent*MParent-MChild*MChild;

  double AmpAtUmax = XEv*XEv/(-2.*(1.-XEv)*umax)+deltaM2*XEv/(-2.*umax*umax)+deltaM2*(MA*MA*(1.-XEv)+XEv*(MChild*MChild+MParent*MParent*(XEv-1.)))/(-3.*umax*umax*umax);
  double AmpAtUmin = XEv*XEv/(-2.*(1.-XEv)*umin)+deltaM2*XEv/(-2.*umin*umin)+deltaM2*(MA*MA*(1.-XEv)+XEv*(MChild*MChild+MParent*MParent*(XEv-1.)))/(-3.*umin*umin*umin);
  double Amp = AmpAtUmax-AmpAtUmin;

  double beta = sqrt(1.-MA*MA/(XEv*XEv*E0*E0));

  double Prefactor = (1.-XEv)*beta;

  double res = Prefactor*Amp;
  // res = Prefactor*Amp*chiIWW(E0)*epsil*epsil*alphaEW*alphaEW/(4*M_PI);
  return res;
}


// Below is cross section obtained by integration of WW ds/dxdTheta
//
double DarkScalarLFC::CrossSectionDSDX_WW(double XEv, double E0)
{
  // ds/dx
  if( XEv*E0 <= MA ) return 0.0; 

  // integration over theta
  gsl_integration_workspace* w1 = gsl_integration_workspace_alloc (1000);
  double result1, error1;
  gsl_function F1;
  BoundParms2 parms = { this, E0, XEv };
  F1.function = _DarkScalarLFCDsDx_integrand;
  F1.params = &parms;
  gsl_integration_qags(&F1, 0.0, ThetaMax, 0, 1e-7, 1000,  w1, &result1, &error1);
  double res = result1;
  gsl_integration_workspace_free(w1);
  return res;
}

// wrapper for DSDXDU, used in electron mode sampling
double DarkScalarLFC::CrossSectionDSDXDU(double XEv, double UThetaEv, double E0)
{
  if(IApprox == 1) return CrossSectionDSDXDU_IWW(XEv, UThetaEv, E0);
  if(IApprox == 2) return CrossSectionDSDXDU_WW(XEv, UThetaEv, E0);
  std::cout << "DarkScalarLFC: wrong value of IApprox, exiting" << std::endl;
  exit(1);
}

// Below is the double-differential cross-section as for Theta expressed
// as u.  UThetaEv = 0.5*theta^2
double DarkScalarLFC::CrossSectionDSDXDU_IWW(double XEv, double UThetaEv, double E0)
{
  // du = d(p^2+k^2-2pk - me^2) = d(-2(xE0^2 - |p||k|costheta)) = 2|p||k|dcos, so divide dsdxdcos with 2|p||k| to get dsdxdu
  // in IWW approximation, this is 2*E0*XEv*E0*beta. If we omit constant prefactors, then we divide dsdxdcos(theta) by XEv*beta to get dsdxdu

  // double prefactor = 1.0/(2*E0)*chiIWW(E0)*E0*epsil*epsil*alphaEW*alphaEW/(2.*M_PI);
  double deltaM2 = MA*MA-MParent*MParent-MChild*MChild;
  double u = -XEv*E0*E0*2*UThetaEv-MA*MA*(1.-XEv)/XEv+MParent*MParent*(1.-XEv)-MChild*MChild;
  double amp = XEv*XEv/(2.-2.*XEv) + deltaM2*(u*XEv+MA*MA*(1. - XEv) + XEv*(MChild*MChild+MParent*MParent*(XEv - 1.)) )/(u*u); // eq. 18
  double res = (1. - XEv)*amp/(u*u);
  return res;
}

// Below is the double-differential cross-section as for Theta expressed
// as u.  UThetaEv = 0.5*theta^2
double DarkScalarLFC::CrossSectionDSDXDU_WW(double XEv, double UThetaEv, double E0)
{
  // du = d(p^2+k^2-2pk - me^2) = d(-2(xE0^2 - |p||k|costheta)) = 2|p||k|dcos, so divide dsdxdcos with 2|p||k| to get dsdxdu
  // in WW approximation, this is 2*E0*XEv*E0*beta. If we omit constant prefactors, then we divide dsdxdcos(theta) by XEv*beta to get dsdxdu

  // double prefactor = 1.0/(2*E0)*chiIWW(E0)*E0*epsil*epsil*alphaEW*alphaEW/(2.*M_PI);
  double deltaM2 = MA*MA-MParent*MParent-MChild*MChild;
  double u = -XEv*E0*E0*2*UThetaEv-MA*MA*(1.-XEv)/XEv+MParent*MParent*(1.-XEv)-MChild*MChild;
  double amp = XEv*XEv/(2.-2.*XEv) + deltaM2*(u*XEv+MA*MA*(1. - XEv) + XEv*(MChild*MChild+MParent*MParent*(XEv - 1.)) )/(u*u); // eq. 18

  // photon flux
  double tmin = u*u/(4.*E0*E0*(1.0-XEv)*(1.0-XEv));
  // double tmax = MA*MA+MParent*MParent;
  double tmax = tMax;
  if(fabs(tMax - 10000.) < 0.001) tmax = E0*E0;
  double aa = 111.*pow(ZNucl,-1./3)/Mel;
  double d = 0.164*pow(ANucl,-2./3);
  double ta = pow(1./aa,2.);
  double td = d;
  double fluxAnalytical = ZNucl*ZNucl*(-((td*td*(((ta - td)*(ta + td + 2.0*tmax)*(tmax - tmin))/((ta + tmax)*(td + tmax)) + (ta + td + 2.0*tmin)*(log(ta + tmax) - log(td + tmax) - log(ta + tmin) + log(td + tmin))))/((ta-td)*(ta-td)*(ta-td))));

  double res = fluxAnalytical*(1. - XEv)*amp/(u*u);
  return res;
}

double DarkScalarLFC::CrossSectionDSDXDPSI(double XEv, double auxpsi, double E0)
{
  if(IApprox == 1) return CrossSectionDSDXDPSI_IWW(XEv, auxpsi, E0);
  if(IApprox == 2) return CrossSectionDSDXDPSI_WW(XEv, auxpsi, E0);
  std::cout << "DarkScalarLFC: wrong value of IApprox, exiting" << std::endl;
  exit(1);
}


// IWW double-differential cross-section for the total cross section, based
// on ds/dydpsi. For sampling reasons, y -> 1.-XEv
double DarkScalarLFC::CrossSectionDSDXDPSI_IWW(double XEv, double auxpsi, double E0)
{
  // In IWW approach we suppose that flux factor \chi doesn't depend on y and psi
  // in that method y=Ee/E0 muon energy fraction and auxpsi=psi^2/2 is an auxiliar variable
  // to sample muon deflection angle psi at first step, the accepted angle will be then psiAcc = sqrt(2.0*auxpsi)
  //       
  //

  //  if(E0*XEv < EThresh) return 0.0;
  double y = 1. - XEv;
  double psi = sqrt(2.0*auxpsi);

  // NOTE: test example of A. Ponten eq. (32) (Master thesis draft)
  double deltaM2 = MA*MA-MParent*MParent-MChild*MChild;
  double t2 = -y*E0*E0*2.0*auxpsi-MChild*MChild*XEv/y+MParent*MParent*XEv;
  double t = MA*MA-t2;
  double Amp2 = (1./(2.*t*t))-deltaM2/(t*t*t)+deltaM2*(MChild*MChild+y*(MParent*MParent*y+deltaM2))/(y*t*t*t*t);

  // NOTE: factor omitted are flux and EW and coupling 
  double beta = sqrt(1.-MChild*MChild/(y*y*E0*E0));
  // double Prefactor = (1.-y)*(1.-y)*(1.-y)*E0*E0*beta;
  double Prefactor = XEv*XEv*XEv*beta;
  double res = Prefactor*Amp2;
  res = sin(psi)*res; // don't forget your change of variables :)
  return res;
}


// WW double-differential cross-section for the total cross section, based
// on ds/dydPsi. For sampling reasons, y -> 1.-XEv
//
double DarkScalarLFC::CrossSectionDSDXDPSI_WW(double XEv, double auxpsi, double E0)
{

  //  if(E0*XEv < EThresh) return 0.0;

  double y = 1. - XEv;
  double psi = sqrt(2.0*auxpsi);

  double deltaM2 = MA*MA-MParent*MParent-MChild*MChild;
  double t2 = -y*E0*E0*2*auxpsi-MChild*MChild*(1.-y)/y+MParent*MParent*(1.-y);
  double t = MA*MA-t2;
  double Amp2 = (1./(2.*t*t))-deltaM2/(t*t*t)+deltaM2*(MChild*MChild+y*(MParent*MParent*y+deltaM2))/(y*t*t*t*t);

  // photon flux
  double tmin = (t/(2.*E0*(1.-y)))*(t/(2.*E0*(1.-y)));
  // double tmax = MA*MA+MParent*MParent;
  double tmax = tMax;
  if(fabs(tMax - 10000.) < 0.001) tmax = E0*E0;
  double aa = 111.*pow(ZNucl,-1./3)/Mel;
  double d = 0.164*pow(ANucl,-2./3);
  double ta = pow(1./aa,2.);
  double td = d;
  double fluxAnalytical = ZNucl*ZNucl*(-((td*td*(((ta - td)*(ta + td + 2.0*tmax)*(tmax - tmin))/((ta + tmax)*(td + tmax)) + (ta + td + 2.0*tmin)*(log(ta + tmax) - log(td + tmax) - log(ta + tmin) + log(td + tmin))))/((ta-td)*(ta-td)*(ta-td))));

  double beta = sqrt(1.-MChild*MChild/(y*y*E0*E0));
  // double Prefactor = epsil*epsil*alphaEW*alphaEW/(2.*M_PI)*(1.-y)*(1.-y)*(1.-y)*E0*E0*beta;
  double Prefactor = (1.-y)*(1.-y)*(1.-y)*beta;
  double res = Prefactor*fluxAnalytical*Amp2;
  res = sin(psi)*res; // don't forget your change of variables :)
  return res;
}


// WW double-differential cross-section for the total cross section, based 
// on ds/dxdTheta
//
double DarkScalarLFC::CrossSectionDSDXDTheta(double XEv, double ThetaEv, double E0)
{
  if(E0*XEv < EThresh) return 0.0;

  // prefactor and amplitude 
  // double prefactor = E0*epsil*epsil*alphaEW*alphaEW/(2.*M_PI);
  double deltaM2 = MA*MA-MParent*MParent-MChild*MChild;
  double momentumOfDP = sqrt(XEv*XEv*E0*E0 - MA*MA);
  double u = -XEv*E0*E0*ThetaEv*ThetaEv-MA*MA*(1.-XEv)/XEv+MParent*MParent*(1.-XEv)-MChild*MChild;
  double amp = XEv*XEv/(2.-2.*XEv) + deltaM2*(u*XEv+MA*MA*(1. - XEv) + XEv*(MChild*MChild+MParent*MParent*(XEv - 1.)) )/(u*u); // eq. 18

  // photon flux
  double tmin = u*u/(4.*E0*E0*(1.0-XEv)*(1.0-XEv));
  // double tmax = MA*MA+MParent*MParent;
  double tmax = tMax;
  if(fabs(tMax - 10000.) < 0.001) tmax = E0*E0;
  double aa = 111.*pow(ZNucl,-1./3)/Mel;
  double d = 0.164*pow(ANucl,-2./3);
  double ta = pow(1./aa,2.);
  double td = d;
  double fluxAnalytical = ZNucl*ZNucl*(-((td*td*(((ta - td)*(ta + td + 2.0*tmax)*(tmax - tmin))/((ta + tmax)*(td + tmax)) + (ta + td + 2.0*tmin)*(log(ta + tmax) - log(td + tmax) - log(ta + tmin) + log(td + tmin))))/((ta-td)*(ta-td)*(ta-td))));

  double res = sin(ThetaEv)*fluxAnalytical*momentumOfDP*(1. - XEv)*amp/(u*u);
  return res;
}

double DarkScalarLFC::Width()
{
  // should I just keep it at 0?
  double width = 0.;
#if 0
  if (MA > MChild+MParent) {
    double p = sqrt(MA*MA-pow(MChild+MParent,2.)*(MA*MA-pow(MChild-MParent,2.)))/MA*0.5;
    width = epsil*epsil*p/(8*M_PI)*(1.0-(MChild*MChild+MParent*MParent)/(MA*MA));
  }
  else {
    G4cout << "DarkScalarLFC: width for the mass below Mel+Mmu is not valid, exiting" << G4endl;
    exit(1);
  }
#endif
  return width;
}

// for plotting purposes
double DarkScalarLFC::CrossSectionDSDX_WW_withprefactor(double XEv, double E0) {
  // multiplies CrossSectionDSDX_WW with its prefactor
  double prefactor = E0*epsil*epsil*alphaEW*alphaEW/(2.*M_PI);
  return prefactor*CrossSectionDSDX_WW(XEv, E0);
}

double DarkScalarLFC::CrossSectionDSDX_IWW_withprefactor(double XEv, double E0) {
  // multiplies CrossSectionDSDX_IWW with its prefactor
  double prefactor = chiIWW(E0)*epsil*epsil*alphaEW*alphaEW/(4.0*M_PI);
  return prefactor*CrossSectionDSDX_IWW(XEv, E0);
}

double DarkScalarLFC::CrossSectionDSDPSI_IWW_noprefactor(double Psi, double E0) {

  double ymin = MChild/E0 + MA*MA*MA*MA/(8.*E0*E0*E0*ANucl);
  double ymax = 1.0 - MA/E0;
  // integration over y
  gsl_integration_workspace* w1 = gsl_integration_workspace_alloc (1000);
  double result1, error1;
  gsl_function F1;
  BoundParms2 parms = { this, E0, Psi };
  F1.function = _DarkScalarLFCDsDpsi_integrand;
  F1.params = &parms;
  gsl_integration_qags(&F1, ymin, ymax, 0, 1e-7, 1000,  w1, &result1, &error1);
  double res = result1;
  gsl_integration_workspace_free(w1);
  return res;

}
