// This is a class to simulate conversion of gammas to ALP in matter gammaN -> AN
// Description to follow
// To be used in a Geant4 application.
//
//
#include "ALP.hh"
#include "Utils.hh"

#include <math.h>
#include <iostream>


ALP::ALP(double MAIn, double EThreshIn, double SigmaNormIn, double ANuclIn, double ZNuclIn, double DensityIn,
         double epsilIn, int IDecayIn)
: DarkMatter(MAIn, EThreshIn, SigmaNormIn, ANuclIn, ZNuclIn, DensityIn, epsilIn, IDecayIn)
{
  DMType = 21;
  ParentPDGID = 22;
  DaughterPDGID = 22;
  PrepareVariables();
  std::cout << "Initialized ALP (gamma conversion to ALP) for material density = " << DensityIn << std::endl;
  std::cout << std::endl;
}


ALP::~ALP()
{;}


double ALP::TotalCrossSectionCalc(double E0)
{
  if(E0 < 2.*MA) return 0.; // TODO: the exact threshold is not so big (note that it is present also in the DarkMatter methods)

  double tmin = MA*MA*MA*MA/(4.*E0*E0);
  //double tmax = MA*MA;

  double GaggBench=epsilBench; // Dimensional coupling of ALP in 1/GeV
  double atomFFcoeff=111.0*pow(ZNucl,-1.0/3.0)/Mel;
  double tAtom=1.0/(atomFFcoeff*atomFFcoeff); // atomic coefficient of form-factor
  double tNucl=0.164*pow(ANucl,-2.0/3.0); // nuclear form-factor coefficient in GeV**2
  double LogFactor=log((tNucl+tmin)/(tAtom+tmin))-2.0; // see Note_ALP.pdf and check it
  double sigmaALPtotal=1.0/8.0*GaggBench*GaggBench*alphaEW*ZNucl*ZNucl*LogFactor*GeVtoPb;

  //G4cout << "Total CS calc, E, M, cs = " << E0 << " " << MA << " " << sigmaALPtotal << G4endl;

  return sigmaALPtotal; 
}


double ALP::GetSigmaTot(double E0)
{
  return TotalCrossSectionCalc(E0);
}


double ALP::CrossSectionDSDX(double XEv, double E0)
{
  if(XEv > 0.999) return 1.;
  return 0.;
}


double ALP::CrossSectionDSDXDU(double XEv, double UThetaEv, double E0)
{
  if(XEv > 0.999) return 1.;
  return 0.;
}


double ALP::Width()
{
  return 1./(64.*3.1415926)*MA*MA*MA*epsil*epsil;
}


/**
 * Differential cross-section of Ngamma->Na w.r.t ALP emission angle
 * WARNING: Missing prefactor since this method is only relevant for sampling
 *
 * For details see e.g. formulas from  https://arxiv.org/pdf/2004.04469.pdf
 *
 * MA is a mass of ALP (GeV)
 * @params E0 energy of incoming photon (GeV)
 * @params ThetaEv is the ALP emission angle
 * @return differential cross-section at given energy and angle values (a.u.)
 *
 */
double ALP::CrossSectionDSDTheta(double ThetaEv, double E0)
{
  double tmin = MA*MA*MA*MA/(4.*E0*E0); //
  double atomFFcoeff=111.0*pow(ZNucl,-1.0/3.0)/Mel;//a-parameter for atomic coefficient of form-factor
  double tAtom=1.0/(atomFFcoeff*atomFFcoeff); // atomic coefficient of form-factor
  double tNucl=0.164*pow(ANucl,-2.0/3.0); // nuclear form-factor coefficient in GeV**2
  double ThetaEv2=ThetaEv*ThetaEv;
  double tAsFuncOfThetaALP= E0*E0*ThetaEv2 +tmin; // see Eq.(11) from 2004.04469, momentum trasnfer squared as a function of ThetaEv
  double FFpref1=tAsFuncOfThetaALP/(1.0+tAsFuncOfThetaALP/tAtom); // atomic Formfactor term without ~ Z a^2 (see Eq.(5) from 2004.04469)
  double FFpref2= 1.0/(1.0+tAsFuncOfThetaALP/tNucl); // nuclear Formfactor term, it seems to be that it does not play important role as soon as ThetaEv < 1.0e-2
  double FFtotal=FFpref1*FFpref2; // total formfactor
  double FFtotalSquared= FFtotal*FFtotal;// Note that it depends on ThetaEv via transfer momentum squared tAsFuncOfThetaALP!!!!
  double deltaALP=MA*MA/(2.0*E0*E0); // see Eq.(16) from 2004.04469
  double deltaALP2=deltaALP*deltaALP;
  double ThetaDependence= ThetaEv2*ThetaEv/((ThetaEv2+deltaALP2)*(ThetaEv2+deltaALP2));
  double DsigmaDthetaWithoutPrefactor=FFtotalSquared*ThetaDependence; // general form of Eq. (16) in 2004.04469 without prefactor (only ThetaEv dependence!!!)

  return DsigmaDthetaWithoutPrefactor;
}


/**
 * Differential cross-section of Ngamma->Na w.r.t ALP emission angle
 * WARNING: Missing prefactor since this method is only relevant for sampling
 *
 * For details see e.g. formulas from  https://arxiv.org/pdf/2004.04469.pdf
 * For thetaMAX see Eq.(17) from 2004.04469, thetaMAX is the ALP emission angle for which the differential cross-section is maximal
 *
 * @params E0 energy of incoming photon (GeV)
 * @return differential cross-section at given energy and angle values (a.u.)
 *
 */
double ALP::CrossSectionDSDThetaMAX(double E0)
{

  double tmin = MA*MA*MA*MA/(4.*E0*E0); // minimum momentum transfer
  double atomFFcoeff=111.0*pow(ZNucl,-1.0/3.0)/Mel; // a-parameter for atomic coefficient of form-factor
  double tAtom=1.0/(atomFFcoeff*atomFFcoeff); // atomic coefficient of form-factor
  double tNucl=0.164*pow(ANucl,-2.0/3.0); // nuclear form-factor coefficient in GeV**2
  double thetaMAX= 1.0/(E0*atomFFcoeff)*sqrt(3.0*(1.0+tmin/tAtom));
  double thetaMAX2=thetaMAX*thetaMAX;
  double tAsFuncOfThetaMAXALP= E0*E0*thetaMAX2 +tmin; // see Eq.(11) from 2004.04469, momentum trasnfer squared as a function of theta
  double FFpref1=tAsFuncOfThetaMAXALP/(1.0+tAsFuncOfThetaMAXALP/tAtom); // atomic Formfactor term without ~ Z a^2 (see Eq.(5) from 2004.04469)
  double FFpref2= 1.0/(1.0+tAsFuncOfThetaMAXALP/tNucl); // nuclear Formfactor term, it seems to be that it does not play important role as soon as t>
  double FFtotal=FFpref1*FFpref2; // total formfactor
  double FFtotalSquared= FFtotal*FFtotal;// Note that it depends on thetaMAX via transfer momentum squared tAsFuncOfThetaMAXALP!!!!
  double deltaALP=MA*MA/(2.0*E0*E0); // see Eq.(16) from 2004.04469
  double deltaALP2=deltaALP*deltaALP;
  double ThetaMAXDependence= thetaMAX2*thetaMAX/((thetaMAX2+deltaALP2)*(thetaMAX2+deltaALP2));
  double DsigmaDthetaMAXWithoutPrefactor=FFtotalSquared*ThetaMAXDependence; // general form of Eq.  (16) in 2004.04469 without prefactor (only theta dependence!!!)

  return DsigmaDthetaMAXWithoutPrefactor;
}
