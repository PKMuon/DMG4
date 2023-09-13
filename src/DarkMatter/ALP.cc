// This is a class to simulate conversion of gammas to ALP in matter gammaN -> AN
// Description to follow
// To be used in a Geant4 application.
//
//
#include "DarkMatter.hh"
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

//This is just a very rough prototype of a NEW METHOD for ALP angle sampling \theta
////it returns just (d\sigma/d\theta) up to arbitrary prefactor (that
////does not depend on \theta). The method CrossSectionDSDTheta  should be added to DM general
//class.
// For details see e.g. formulas from  https://arxiv.org/pdf/2004.04469.pdf
double ALP::CrossSectionDSDTheta( double E0)
{
  double tmin = MA*MA*MA*MA/(4.*E0*E0); //
  // E0 is a energy of incoming photon (GeV)
  // MA is a mass of ALP (GeV)
  // theta is an ALP emission angle
  double atomFFcoeff=111.0*pow(ZNucl,-1.0/3.0)/Mel;//a-parameter for atomic coefficient of form-factor
  double tAtom=1.0/(atomFFcoeff*atomFFcoeff); // atomic coefficient of form-factor
  double tNucl=0.164*pow(ANucl,-2.0/3.0); // nuclear form-factor coefficient in GeV**2
  double tAsFuncOfThetaALP= E0*E0*theta*theta +tmin; // see Eq.(11) from 2004.04469, momentum trasnfer squared as a function of theta
  double FFpref1=tAsFuncOfThetaALP/(1.0+tAsFuncOfThetaALP/tAtom); // atomic Formfactor term without ~ Z a^2 (see Eq.(5) from 2004.04469)
  double FFpref2= 1.0/(1.0+tAsFuncOfThetaALP/tNucl); // nuclear Formfactor term, it seems to be that it does not play important role as soon as theta < 1.0e-2
  double FFtotal=FFpref1*FFpref2; // total formfactor
  double FFtotalSquared= FFtotal*FFtotal;// Note that it depends on theta via transfer momentum squared tAsFuncOfThetaALP!!!!
  double deltaALP=MA*MA/(2.0*E0*E0); // see Eq.(16) from 2004.04469
  double deltaALP2=deltaALP*deltaALP;
  double theta2=theta*theta;
  double MA2=MA*MA;
  double ThetaDependence= theta2*theta/((theta2+deltaALP2)*(theta2+deltaALP2));
  double DsigmaDthetaWithoutPrefactor=FFtotalSquared*ThetaDependence; // general form of Eq. (16) in 2004.04469 without prefactor (only theta dependence!!!)

  // Misha this is just a link between X and theta: Xaccepted = Ea/Egamma from Eq. (10) in 2004.04469
  double Xaccepted =1.-E0*theta2/(2.0*ANucl)-MA2*MA2/(8.0*ANucl*E0*E0*E0); // we assume here that Mnucleus = ANucl

  return DsigmaDthetaWithoutPrefactor;
}




// we perform below the maximum value of the DsigmaDthetaWithoutPrefactor function
// that implies substituted thetaMAX (see Eq.(17) from 2004.04469)
double ALP::CrossSectionDSDThetaMAX( double E0)
{

  double tmin = MA*MA*MA*MA/(4.*E0*E0); //
  // E0 is a energy of incoming photon (GeV)
  // MA is a mass of ALP (GeV)
  // thetaMAX is a ALP emission angle that provides maximum value of the differential cross section
  double atomFFcoeff=111.0*pow(ZNucl,-1.0/3.0)/Mel;//a-parameter for atomic coefficient of form-factor
  double tAtom=1.0/(atomFFcoeff*atomFFcoeff); // atomic coefficient of form-factor
  double tNucl=0.164*pow(ANucl,-2.0/3.0); // nuclear form-factor coefficient in GeV**2
  double thetaMAX= 1.0/(E0*atomFFcoeff)*sqrt(3.0*(1.0+tmin/tAtom)); // (see Eq.(17) from 2004.04469)
  double tAsFuncOfThetaMAXALP= E0*E0*thetaMAX*thetaMAX +tmin; // see Eq.(11) from 2004.04469, momentum trasnfer squared as a function of theta
  double FFpref1=tAsFuncOfThetaMAXALP/(1.0+tAsFuncOfThetaMAXALP/tAtom); // atomic Formfactor term without ~ Z a^2 (see Eq.(5) from 2004.04469)
  double FFpref2= 1.0/(1.0+tAsFuncOfThetaMAXALP/tNucl); // nuclear Formfactor term, it seems to be that it does not play important role as soon as t>
  double FFtotal=FFpref1*FFpref2; // total formfactor
  double FFtotalSquared= FFtotal*FFtotal;// Note that it depends on thetaMAX via transfer momentum squared tAsFuncOfThetaMAXALP!!!!
  double deltaALP=MA*MA/(2.0*E0*E0); // see Eq.(16) from 2004.04469
  double deltaALP2=deltaALP*deltaALP;
  double thetaMAX2=thetaMAX*thetaMAX;
  double MA2=MA*MA;
  double ThetaMAXDependence= thetaMAX2*thetaMAX/((thetaMAX2+deltaALP2)*(thetaMAX2+deltaALP2));
  double DsigmaDthetaMAXWithoutPrefactor=FFtotalSquared*ThetaMAXDependence; // general form of Eq.  (16) in 2004.04469 without prefactor (only theta dep>

  return DsigmaDthetaMAXWithoutPrefactor;
}


