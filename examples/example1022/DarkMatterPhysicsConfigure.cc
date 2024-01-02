#include "DarkMatterPhysics.hh"

#include "DarkMatter.hh"
#include "DarkPhotons.hh"
#include "DarkScalars.hh"
#include "DarkPseudoScalars.hh"
#include "DarkAxials.hh"
#include "ALP.hh"
#include "DarkZ.hh"
#include "DarkPhotonsAnnihilation.hh"
#include "DarkScalarsAnnihilation.hh"
#include "DarkPseudoScalarsAnnihilation.hh"
#include "DarkAxialsAnnihilation.hh"


#include "DarkMatterParametersFactory.hh"

#include "G4SystemOfUnits.hh"


// BiasSigmaFactor Invisible mode Vector EThresh=35
// 900.  9.e12
// 16.7  8.e8
//  5.   1.75e8
//  2.   3.5e7
// 0.5   1.2e7
// 0.1   3.6e6
// 0.01  1.55e6
// 0.002 1.38e6

// BiasSigmaFactor Invisible mode Scalar EThresh=35
// 16.7  2.3e9

// BiasSigmaFactor Visible mode Vector EThresh=18
// 16.7  3.4e8

bool DarkMatterPhysics::DarkMatterPhysicsConfigure() 
{
  //call an instance of the class
  DarkMatterParametersFactory* DMpar = DarkMatterParametersFactory::GetInstance();  
  
  DMpar->RegisterNewParam("BiasSigmaFactor0", 6.e9);
  DMpar->RegisterNewParam("EThresh", 18.*GeV); // for sensitivity calculations invisible mode
  //G4double EThresh = 18.; // for sensitivity calculations visible mode
  //G4double EThresh = 1.; // for shape studies
  //G4double EThresh = 2000.; // to turn off A emissions  

  //select particle type and details
  DMpar->RegisterNewParam("DMProcessType", 21.);
  DMpar->RegisterNewParam("DMMass", 0.03*GeV);
  DMpar->RegisterNewParam("Epsilon", 0.0009);

  // Initialize for Pb
  DMpar->RegisterNewParam("ANucl"      ,207.   );
  DMpar->RegisterNewParam("ZNucl"      ,82.    );
  DMpar->RegisterNewParam("Density"    ,11.35*(g/cm3) );

  // Initialize for W
//  DMpar->RegisterNewParam("ANucl"   ,184.   );
//  DMpar->RegisterNewParam("ZNucl"   ,74.    );
//  DMpar->RegisterNewParam("Density" ,19.25*(g/cm3) );

  DMpar->RegisterNewParam("DecayType", 2.); // 0 invisible, 2 visible

  // additional parameters for annihilation
//  DMpar->RegisterNewParam("RDM", 1./3.);
//  DMpar->RegisterNewParam("AlphaD", 0.5);
  
  return true;
}
