#include "DarkMatterPhysics.hh"

#include "DarkMatter.hh"
#include "DarkPhotons.hh"
#include "DarkZ.hh"
#include "ALP.hh"
#include "DarkPhotonsAnnihilation.hh"
#include "DarkScalarsAnnihilation.hh"
#include "DarkScalars.hh"
#include "DarkPseudoScalars.hh"
#include "DarkAxials.hh"

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
  
  G4double BiasSigmaFactor0 = DMpar->GetRegisteredParam("Bias", 9.e8);
  G4double EThresh = DMpar->GetRegisteredParam("EThresh", 35.); // for sensitivity calculations invisible mode
  //G4double EThresh = 18.; // for sensitivity calculations visible mode
  //G4double EThresh = 1.; // for shape studies
  //G4double EThresh = 2000.; // to turn off A emissions  

  //select particle type and details
  G4int DMType    = DMpar->GetRegisteredParam("DMType"     ,0.     );
  double DMMass   = DMpar->GetRegisteredParam("MassA", 0.0167);
  double Epsilon  = DMpar->GetRegisteredParam("Epsil", 0.0001);
  // Initialize by default for Pb with eps=0.0001
  G4double ANucl     = DMpar->GetRegisteredParam("ANucl"      ,207.   );
  G4double ZNucl     = DMpar->GetRegisteredParam("ZNucl"      ,82.    );
  G4double Density   = DMpar->GetRegisteredParam("Density"    ,11.35  );
  G4int DecayType = DMpar->GetRegisteredParam("DecayType"  ,0     ); //0 invisible, 2 visible

  //additional pareameters
  G4double AlphaD = DMpar->GetRegisteredParam("AlphaD", 0.1);
  
  switch(DMType)
    {
    case 0:
      G4cout << "Initialize DarkPhoton\n";
      myDarkMatter = new DarkPhotons(DMMass, EThresh, 1., ANucl, ZNucl, Density, Epsilon, DecayType); 
      break;
    case 1:
      G4cout << "Initialize DarkScalar\n";
      myDarkMatter = new DarkScalars(DMMass, EThresh, 1., ANucl, ZNucl, Density, Epsilon, DecayType);
      break;
    case 2:
      G4cout << "Initialize DarkPseudoScalars\n";
      myDarkMatter = new DarkPseudoScalars(DMMass, EThresh, 1., ANucl, ZNucl, Density,  Epsilon, DecayType);
    case 3:
      G4cout << "Initialize DarkAxials\n";
      myDarkMatter = new DarkAxials(DMMass, EThresh, 1., ANucl, ZNucl, Density,  Epsilon, DecayType);
    case 4:
      G4cout << "Initialize ALP\n";
      myDarkMatter = new ALP(DMMass, EThresh, 1., ANucl, ZNucl, Density,  Epsilon, DecayType);
    case 5:
      G4cout << "Initialize DarkZ\n";
      myDarkMatter = new DarkZ(DMMass, EThresh, 1., ANucl, ZNucl, Density,  Epsilon, DecayType);
    default:
      G4cout << "Initialize DarkPhoton\n";
      myDarkMatter = new DarkPhotons(DMMass, EThresh, 1., ANucl, ZNucl, Density, Epsilon, DecayType); 
      break;      
    }

  // End of code for the visible mode

  BiasSigmaFactor = BiasSigmaFactor0 * 0.0001 * 0.0001 / (myDarkMatter->Getepsil()*myDarkMatter->Getepsil());

  if(!myDarkMatter) return false;
  return true;
}
