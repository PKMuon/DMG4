#include "DMParticleAPrime.hh"
#include "DarkMatterParametersFactory.hh"

#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DalitzDecayChannel.hh"
#include "G4DecayTable.hh"
#include "G4MuonMinus.hh"
#include "G4PionZero.hh"

DMParticleAPrime * DMParticleAPrime::theInstance = nullptr;

DMParticleAPrime* DMParticleAPrime::Definition()
{
  if( theInstance ) {
    return theInstance;
  }
  //get parameters from factory (NOTE: mass is parsed in GeV)
  DarkMatterParametersFactory* DMpar = DarkMatterParametersFactory::GetInstance();
  G4double MassIn    = DMpar->GetRegisteredParam("DMMass")*GeV;
  G4double epsilIn   = DMpar->GetRegisteredParam("Epsilon");
  G4double DecayType = DMpar->GetRegisteredParam("DecayType");  
  G4double BranchingType = DMpar->GetRegisteredParam("BranchingType", 0);

  const G4String name = "DMParticleAPrime";
  // search in particle table]
  G4ParticleTable * pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition * anInstance = pTable->FindParticle(name);
  const G4double muMass = G4MuonMinus::MuonMinusDefinition()->GetPDGMass();
  const G4double pi0Mass = G4PionZero::PionZeroDefinition()->GetPDGMass();
  G4double RatioEA2 = electron_mass_c2*electron_mass_c2/(MassIn*MassIn);
  G4double massRatio2 = muMass*muMass/(MassIn*MassIn);
  G4bool isStable = true;
  G4double WidthIn = 0.;
  G4double nuWidth = 0.;
  G4double eWidth = 0.;
  G4double muWidth = 0.;
  G4double hWidth = 0.;
  G4double nuBrRatio = 0.;
  G4double eBrRatio = 0.;
  G4double muBrRatio = 0.;
  G4double hBrRatio = 0.;
  G4int IDPDG = 5500022; // Totally invisible A' PDG ID, can be redefined below for different decays
                         // https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
  if(DecayType) {
    isStable = false;
    if (BranchingType == 0) { // X boson with visible decays
      if(MassIn < 2.001*electron_mass_c2) isStable = true;
      WidthIn =
        isStable ? 0 : (1./3.)*CLHEP::fine_structure_const*MassIn*epsilIn*epsilIn*sqrt(1.-4.*RatioEA2)*(1.+2.*RatioEA2);
      IDPDG = 5500122;
    } else if (BranchingType == 1) { // B-L Z' boson with coupling to all SM particles
      if(MassIn > 600.) {G4cout << "Branching ratios for this BranchingType and this mass are not yet implemented, exiting" << G4endl; exit(1);}
      nuWidth = epsilIn*epsilIn*CLHEP::fine_structure_const*MassIn;
      if(MassIn > 2.*electron_mass_c2) eWidth = (1./3.)*CLHEP::fine_structure_const*MassIn*epsilIn*epsilIn*sqrt(1.-4.*RatioEA2)*(1.+2.*RatioEA2);
      if(MassIn > 2.*muMass) muWidth = (1./3.)*CLHEP::fine_structure_const*MassIn*epsilIn*epsilIn*sqrt(1.-4.*massRatio2)*(1.+2.*massRatio2);
      if(MassIn > pi0Mass) {
        hWidth = (CLHEP::fine_structure_const*fine_structure_const*epsilIn*epsilIn*MassIn*MassIn*MassIn) /
                 (96.*3.141*3.141*3.141*0.93*0.93*pi0Mass*pi0Mass);
        hWidth *= (1. - pi0Mass*pi0Mass/(MassIn*MassIn));
        G4double a = 1. - (MassIn*MassIn)/(782.66*782.66);
        G4double b = 12.3/782.66;
        G4double mod2 = 1./(a*a + b*b);
        hWidth *= mod2;
      }
      WidthIn = nuWidth + eWidth + muWidth + hWidth;
      nuBrRatio = nuWidth/WidthIn;
      eBrRatio = eWidth/WidthIn;
      muBrRatio = muWidth/WidthIn;
      hBrRatio = hWidth/WidthIn;
      IDPDG = 5500222;
    } else {
      G4cout << "BranchingType = " << BranchingType << " is not implemented, exiting" << G4endl;
      exit(1);
    }
  }

  if( !anInstance ) {
    anInstance = new G4ParticleDefinition(
        /* Name ..................... */ name,
        /* Mass ..................... */ MassIn,
        /* Decay width .............. */ WidthIn,
        /* Charge ................... */ 0.*eplus,
        /* 2*spin ................... */ 2,
        /* parity ................... */ +1,
        /* C-conjugation ............ */ 0,
        /* 2*Isospin ................ */ 0,
        /* 2*Isospin3 ............... */ 0,
        /* G-parity ................. */ 0,
        /* type ..................... */ "boson",
        /* lepton number ............ */ 0,
        /* baryon number ............ */ 0,
        /* PDG encoding ............. */ IDPDG,
        /* stable ................... */ isStable,
        /* lifetime.................. */ 0,
        /* decay table .............. */ NULL,
        /* shortlived ............... */ false,
        /* subType .................. */ "DMParticleAPrime",
        /* anti particle encoding ... */ IDPDG
          );

    if(!isStable)
    {
      // Life time is given from width
      ((DMParticle*)anInstance)->CalculateLifeTime();

      //create Decay Table
      G4DecayTable* table = new G4DecayTable();

      if (BranchingType == 0) { // X boson with visible decays

        // create a decay channel
        // X -> e+ + e-
        G4VDecayChannel* mode = new G4PhaseSpaceDecayChannel(name, 1., 2, "e-", "e+");

        table->Insert(mode);
      }

      if (BranchingType == 1) { // B-L Z' boson with coupling to all SM particles

        G4VDecayChannel** mode = new G4VDecayChannel*[6];
        // DMParticleZPrime -> nu_e + anti_nu_e
        mode[0] = new G4PhaseSpaceDecayChannel(name, nuBrRatio/3., 2, "anti_nu_e", "nu_e");
        // DMParticleZPrime -> nu_mu + anti_nu_mu
        mode[1] = new G4PhaseSpaceDecayChannel(name, nuBrRatio/3., 2, "anti_nu_mu", "nu_mu");
        // DMParticleZPrime -> nu_tau + anti_nu_tau
        mode[2] = new G4PhaseSpaceDecayChannel(name, nuBrRatio/3., 2, "anti_nu_tau", "nu_tau");
        // DMParticleZPrime -> e+ + e-
        mode[3] = new G4PhaseSpaceDecayChannel(name, eBrRatio, 2, "e+", "e-");
        // DMParticleZPrime -> mu+ + mu-
        mode[4] = new G4PhaseSpaceDecayChannel(name, muBrRatio, 2, "mu+", "mu-");
        // DMParticleZPrime -> pi0 + gamma
        mode[5] = new G4PhaseSpaceDecayChannel(name, hBrRatio, 2, "pi0", "gamma");

       for (G4int index = 0; index < 6; index++) table->Insert(mode[index]);
       delete [] mode;
      }

      anInstance->SetDecayTable(table);
      anInstance->DumpTable();
    }
  }
  theInstance = reinterpret_cast<DMParticleAPrime*>(anInstance);
  return theInstance;
}
