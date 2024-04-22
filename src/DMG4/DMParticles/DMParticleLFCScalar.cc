#include "DMParticleLFCScalar.hh"
#include "DarkMatterParametersFactory.hh"

#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4PhysicalConstants.hh"
#include "G4DalitzDecayChannel.hh"
#include "G4DecayTable.hh"
#include "G4Electron.hh"
#include "G4MuonMinus.hh"

DMParticleLFCScalar * DMParticleLFCScalar::theInstance = nullptr;

DMParticleLFCScalar* DMParticleLFCScalar::Definition()
{
  if( theInstance ) {
    return theInstance;
  }
  //get parameters from factory (NOTE: mass is parsed in GeV)
  DarkMatterParametersFactory* DMpar = DarkMatterParametersFactory::GetInstance();
  G4double MassIn    = DMpar->GetRegisteredParam("DMMass");
  G4double epsilIn   = DMpar->GetRegisteredParam("Epsilon");
  G4double DecayType = DMpar->GetRegisteredParam("DecayType");
  G4double BranchingType = DMpar->GetRegisteredParam("BranchingType", 0);

  G4String name = "DMParticleLFCScalar";
  const G4String nameSubType = "DMParticleLFCScalar";
  // search in particle table
  G4ParticleTable * pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition * anInstance = pTable->FindParticle(name);
  const G4double muMass = G4MuonMinus::MuonMinusDefinition()->GetPDGMass();
  const G4double elMass = G4Electron::Definition()->GetPDGMass();
  G4double RatioEA2 = muMass*muMass/(MassIn*MassIn);
  G4bool isStable = true;
  G4double WidthIn = 0.; // in MeV
  G4int IDPDG = 5400023; // https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf

  // compute width
  // TODO: THIS PART IS NOT IMPLEMENTED YET!
  if(DecayType) {
    isStable = false;
    if (BranchingType == 0) {
      IDPDG = 5400123; // https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
      name = "DMParticleLFCScalar";
      if(MassIn > 2.001*muMass) {
        WidthIn = (1./2.)*CLHEP::fine_structure_const*MassIn*epsilIn*epsilIn*sqrt(1.-4.*RatioEA2)*(1.-4.*RatioEA2);
      }
    } else if (BranchingType == 1) {
      IDPDG = 5400199; // https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
      if (MassIn > elMass+muMass) {
        // width is implemented according to https://arxiv.org/abs/2202.04410
        G4double p = sqrt(MassIn*MassIn-pow(elMass+muMass,2.)*(MassIn*MassIn-pow(elMass-muMass,2.)))/MassIn*0.5;
        WidthIn = epsilIn*epsilIn*p/(8*M_PI)*(1.0-(elMass*elMass+muMass*muMass)/(MassIn*MassIn));
      }
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
        /* 2*spin ................... */ 0,
        /* parity ................... */ 0,
        /* C-conjugation ............ */ 0,
        /* 2*Isospin ................ */ 0,
        /* 2*Isospin3 ............... */ 0,
        /* G-parity ................. */ 0,
        /* type ..................... */ "boson",
        /* lepton number ............ */ 0,
        /* baryon number ............ */ 0,
        /* PDG encoding ............. */ IDPDG, // https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
        /* stable ................... */ isStable,
        /* lifetime.................. */ 0,
        /* decay table .............. */ NULL,
        /* shortlived ............... */ false,
        /* subType .................. */ nameSubType,
        /* anti particle encoding ... */ IDPDG
          );

    if(!isStable)
    {
      // Life time is given from width
      ((DMParticle*)anInstance)->CalculateLifeTime();

      //create Decay Table
      G4DecayTable* table = new G4DecayTable();

      // Define decay channels according to branching type
      if (BranchingType == 0) {
        // X -> mu- + e+
        G4VDecayChannel* mode = new G4PhaseSpaceDecayChannel(name, 1, 2, "mu-", "e+");
        table->Insert(mode);
        delete mode;
      } else if (BranchingType == 1) {
        // X -> mu- + e+
        G4VDecayChannel* mode = new G4PhaseSpaceDecayChannel(name, 1, 2, "mu-", "e+");
        table->Insert(mode);
        delete mode;
      } else if (BranchingType == 2) {
        // X -> tau- + mu+
        G4VDecayChannel* mode = new G4PhaseSpaceDecayChannel(name, 1., 2, "tau-", "mu+");
        table->Insert(mode);
        delete mode;
      } else {
        G4cout << "BranchingType = " << BranchingType << " is not implemented, exiting" << G4endl;
        exit(1);
      }
      anInstance->SetDecayTable(table);
      anInstance->DumpTable();
    }

  }
  theInstance = reinterpret_cast<DMParticleLFCScalar*>(anInstance);
  return theInstance;
}
