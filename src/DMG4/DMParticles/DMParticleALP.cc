#include "DMParticleALP.hh"
#include "DarkMatterParametersRegistry.hh"

#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DalitzDecayChannel.hh"
#include "G4DecayTable.hh"

DMParticleALP * DMParticleALP::theInstance = nullptr;

DMParticleALP* DMParticleALP::Definition()
{
  if( theInstance ) {
    return theInstance;
  }
  //get parameters from registry (NOTE: mass is parsed in GeV)
  DarkMatterParametersRegistry* DMpar = DarkMatterParametersRegistry::GetInstance();
  G4double MassIn    = DMpar->GetRegisteredParam("DMMass");
  G4double epsilIn   = DMpar->GetRegisteredParam("Epsilon");
  G4double DecayType = DMpar->GetRegisteredParam("DecayType");

  G4String name = "DMParticleDarkALP";
  // search in particle table]
  G4ParticleTable * pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition * anInstance = pTable->FindParticle(name);

  G4bool isStable = DecayType > 0 ? false : true;
  G4int IDPDG = 5300022; // https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
  G4double WidthIn = 0.;
  if(!isStable) {
    WidthIn = 1./(64.*pi)*MassIn*MassIn*MassIn*epsilIn*epsilIn;
    IDPDG = 5300122;
    name = "DMParticleALP";
  }
  if( !anInstance ) {
    anInstance = new G4ParticleDefinition(
        /* Name ..................... */ name,
        /* Mass ..................... */ MassIn,
        /* Decay width .............. */ WidthIn,
        /* Charge ................... */ 0.,
        /* 2*spin ................... */ 0,
        /* parity ................... */ -1,
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
        /* subType .................. */ "DMParticleALP",
        /* anti particle encoding ... */ IDPDG
          );

    if(!isStable)
    {
      // Life time is given from width
      ((DMParticle*)anInstance)->CalculateLifeTime();

      //create Decay Table
      G4DecayTable* table = new G4DecayTable();

      // create a decay channel
      G4VDecayChannel* mode;
      // ALP -> gamma + gamma
      mode = new G4PhaseSpaceDecayChannel(name, 1., 2, "gamma", "gamma");
      table->Insert(mode);

      anInstance->SetDecayTable(table);
    }
  }
  theInstance = reinterpret_cast<DMParticleALP*>(anInstance);
  return theInstance;
}
