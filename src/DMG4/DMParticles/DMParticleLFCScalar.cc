#include "DMParticleLFCScalar.hh"
#include "DarkMatterParametersFactory.hh"

#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4PhysicalConstants.hh"
#include "G4DalitzDecayChannel.hh"
#include "G4DecayTable.hh"
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

  G4String name = "DMParticleLFCScalar";
  const G4String nameSubType = "DMParticleLFCScalar";
  // search in particle table
  G4ParticleTable * pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition * anInstance = pTable->FindParticle(name);
  const G4double muMass = G4MuonMinus::MuonMinusDefinition()->GetPDGMass();
  G4double RatioEA2 = muMass*muMass/(MassIn*MassIn);
  G4bool isStable = DecayType > 0 ? false : true;
  if(MassIn < 2.001*muMass*GeV) isStable = true;

  G4int IDPDG = 5400023; // https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
  G4double WidthIn = 0.;
  if(!isStable) {
    WidthIn = (1./2.)*CLHEP::fine_structure_const*MassIn*epsilIn*epsilIn*sqrt(1.-4.*RatioEA2)*(1.-4.*RatioEA2);
    IDPDG = 5400123;
    name = "DMParticleLFCScalar";
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
        /* PDG encoding ............. */ IDPDG,
        /* stable ................... */ isStable,
        /* lifetime.................. */ 0,
        /* decay table .............. */ NULL,
        /* shortlived ............... */ false,
        /* subType .................. */ nameSubType,
        /* anti particle encoding ... */ IDPDG
          );
  }
  if(!isStable)
  {
    // Life time is given from width
    ((DMParticle*)anInstance)->CalculateLifeTime();

    //create Decay Table
    G4DecayTable* table = new G4DecayTable();

    // create a decay channel
    // X -> mu+ + mu- + anti_nu_mu + nu_tau
    G4VDecayChannel* mode = new G4PhaseSpaceDecayChannel(name, 1., 2, "tau-", "mu+");

    table->Insert(mode);
    anInstance->SetDecayTable(table);
  }
  theInstance = reinterpret_cast<DMParticleLFCScalar*>(anInstance);
  return theInstance;
}
