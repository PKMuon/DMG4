#include "DMParticleAxial.hh"
#include "DarkMatterParametersFactory.hh"

#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DalitzDecayChannel.hh"
#include "G4DecayTable.hh"

DMParticleAxial * DMParticleAxial::theInstance = nullptr;

DMParticleAxial* DMParticleAxial::Definition()
{
  if( theInstance ) {
    return theInstance;
  }
  //get parameters from factory (NOTE: mass is parsed in GeV)
  DarkMatterParametersFactory* DMpar = DarkMatterParametersFactory::GetInstance();
  G4double MassIn    = DMpar->GetRegisteredParam("DMMass")*GeV;
  G4double epsilIn   = DMpar->GetRegisteredParam("Epsilon");
  G4double DecayType = DMpar->GetRegisteredParam("DecayType");

  const G4String name = "DMParticleAxial";
  // search in particle table]
  G4ParticleTable * pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition * anInstance = pTable->FindParticle(name);
  G4double RatioEA2 = electron_mass_c2*electron_mass_c2/(MassIn*MassIn);
  G4bool isStable = DecayType > 0 ? false : true;
  G4double WidthIn =
    isStable ? 0 : (1./3.)*CLHEP::fine_structure_const*MassIn*epsilIn*epsilIn*sqrt(1.-4.*RatioEA2)*(1.-4.*RatioEA2);
  if( !anInstance ) {
    anInstance = new G4ParticleDefinition(
        /* Name ..................... */ name,
        /* Mass ..................... */ MassIn,
        /* Decay width .............. */ WidthIn,
        /* Charge ................... */ 0.*eplus,
        /* 2*spin ................... */ 2,
        /* parity ................... */ -1,
        /* C-conjugation ............ */ 0,
        /* 2*Isospin ................ */ 0,
        /* 2*Isospin3 ............... */ 0,
        /* G-parity ................. */ 0,
        /* type ..................... */ "boson",
        /* lepton number ............ */ 0,
        /* baryon number ............ */ 0,
        /* PDG encoding ............. */ 5510022, // https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
        /* stable ................... */ true,
        /* lifetime.................. */ 0,
        /* decay table .............. */ NULL,
        /* shortlived ............... */ false,
        /* subType .................. */ "DMParticleAxial",
        /* anti particle encoding ... */ 5510022
          );
  }
  theInstance = reinterpret_cast<DMParticleAxial*>(anInstance);
  return theInstance;
}
