#include "DMParticleChi1.hh"
#include "DarkMatterParametersRegistry.hh"

#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DecayTable.hh"

DMParticleChi1 * DMParticleChi1::theInstance = nullptr;

DMParticleChi1* DMParticleChi1::Definition()
{
  if( theInstance ) {
    return theInstance;
  }
  //get parameters from registry (NOTE: mass is parsed in GeV)
  DarkMatterParametersRegistry* DMpar = DarkMatterParametersRegistry::GetInstance();
  double MassChi1 = (DMpar->GetRegisteredParam("DMMass")) * DMpar->GetRegisteredParam("RDM", 1./3.);

  const G4String name = "DMParticleChi1";
  // search in particle table
  G4ParticleTable * pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition * anInstance = pTable->FindParticle(name);

  if( !anInstance ) {
    anInstance = new G4ParticleDefinition(
        /* Name ..................... */ name,
        /* Mass ..................... */ MassChi1,
        /* Decay width .............. */ 0.,
        /* Charge ................... */ 0.*eplus,
        /* 2*spin ................... */ 1,
        /* parity ................... */ +1,
        /* C-conjugation ............ */ 0,
        /* 2*Isospin ................ */ 0,
        /* 2*Isospin3 ............... */ 0,
        /* G-parity ................. */ 0,
        /* type ..................... */ "fermion",
        /* lepton number ............ */ 0,
        /* baryon number ............ */ 0,
        /* PDG encoding ............. */ 5200014, // https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
        /* stable ................... */ true,
        /* lifetime.................. */ 0,
        /* decay table .............. */ NULL,
        /* shortlived ............... */ false,
        /* subType .................. */ "DMParticle",
        /* anti particle encoding ... */ -5200014
          );
  }
  theInstance = reinterpret_cast<DMParticleChi1*>(anInstance);
  G4cout << "The particle: " << theInstance->GetParticleName() << " mass in GeV is: " << theInstance->GetPDGMass()/GeV << "\n";
  return theInstance;
}
