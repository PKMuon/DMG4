#include "DMParticleChi2.hh"
#include "DarkMatterParametersFactory.hh"

#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DecayTable.hh"

DMParticleChi2 * DMParticleChi2::theInstance = nullptr;

DMParticleChi2* DMParticleChi2::Definition()
{
  if( theInstance ) {
    return theInstance;
  }
  //get parameters from factory (NOTE: mass is parsed in GeV)
  DarkMatterParametersFactory* DMpar = DarkMatterParametersFactory::GetInstance();
  const G4double DMMass = DMpar->GetRegisteredParam("DMMass");
  const G4double MassChi1 = (DMpar->GetRegisteredParam("DMMass")) * DMpar->GetRegisteredParam("RDM", 1./3.);
  const G4double MassChi2 = (1. + DMpar->GetRegisteredParam("Ffactor", 0.4)) * MassChi1;

  const G4String name = "DMParticleChi2";
  // search in particle table
  G4ParticleTable * pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition * anInstance = pTable->FindParticle(name);

  // parameters for the width
  G4double Epsilon = DMpar->GetRegisteredParam("Epsilon");
  G4double AlphaD = DMpar->GetRegisteredParam("AlphaD");
  G4double Theta = DMpar->GetRegisteredParam("Theta");
  G4double Splitting = MassChi2 - MassChi1;

  std::cout << "Mass Chi2: " << MassChi2 <<std::endl;
  std::cout << "Mass Chi1: " << MassChi1 <<std::endl;
  std::cout << "Splitting: " << Splitting<<std::endl;
  // calculate the width
  const G4double K = 0.640;
  G4double WidthIn = K*4.*Epsilon*Epsilon*fine_structure_const*AlphaD*pow(Splitting,5.)/(15.*CLHEP::pi*pow(DMMass,4.));
  std::cout << "Mass Chi2: " << MassChi2/MeV << " MeV " <<std::endl;
  std::cout << "===> Width Chi2->chi1ee " << WidthIn/MeV <<" MeV "<< std::endl;

  if( !anInstance ) {
    anInstance = new G4ParticleDefinition(
        /* Name ..................... */ name,
        /* Mass ..................... */ MassChi2,
        /* Decay width .............. */ WidthIn,
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
        /* PDG encoding ............. */ 5200013, // https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
        /* stable ................... */ false,
        /* lifetime.................. */ 0,
        /* decay table .............. */ NULL,
        /* shortlived ............... */ false,
        /* subType .................. */ "DMParticle",
        /* anti particle encoding ... */ -5200013
          );

    // Life time is given from width
    ((DMParticle*)anInstance)->CalculateLifeTime();

    double lifetime =((DMParticle*)anInstance)->GetPDGLifeTime();
    std::cout << "===> Lifetime Chi2 " << lifetime/s <<" s, c*tau: "<<(lifetime*CLHEP::c_light)/cm <<" cm "<< std::endl;

    // create decay table and add modes
    G4DecayTable* table = new G4DecayTable();
    G4VDecayChannel** mode = new G4VDecayChannel*[1];
    mode[0] = new G4PhaseSpaceDecayChannel(name, 1., 3, "DMParticleChi1", "e+", "e-");
    for (G4int index = 0; index < 1; index++) table->Insert(mode[index]);
    delete [] mode;

    anInstance->SetDecayTable(table);
    anInstance->DumpTable();

  }
  theInstance = reinterpret_cast<DMParticleChi2*>(anInstance);
  G4cout << "The particle: " << theInstance->GetParticleName() << " mass is: " << theInstance->GetPDGMass()/GeV << "\n";
  return theInstance;
}
