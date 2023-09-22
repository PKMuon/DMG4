#include "DMParticleScalarLFC.hh"
#include "DarkMatterParametersFactory.hh"

#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DalitzDecayChannel.hh"
#include "G4DecayTable.hh"
#include "G4Electron.hh"
#include "G4MuonMinus.hh"

DMParticleScalarLFC * DMParticleScalarLFC::theInstance = nullptr;

DMParticleScalarLFC* DMParticleScalarLFC::Definition()
{
  if( theInstance ) {
    return theInstance;
  }
  //get parameters from factory (NOTE: mass is parsed in GeV)
  DarkMatterParametersFactory* DMpar = DarkMatterParametersFactory::GetInstance();
  G4double MassIn    = DMpar->GetRegisteredParam("DMMass");
  G4double epsilIn   = DMpar->GetRegisteredParam("Epsilon");
  G4double DecayType = DMpar->GetRegisteredParam("DecayType"); 

  const G4String name = "DMParticleScalarLFC";
  // search in particle table]
  G4ParticleTable * pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition * anInstance = pTable->FindParticle(name);

  G4int IDPDG = 5400099;
  G4double Mel = G4Electron::Definition()->GetPDGMass()/GeV;
  G4double Mmu = G4MuonMinus::Definition()->GetPDGMass()/GeV;
  if (MassIn < Mel+Mmu) {
    G4cout << "Invalid mass for " << name << ", aborting..." << G4endl;
    exit(1);
  }

  // compute width
  // width is implemented according to https://arxiv.org/abs/2202.04410
  G4double WidthIn = 0.;
  G4bool isStable = DecayType > 0 ? false : true;
  if (DecayType > 0) IDPDG = 5400199;
  if (MassIn > Mel+Mmu) {
    G4double p = sqrt(MassIn*MassIn-pow(Mel+Mmu,2.)*(MassIn*MassIn-pow(Mel-Mmu,2.)))/MassIn*0.5;
    WidthIn = isStable ? 0. : epsilIn*epsilIn*p/(8*M_PI)*(1.0-(Mel*Mel+Mmu*Mmu)/(MassIn*MassIn));
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
        /* subType .................. */ "DMParticleScalarLFC",
        /* anti particle encoding ... */ IDPDG
          );

    if (MassIn > Mel+Mmu && !isStable) {
      // set the decay channels
      ((DMParticle*)anInstance)->CalculateLifeTime();
      // create decay table and add modes
      G4DecayTable* table = new G4DecayTable();
      G4VDecayChannel** mode = new G4VDecayChannel*[1];
      mode[0] = new G4PhaseSpaceDecayChannel(name, 1, 2, "mu-", "e+");
      for (G4int index = 0; index < 1; index++) table->Insert(mode[index]);
      delete [] mode;

      anInstance->SetDecayTable(table);
      anInstance->DumpTable();
    }
    
  }
  theInstance = reinterpret_cast<DMParticleScalarLFC*>(anInstance);
  return theInstance;
}
