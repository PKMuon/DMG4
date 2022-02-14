#include "DMParticleZPrime.hh"
#include "DarkMatterParametersFactory.hh"

#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DalitzDecayChannel.hh"
#include "G4DecayTable.hh"
#include "G4MuonMinus.hh"

DMParticleZPrime * DMParticleZPrime::theInstance = nullptr;

DMParticleZPrime* DMParticleZPrime::Definition()
{
  if( theInstance ) {
    return theInstance;
  }
  //get parameters from factory (NOTE: mass is parsed in GeV)
  DarkMatterParametersFactory* DMpar = DarkMatterParametersFactory::GetInstance();
  G4double MassIn    = DMpar->GetRegisteredParam("MassA")*GeV;
  G4double epsilIn   = DMpar->GetRegisteredParam("Epsil");
  G4double DecayType = DMpar->GetRegisteredParam("DecayType");

  const G4String name = "DMParticleZPrime";
  // search in particle table]
  G4ParticleTable * pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition * anInstance = pTable->FindParticle(name);

  // determine decay rates
  const G4double muMass = G4MuonMinus::MuonMinusDefinition()->GetPDGMass()/MeV;
  G4double massRatio2 = muMass*muMass/(MassIn*MassIn);
  G4double WidthIn = 0.; // in MeV
  G4double nuWidth = 0.;
  G4double muWidth = 0.;
  G4double nuBrRatio = 0.;
  G4double muBrRatio = 0.;
  nuWidth = epsilIn*epsilIn*CLHEP::fine_structure_const*(1./3.)*MassIn;
  if (MassIn > 2.*muMass) {
    G4double factor = (1.+2.*massRatio2)*sqrt(1.-4.*massRatio2);
    muWidth = nuWidth*factor;
    muBrRatio = factor/(1.+factor);
  }
  nuBrRatio = 1. - muBrRatio;
  G4bool isStable = DecayType > 0 ? false : true;
  WidthIn = isStable ? 0 : nuWidth + muWidth;

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
        /* PDG encoding ............. */ 5500023, // https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
        /* stable ................... */ isStable,
        /* lifetime.................. */ 0,
        /* decay table .............. */ NULL,
        /* shortlived ............... */ false,
        /* subType .................. */ "DMParticleZPrime",
        /* anti particle encoding ... */ 5500023
          );
    if(!isStable)
    {
      // Life time is given from width
      ((DMParticle*)anInstance)->CalculateLifeTime();

      // create decay table and add modes
      G4DecayTable* table = new G4DecayTable();
      G4VDecayChannel** mode = new G4VDecayChannel*[3];
      // DMParticleZPrime -> nu_mu + anti_nu_mu
      mode[0] = new G4PhaseSpaceDecayChannel(name, nuBrRatio/2., 2, "anti_nu_mu", "nu_mu");
      // DMParticleZPrime -> nu_tau + anti_nu_tau
      mode[1] = new G4PhaseSpaceDecayChannel(name, nuBrRatio/2., 2, "anti_nu_tau", "nu_tau");
      // DMParticleZPrime -> mu+ + mu-
      mode[2] = new G4PhaseSpaceDecayChannel(name, muBrRatio, 2, "mu+", "mu-");
      for (G4int index = 0; index < 3; index++) table->Insert(mode[index]);
      delete [] mode;

      anInstance->SetDecayTable(table);
      //anInstance->DumpTable();
    }

  }
  theInstance = reinterpret_cast<DMParticleZPrime*>(anInstance);
  return theInstance;
}
