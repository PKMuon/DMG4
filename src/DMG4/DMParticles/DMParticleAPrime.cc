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
  G4double MassIn    = DMpar->GetRegisteredParam("DMMass");
  G4double epsilIn   = DMpar->GetRegisteredParam("Epsilon");
  G4double DecayType = DMpar->GetRegisteredParam("DecayType");
  G4double BranchingType = DMpar->GetRegisteredParam("BranchingType", 0);

  G4String name = "DMParticleAPrime";
  const G4String nameSubType = "DMParticleAPrime";
  // search in particle table
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
  G4double Chi12Width = 0.;
  G4double Chi11Width = 0.;
  G4double Chi22Width = 0.;
  G4double nuBrRatio = 0.;
  G4double eBrRatio = 0.;
  G4double muBrRatio = 0.;
  G4double hBrRatio = 0.;
  G4double Chi12BrRatio = 0.;
  G4double Chi11BrRatio = 0.;
  G4double Chi22BrRatio = 0.;
  G4int IDPDG = 5500022; // Totally invisible A' PDG ID, can be redefined below for different decays
                         // https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf
  if(DecayType) {
    isStable = false;
    if (BranchingType == 0) { // X boson with visible decays
      if(MassIn < 2.001*electron_mass_c2) isStable = true;
      if(!isStable) {
        WidthIn = (1./3.)*CLHEP::fine_structure_const*MassIn*epsilIn*epsilIn*sqrt(1.-4.*RatioEA2)*(1.+2.*RatioEA2);
        IDPDG = 5500122;
        name = "DMParticleXBoson";
      }
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
      name = "DMParticleB-LBoson";
    } else if (BranchingType == 2) { // Inelastic DM: decay to Chi2 + Chi1
      const G4double MChi1 = (DMpar->GetRegisteredParam("DMMass")) * DMpar->GetRegisteredParam("RDM", 1./3.);
      const G4double MChi2 = (1. + DMpar->GetRegisteredParam("Ffactor", 0.4)) * MChi1;
      const G4double AlphaD = DMpar->GetRegisteredParam("AlphaD");
      const G4double Delta = MChi2 - MChi1;
      if(MassIn > 2.*electron_mass_c2) eWidth = (1./3.)*CLHEP::fine_structure_const*MassIn*epsilIn*epsilIn*sqrt(1.-4.*RatioEA2)*(1.+2.*RatioEA2);
      if(MassIn > MChi1+MChi2) {
        Chi12Width = AlphaD*MassIn/6.*sqrt(1.+MChi1*MChi1/(MassIn*MassIn)*(Delta*Delta/(MassIn*MassIn)*(Delta/MChi1+2.)*(Delta/MChi1+2.)-2.*
                     (1.+(1.+Delta/MChi1)*(1.+Delta/MChi1))))*(2.-MChi1*MChi1/(MassIn*MassIn)*(1.+Delta*Delta/(MassIn*MassIn)*(2.+Delta/MChi1)*
                     (2.+Delta/MChi1)+(1.+Delta/MChi1)*(1.+Delta/MChi1)-6.*(1.+Delta/MChi1)));
      }
      WidthIn = eWidth + Chi12Width;
      if(WidthIn == 0.) isStable = true;
      eBrRatio = eWidth/WidthIn;
      Chi12BrRatio = Chi12Width/WidthIn;
      IDPDG = 5500322;
      name = "DMParticleInelasticBoson";
    } else if (BranchingType == 3) { // Dirac inelastic DM: decay to Chi2 + Chi1, Chi1 + Chi1, Chi2 + Chi2                                                                  
      const G4double MChi1 = (DMpar->GetRegisteredParam("DMMass")) * DMpar->GetRegisteredParam("RDM", 1./3.);
      const G4double MChi2 = (1. + DMpar->GetRegisteredParam("Ffactor", 0.4)) * MChi1;
      const G4double AlphaD = DMpar->GetRegisteredParam("AlphaD");
      const G4double Delta = MChi2 - MChi1;
      const G4double Theta = DMpar->GetRegisteredParam("Theta");
      if(MassIn > 2.*electron_mass_c2) eWidth = (1./3.)*CLHEP::fine_structure_const*MassIn*epsilIn*epsilIn*sqrt(1.-4.*RatioEA2)*(1.+2.*RatioEA2);
      if (MassIn > MChi1+MChi2) {
        Chi12Width = AlphaD*MassIn/6.*sqrt(1.+MChi1*MChi1/(MassIn*MassIn)*(Delta*Delta/(MassIn*MassIn)*(Delta/MChi1+2.)*(Delta/MChi1+2.)-2.*
                    (1.+(1.+Delta/MChi1)*(1.+Delta/MChi1))))*(2.-MChi1*MChi1/(MassIn*MassIn)*(1.+Delta*Delta/(MassIn*MassIn)*(2.+Delta/MChi1)*
                    (2.+Delta/MChi1)+(1.+Delta/MChi1)*(1.+Delta/MChi1)-6.*(1.+Delta/MChi1)));
      }
      if (MassIn > 2.*MChi1) {
        Chi11Width = AlphaD/3.*MassIn*(1.+2*MChi1*MChi1/(MassIn*MassIn))*sqrt(1.-4*MChi1*MChi1/(MassIn*MassIn));
      }
      if (MassIn > 2.*MChi2) {
        Chi22Width = AlphaD/3.*MassIn*(1.+2*MChi2*MChi2/(MassIn*MassIn))*sqrt(1.-4*MChi2*MChi2/(MassIn*MassIn));
      }
      WidthIn = eWidth + pow(sin(2*Theta),2.)*Chi12Width+pow(sin(Theta),4.)*Chi11Width+pow(cos(Theta),4.)*Chi22Width;
      if(WidthIn == 0.) isStable = true;
      // branching ratio
      eBrRatio = eWidth/WidthIn;
      Chi12BrRatio = Chi12Width/WidthIn;
      Chi11BrRatio = Chi11Width/WidthIn;
      Chi22BrRatio = Chi22Width/WidthIn;
      IDPDG = 5500322;
      name = "DMParticleInelasticBoson";

      std::cout << "===> Width dark photon to chi2chi2 " << Chi22Width << std::endl;
      std::cout << "===> Width dark photon to chi1chi1 " << Chi11Width << std::endl;
      std::cout << "===> Width dark photon to chi1chi2 " << Chi12Width << std::endl;

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
        /* subType .................. */ nameSubType,
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

      if (BranchingType == 2) { // Inelastic DM: decay to Chi2 + Chi1

        G4VDecayChannel** mode = new G4VDecayChannel*[2];
        // DMParticleZPrime -> e+ + e-
        mode[0] = new G4PhaseSpaceDecayChannel(name, eBrRatio, 2, "e+", "e-");
        mode[1] = new G4PhaseSpaceDecayChannel(name, Chi12BrRatio, 2, "DMParticleChi1", "DMParticleChi2");
        for (G4int index = 0; index < 2; index++) table->Insert(mode[index]);
        delete [] mode;
      }

      if (BranchingType == 3) { // Inelastic DM: decay to Chi1 + Chi1, Chi1 + Chi2, Chi2 + Chi2
        // TODO: implement the correct branching ratio, for now we only assume
        // that A'->chi1+chi2, so Br(A' -> everything else) = 0.

        G4VDecayChannel** mode = new G4VDecayChannel*[4];
        // DMParticleZPrime -> e+ + e-
        mode[0] = new G4PhaseSpaceDecayChannel(name, eBrRatio, 2, "e+", "e-");
        mode[1] = new G4PhaseSpaceDecayChannel(name, Chi12BrRatio, 2, "DMParticleChi1", "DMParticleChi2");
        mode[2] = new G4PhaseSpaceDecayChannel(name, Chi11BrRatio, 2, "DMParticleChi1", "DMParticleChi1");
        mode[3] = new G4PhaseSpaceDecayChannel(name, Chi22BrRatio, 2, "DMParticleChi2", "DMParticleChi2");
        for (G4int index = 0; index < 4; index++) table->Insert(mode[index]);
        delete [] mode;
      }

      anInstance->SetDecayTable(table);
      anInstance->DumpTable();
    }
  }
  theInstance = reinterpret_cast<DMParticleAPrime*>(anInstance);
  return theInstance;
}
