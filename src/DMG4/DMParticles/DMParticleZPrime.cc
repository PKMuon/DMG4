#include "DMParticleZPrime.hh"
#include "DarkMatterParametersRegistry.hh"

#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4DalitzDecayChannel.hh"
#include "G4DecayTable.hh"
#include "G4MuonMinus.hh"
#include "G4PionZero.hh"

DMParticleZPrime * DMParticleZPrime::theInstance = nullptr;

DMParticleZPrime* DMParticleZPrime::Definition()
{
  if( theInstance ) {
    return theInstance;
  }

  // Function to calculate widths for dark photon
  auto ZPrimeWidth = [](double mass1, double mass2, double MassIn) {
    // Reduce the formula to the case where the two masses are equal
    if (mass1 == mass2)
      return MassIn/3.*
        sqrt(1.-4.*mass1*mass1/(MassIn*MassIn))*
        (1.+2.*mass1*mass1/(MassIn*MassIn));
    // Full expression for the general case
    return MassIn/3.*
      sqrt((1.-pow((mass1+mass2),2.)/(MassIn*MassIn))*(1.-pow((mass1-mass2),2.)/(MassIn*MassIn)))*
      (1.-(pow((mass1-mass2),2.)-4.*mass1*mass2)/(2.*MassIn*MassIn)-(pow(mass1*mass1-mass2*mass2,2.)/(2.*pow(MassIn,4.))));
  };

  //get parameters from registry (NOTE: mass is parsed in GeV)
  DarkMatterParametersRegistry* DMpar = DarkMatterParametersRegistry::GetInstance();
  G4double MassIn    = DMpar->GetRegisteredParam("DMMass");
  G4double epsilIn   = DMpar->GetRegisteredParam("Epsilon");
  G4double alphaD = DMpar->GetRegisteredParam("AlphaD", 0.5);
  G4double DecayType = DMpar->GetRegisteredParam("DecayType");
  G4double BranchingType = DMpar->GetRegisteredParam("BranchingType", 0);
  const G4double mChi = (DMpar->GetRegisteredParam("DMMass")) * DMpar->GetRegisteredParam("RDM", 1./3.);

  G4String name = "DMParticleZPrime";
  const G4String nameSubType = "DMParticleZPrime";
  // search in particle table
  G4ParticleTable * pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition * anInstance = pTable->FindParticle(name);

  // determine decay rates
  const G4double muMass = G4MuonMinus::MuonMinusDefinition()->GetPDGMass();
  const G4double pi0Mass = G4PionZero::PionZeroDefinition()->GetPDGMass();
  const G4double omegaMass = 782.66*MeV; // omega meson mass from pdg:  https://pdglive.lbl.gov/Particle.action?init=0&node=M001&home=MXXX005
  G4bool isStable = true;
  G4double massRatio2 = muMass*muMass/(MassIn*MassIn);
  G4double WidthIn = 0.; // in MeV
  G4double nuWidth = 0.;
  G4double eWidth = 0.;
  G4double muWidth = 0.;
  G4double chiWidth = 0.;
  G4double hWidth = 0.;
  G4double nuBrRatio = 0.;
  G4double eBrRatio = 0.;
  G4double muBrRatio = 0.;
  G4double hBrRatio = 0.;
  G4double chiBrRatio = 0.;
  G4int IDPDG = 5500023; // https://pdg.lbl.gov/2019/reviews/rpp2019-rev-monte-carlo-numbering.pdf

  if(DecayType) {
    isStable = false;
    // Lmu-Ltau, vanilla, neutrino decay
    if (BranchingType == 0) {
      // Change name accordingly
      IDPDG = 5500123;
      name = "DMParticleZPrimeVis";
      // Calculate widths
      nuWidth = epsilIn*epsilIn*(1./3.)*MassIn/(4*M_PI);
      if (MassIn > 2.*muMass) muWidth = nuWidth*(1.+2.*massRatio2)*sqrt(1.-4.*massRatio2);
      // Define total width and branching ratios
      WidthIn = nuWidth + muWidth;
      nuBrRatio = nuWidth/WidthIn;
      muBrRatio = muWidth/WidthIn;

    } else if (BranchingType == 1) {
      // Change name accordingly
      IDPDG = 5500123;
      name = "DMParticleZPrimeVis";
      // Calculate widths
      nuWidth = epsilIn*epsilIn*(1./3.)*MassIn/(4*M_PI);
      if (MassIn > 2.*muMass) muWidth = nuWidth*(1.+2.*massRatio2)*sqrt(1.-4.*massRatio2);
      if (MassIn > 2.*mChi) chiWidth= MassIn/12*alphaD*pow((1-4*mChi*mChi/(MassIn*MassIn)),3./2); //Z->DM DM
      // Define total width and branching ratios
      WidthIn = nuWidth + muWidth + chiWidth;
      nuBrRatio = nuWidth/WidthIn;
      muBrRatio = muWidth/WidthIn;
      chiBrRatio = chiWidth/WidthIn;

    } else if (BranchingType == 10) {
      // Change name accordingly
      IDPDG = 5500222;
      name = "DMParticleB-LBoson";
      // Calculate widths
      if(MassIn > 600.*MeV) {G4cout << "Branching ratios for this BranchingType and this mass are not yet implemented, exiting" << G4endl; exit(1);}
      nuWidth = epsilIn*epsilIn*MassIn/(4*M_PI);
      if(MassIn > 2.*CLHEP::electron_mass_c2) eWidth = epsilIn * epsilIn / (4*M_PI) * ZPrimeWidth(CLHEP::electron_mass_c2, CLHEP::electron_mass_c2, MassIn);
      if(MassIn > 2.*muMass) muWidth = epsilIn * epsilIn / (4*M_PI) * ZPrimeWidth(muMass, muMass, MassIn);
      if(MassIn > pi0Mass) {
        hWidth = (CLHEP::fine_structure_const*epsilIn*epsilIn/(4*M_PI)*MassIn*MassIn*MassIn) /
                 (96.*M_PI*M_PI*M_PI*0.93*0.93*pi0Mass*pi0Mass);
        hWidth *= pow((1. - pi0Mass*pi0Mass/(MassIn*MassIn)),3);
        G4double a = 1. - (MassIn*MassIn)/(omegaMass*omegaMass);
        G4double b = 12.3/omegaMass;
        G4double mod2 = 1./(a*a + b*b);
        hWidth *= mod2;
      }
      // Define total width and branching ratios
      WidthIn = nuWidth + eWidth + muWidth + hWidth;
      nuBrRatio = nuWidth/WidthIn;
      eBrRatio = eWidth/WidthIn;
      muBrRatio = muWidth/WidthIn;
      hBrRatio = hWidth/WidthIn;

    } else if (BranchingType == 11) {
      // Change name accordingly
      IDPDG = 5500222;
      name = "DMParticleB-LBoson";
      // Calculate widths
      if(MassIn > 600.*MeV) {G4cout << "Branching ratios for this BranchingType and this mass are not yet implemented, exiting" << G4endl; exit(1);}
      nuWidth = epsilIn*epsilIn*MassIn/(4*M_PI);
      if(MassIn > 2.*CLHEP::electron_mass_c2) eWidth = epsilIn * epsilIn / (4*M_PI) * ZPrimeWidth(CLHEP::electron_mass_c2, CLHEP::electron_mass_c2, MassIn);
      if(MassIn > 2.*muMass) muWidth = epsilIn * epsilIn / (4*M_PI) * ZPrimeWidth(muMass, muMass, MassIn);
      if(MassIn > pi0Mass) {
        hWidth = (CLHEP::fine_structure_const*epsilIn*epsilIn/(4*M_PI)*MassIn*MassIn*MassIn) /
                 (96.*M_PI*M_PI*M_PI*0.93*0.93*pi0Mass*pi0Mass);
        hWidth *= pow((1. - pi0Mass*pi0Mass/(MassIn*MassIn)),3);
        G4double a = 1. - (MassIn*MassIn)/(782.66*782.66); // 782.66 MeV is the mass of omega meson
        G4double b = 12.3/782.66;
        G4double mod2 = 1./(a*a + b*b);
        hWidth *= mod2;
      }
      if (MassIn > 2.*mChi) chiWidth= MassIn/12*alphaD*pow((1-4*mChi*mChi/(MassIn*MassIn)),3./2); //Z->DM DM
      // Define total width and branching ratios
      WidthIn = nuWidth + eWidth + muWidth + hWidth + chiWidth;
      nuBrRatio = nuWidth/WidthIn;
      eBrRatio = eWidth/WidthIn;
      muBrRatio = muWidth/WidthIn;
      hBrRatio = hWidth/WidthIn;
      chiBrRatio = chiWidth/WidthIn;

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

      // create decay table and add modes
      G4DecayTable* table = new G4DecayTable();

      // Define decay channels according to branching type
      if (BranchingType == 0) {
        G4VDecayChannel** mode = new G4VDecayChannel*[3];
        // DMParticleZPrime -> nu_mu + anti_nu_mu
        mode[0] = new G4PhaseSpaceDecayChannel(name, nuBrRatio/2., 2, "anti_nu_mu", "nu_mu");
        // DMParticleZPrime -> nu_tau + anti_nu_tau
        mode[1] = new G4PhaseSpaceDecayChannel(name, nuBrRatio/2., 2, "anti_nu_tau", "nu_tau");
        // DMParticleZPrime -> mu+ + mu-
        mode[2] = new G4PhaseSpaceDecayChannel(name, muBrRatio, 2, "mu+", "mu-");
        for (G4int index = 0; index < 3; index++) table->Insert(mode[index]);
        delete [] mode;

      } else if (BranchingType == 1) { // Lmu-Ltau Z' boson
        G4VDecayChannel** mode = new G4VDecayChannel*[4];
        // DMParticleZPrime -> nu_mu + anti_nu_mu
        mode[0] = new G4PhaseSpaceDecayChannel(name, nuBrRatio/2., 2, "anti_nu_mu", "nu_mu");
        // DMParticleZPrime -> nu_tau + anti_nu_tau
        mode[1] = new G4PhaseSpaceDecayChannel(name, nuBrRatio/2., 2, "anti_nu_tau", "nu_tau");
        // DMParticleZPrime -> mu+ + mu-
        mode[2] = new G4PhaseSpaceDecayChannel(name, muBrRatio, 2, "mu+", "mu-");
        // DMParticleZPrime -> chi + chi
        mode[3] = new G4PhaseSpaceDecayChannel(name, chiBrRatio, 2, "DMParticleChiScalar", "DMParticleChiScalar");
        for (G4int index = 0; index < 4; index++) table->Insert(mode[index]);
        delete [] mode;

      } else if (BranchingType == 10) { // B-L Z' boson with coupling to all SM particles
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

      } else if (BranchingType == 11) { // B-L Z' boson with coupling to all SM particles + Scalar Chi
        G4VDecayChannel** mode = new G4VDecayChannel*[7];
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
        // DMParticleZPrime -> chi + chi
        mode[6] = new G4PhaseSpaceDecayChannel(name, chiBrRatio, 2, "DMParticleChiScalar", "DMParticleChiScalar");

        for (G4int index = 0; index < 7; index++) table->Insert(mode[index]);
        delete [] mode;
      }
      anInstance->SetDecayTable(table);
      anInstance->DumpTable();
    }

  }
  theInstance = reinterpret_cast<DMParticleZPrime*>(anInstance);
  return theInstance;
}
