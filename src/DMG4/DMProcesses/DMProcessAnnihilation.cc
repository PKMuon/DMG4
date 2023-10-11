#include "DMProcessAnnihilation.hh"

#include "DarkMatter.hh"
#include "DarkMatterAnnihilation.hh"
#include "DarkPhotons.hh"

#include "DMParticleAPrime.hh"

#include "DMParticleChi.hh"
#include "DMParticleChi1.hh"
#include "DMParticleChi2.hh"

#include "G4ProcessType.hh"
#include "G4EmProcessSubType.hh"
#include "G4SystemOfUnits.hh"

#include "DarkMatterParametersFactory.hh"

#define EDEP_ALONG_STEP

DMProcessAnnihilation::DMProcessAnnihilation(DarkMatterAnnihilation *DarkMatterPointerIn, G4ParticleDefinition *theDMParticlePtrIn, G4double BiasSigmaFactorIn) :
    G4VDiscreteProcess("DMProcessAnnihilation", fUserDefined),  // fElectromagnetic
    myDarkMatterAnnihilation(DarkMatterPointerIn), theDMParticlePtr(theDMParticlePtrIn), BiasSigmaFactor(BiasSigmaFactorIn), DMpar(0), iBranchingType(0), mChi(0), mChi1(0), mChi2(
        0) {
  SetProcessSubType(1); //fBremsstrahlung? // TODO: verify this

  DMpar = DarkMatterParametersFactory::GetInstance();
  if (DMpar) {
    iBranchingType = (int) DMpar->GetRegisteredParam("BranchingType");

    if (iBranchingType == 0) {
      mChi = DMpar->GetRegisteredParam("RDM") * myDarkMatterAnnihilation->GetMA(); //A.C. by default, DarkMatter units are GeV
    } else {
      double r = DMpar->GetRegisteredParam("RDM", 1. / 3);
      mChi1 = myDarkMatterAnnihilation->GetMA() * r;
      mChi2 = (1. + DMpar->GetRegisteredParam("Ffactor")) * mChi1;
    }
  }
  CrossSectionStepE=0;
  CrossSectionStepVal=0;

  //Init Xi
  //see documentation: https://gitlab.cern.ch/P348/DMG4/-/issues/14
  G4double deltaE = 0.005*MeV; //according to B. Banto studies, documentation before

  G4double Emax=(myDarkMatterAnnihilation->GetMA()*myDarkMatterAnnihilation->GetMA())/(2*Mel); //this is resonant energy in GeV
  if (DMpar->ExistsRegisteredParam("dEmaxPerStep")){
    deltaE=DMpar->GetRegisteredParam("dEmaxPerStep")/GeV;
  }
  xi=Emax/(Emax+deltaE);

  G4cout<<"DMProcessAnnihilation, init xi: "<<xi<<G4endl;
}

G4bool DMProcessAnnihilation::IsApplicable(const G4ParticleDefinition &pDef) {
  return ("e+" == pDef.GetParticleName());
}

G4double DMProcessAnnihilation::GetMeanFreePath(const G4Track &aTrack, G4double, /*previousStepSize*/
G4ForceCondition* /*condition*/) {
  G4double DensityMat = aTrack.GetMaterial()->GetDensity() / (g / cm3);
  G4double ekin = aTrack.GetKineticEnergy() / GeV; //this is the energy of the positron at the beginning of the step

  if (myDarkMatterAnnihilation->EmissionAllowed(ekin, DensityMat)) {
    G4double CrossSection =myDarkMatterAnnihilation->GetSigmaTot(ekin); //keep this line here

#ifdef EDEP_ALONG_STEP
    /*
     * This part takes into account the energy dependence of the cross section along the step.
     * See: https://gitlab.cern.ch/P348/DMG4/-/issues/14
     */
    G4double Emax=(myDarkMatterAnnihilation->GetMA()*myDarkMatterAnnihilation->GetMA())/(2*Mel); //this is in GeV


    //First case, the energy at the beginning of the step is smaller than the resonant energy
    if (ekin < Emax){
      this->CrossSectionStepE=ekin;
      this->CrossSectionStepVal=myDarkMatterAnnihilation->GetSigmaTot(this->CrossSectionStepE);
    }
    //Second case, the energy at the beginning of the step is larger than the resonant energy
    else if (ekin > (Emax/xi)){//take the largest value between the resonant energy and the value xi*E
      if (Emax > xi * ekin){
        this->CrossSectionStepE=Emax;
        this->CrossSectionStepVal=myDarkMatterAnnihilation->GetTotalCrossSectionMax();
      }
      else{
        this->CrossSectionStepE=xi*ekin;
        this->CrossSectionStepVal=myDarkMatterAnnihilation->GetSigmaTot(this->CrossSectionStepE);
     }
    }
    //Third case, the energy at the beginning of the step is larger than the resonant energy but smaller than Emax/xi
    else{
      this->CrossSectionStepE=Emax;
      this->CrossSectionStepVal=myDarkMatterAnnihilation->GetTotalCrossSectionMax();
    }

    //avoid numerical manipulations of this->CrossSectionStepVal
    CrossSection = this->CrossSectionStepVal;
#endif

    CrossSection *= picobarn;
    //The DarkMatterAnnihilation classes compute the cross section for eps = epsilBench. Here, we revert back to epsilon
    CrossSection *= (myDarkMatterAnnihilation->Getepsil() * myDarkMatterAnnihilation->Getepsil()) / (myDarkMatterAnnihilation->GetepsilBench() * myDarkMatterAnnihilation->GetepsilBench());
    CrossSection /= myDarkMatterAnnihilation->GetSigmaNorm();

    G4double n = aTrack.GetMaterial()->GetTotNbOfAtomsPerVolume(); //The annihilation cross section already contains a multiplicative factor "Z".
    G4double XMeanFreePath = 1. / (n * CrossSection);
    XMeanFreePath /= BiasSigmaFactor;

    return XMeanFreePath;

  }
  return DBL_MAX;
}

//This method is called if this process is selected as that responsible for the step limitation
G4VParticleChange* DMProcessAnnihilation::PostStepDoIt(const G4Track &aTrack, const G4Step &aStep) {
  const G4double incidentE = aTrack.GetKineticEnergy(); //this is the energy at the end of the step
  const G4double initialE =  aStep.GetPreStepPoint()->GetKineticEnergy(); //this is the energy at the beginning of the step

  G4double DMMass = myDarkMatterAnnihilation->GetMA()*GeV; // in MeV
  const G4double initialCrossSection=myDarkMatterAnnihilation->GetSigmaTot(initialE/GeV); //this is the cross section at the beginning of the step
  const G4double finalCrossSection=myDarkMatterAnnihilation->GetSigmaTot(incidentE/GeV); //this is the cross section at the end of the step
  G4double diffE = (initialE - incidentE)/GeV; // in GeV
  /*
   * Here are the calculations needed for the average cross-section
   *
  G4double ssin = 2 * initialE * CLHEP::electron_mass_c2; // in MeV^2
  G4double ssout = 2 * incidentE * CLHEP::electron_mass_c2; // in MeV^2
  G4double avgE = (incidentE + initialE)/2.; // in MeV
  G4double width = myDarkMatterAnnihilation->Width()*GeV; // in GeV
  G4double prefactor = myDarkMatterAnnihilation->PreFactor(avgE/GeV);
  //if (incidentE < 250.*250./2./CLHEP::electron_mass_c2) G4cout << "PASSED sigma_max! "  << std::endl;
  G4double arctanin = atan((ssin-DMMass*DMMass)/(DMMass*width));
  G4double arctanout = atan((ssout-DMMass*DMMass)/(DMMass*width));
  G4double avgCrossSection = prefactor/(diffE)/(2*CLHEP::electron_mass_c2*width*DMMass/GeV/GeV/GeV)*(arctanin-arctanout);
  const G4double prob=avgCrossSection/this->CrossSectionStepVal;
  G4cout << "PostStepDoIt: Avg.  CS = " << avgCrossSection << std::endl;
  G4cout << "PostStepDoIt: Step  CS = " << this->CrossSectionStepVal << std::endl;
  G4cout << "PostStepDoIt: Final CS = " << finalCrossSection << std::endl;
  */
  const G4double prob=finalCrossSection/this->CrossSectionStepVal;

#ifdef EDEP_ALONG_STEP
 // G4cout<<"PostStepDoIt1 "<<aStep.GetPreStepPoint()->GetKineticEnergy()/GeV<<" "<<aStep.GetPostStepPoint()->GetKineticEnergy()/GeV<<G4endl;
 // G4cout<<"PostStepDoIt1a "<<aTrack.GetTrackID()<<" "<<aTrack.GetParentID()<<G4endl;
 // G4cout<<"PostStepDoIt2 "<<initialCrossSection<<" "<<finalCrossSection<<G4endl;
 // G4cout<<"PostStepDoIt3 "<<this->CrossSectionStepE<<" "<<this->CrossSectionStepVal<<" "<<prob<<G4endl;

  G4double p=G4UniformRand();
  //Do nothing if p>prob
  if (p>prob || diffE == 0){
    aStep.GetPostStepPoint()->SetProcessDefinedStep(0); //important for the G4SteppingAction
    aParticleChange.Initialize(aTrack);
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
#endif



  //const G4double DMMass = theDMParticleAPrimePtr->GetPDGMass();
  G4ThreeVector incidentDir = aTrack.GetMomentumDirection();

  //For the e+e- --> R --> ff process, the only relevant variable is the cosine of the angle of the f in the CM frame.

  //Do not simulate the decay.
  if (myDarkMatterAnnihilation->Decay() == 0) {
    G4ThreeVector DMDirection = incidentDir;
    G4double DMM = myDarkMatterAnnihilation->GetMA() * GeV; //A.C. Dark Matter units are, by default, GeV
    G4double DMKinE = incidentE - DMM;

    G4DynamicParticle *movingDM = new G4DynamicParticle(theDMParticlePtr, DMDirection, DMKinE);
    aParticleChange.Initialize(aTrack);

    // Set DM:
    aParticleChange.SetNumberOfSecondaries(1);
    aParticleChange.AddSecondary(movingDM);
    // Kill projectile:
    aParticleChange.ProposeEnergy(0.);
    aParticleChange.ProposeTrackStatus(fStopAndKill);

    std::cout << "DM PDG ID = " << theDMParticlePtr->GetPDGEncoding() << " emitted by " << aTrack.GetDefinition()->GetParticleName() << " with energy = "
        << incidentE / GeV << " GeV, DM energy = " << incidentE / GeV << " GeV " << std::endl;
    std::cout << "DM cross section=" <<initialCrossSection<<" [MAX VALUE: "<<myDarkMatterAnnihilation->GetTotalCrossSectionMax()<<"]"<<std::endl;

    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  } else { //simulate the decay e+e- -->A' -->ff

    if (myDarkMatterAnnihilation->GetDMType() == 1) { //dark photon

      //1: Get the cosine of the final state f in the CM frame.
      G4double DMeta_CM = myDarkMatterAnnihilation->SimulateEmissionResonant(incidentE);
      //2: Get the phi angle of the final state f in the CM frame
      G4double DMphi_CM = G4UniformRand() * 2 * CLHEP::pi;
      //3: Define the Lorentz Vector of the CM (e+ + e-)
      G4LorentzVector vCM(aTrack.GetMomentum(), aTrack.GetTotalEnergy() + CLHEP::electron_mass_c2);

      G4double ss = 2 * (incidentE) * CLHEP::electron_mass_c2 + 2 * CLHEP::electron_mass_c2 * CLHEP::electron_mass_c2;
      G4double Pcm, Px_cm, Py_cm, Pz_cm;
      switch (iBranchingType) {
      case 0: //fermionic and scalar elastic
      case 1: {
        G4double Ecm = sqrt(ss) / 2;
        Pcm = sqrt(Ecm * Ecm - mChi * mChi);

        Px_cm = Pcm * sqrt(1 - DMeta_CM * DMeta_CM) * sin(DMphi_CM);
        Py_cm = Pcm * sqrt(1 - DMeta_CM * DMeta_CM) * cos(DMphi_CM);
        Pz_cm = Pcm * DMeta_CM;

        //The two lorentz vectors in the CM frame
        G4LorentzVector v1(Px_cm, Py_cm, Pz_cm, Ecm);
        G4LorentzVector v2(-Px_cm, -Py_cm, -Pz_cm, Ecm);

        v1.boost(vCM.boostVector());
        v2.boost(vCM.boostVector());

        G4DynamicParticle *movingDM1 = new G4DynamicParticle(DMParticleChi::Definition(), v1.vect());
        G4DynamicParticle *movingDM2 = new G4DynamicParticle(DMParticleChi::Definition(), v2.vect());

        aParticleChange.Initialize(aTrack);

        // Set DM:
        aParticleChange.SetNumberOfSecondaries(2);
        aParticleChange.AddSecondary(movingDM1);
        aParticleChange.AddSecondary(movingDM2);

        // Kill projectile:
        aParticleChange.ProposeEnergy(0.);
        aParticleChange.ProposeTrackStatus(fStopAndKill);

        std::cout << "DM PDG ID = " << theDMParticlePtr->GetPDGEncoding() << " emitted by " << aTrack.GetDefinition()->GetParticleName() << " with energy = "
            << incidentE / GeV << "GeV, DM energy = " << incidentE / GeV << " GeV " << std::endl;

        std::cout << "Daugther1: " << DMParticleChi::Definition()->GetPDGEncoding() << " energy= " << v1.e() / GeV << "GeV " << std::endl;
        std::cout << "Daugther2: " << DMParticleChi::Definition()->GetPDGEncoding() << " energy= " << v2.e() / GeV << "GeV " << std::endl;
        std::cout << "cos(angle) CM wrt e+: "<<DMeta_CM<<std::endl;

        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

        break;
      }
      case 2: //asymmetric DM
      {

        G4double Ecm1 = (ss + mChi1 * mChi1 - mChi2 * mChi2) / (2 * sqrt(ss));
        G4double Ecm2 = (ss - mChi1 * mChi1 + mChi2 * mChi2) / (2 * sqrt(ss));

        Pcm = sqrt(Ecm1 * Ecm1 - mChi1 * mChi1);
        Px_cm = Pcm * sqrt(1 - DMeta_CM * DMeta_CM) * sin(DMphi_CM);
        Py_cm = Pcm * sqrt(1 - DMeta_CM * DMeta_CM) * cos(DMphi_CM);
        Pz_cm = Pcm * DMeta_CM;

        //The two lorentz vectors in the CM frame
        G4LorentzVector v1(Px_cm, Py_cm, Pz_cm, Ecm1);
        G4LorentzVector v2(-Px_cm, -Py_cm, -Pz_cm, Ecm2);

        v1.boost(vCM.boostVector());
        v2.boost(vCM.boostVector());

        G4DynamicParticle *movingDM1 = new G4DynamicParticle(DMParticleChi1::Definition(), v1.vect());
        G4DynamicParticle *movingDM2 = new G4DynamicParticle(DMParticleChi2::Definition(), v2.vect());

        aParticleChange.Initialize(aTrack);

        // Set DM:
        aParticleChange.SetNumberOfSecondaries(2);
        aParticleChange.AddSecondary(movingDM1);
        aParticleChange.AddSecondary(movingDM2);

        // Kill projectile:
        aParticleChange.ProposeEnergy(0.);
        aParticleChange.ProposeTrackStatus(fStopAndKill);

        std::cout << "DM PDG ID = " << theDMParticlePtr->GetPDGEncoding() << " emitted by " << aTrack.GetDefinition()->GetParticleName() << " with energy = "
            << incidentE / GeV << "GeV, DM total energy = " << incidentE / GeV << std::endl;

        std::cout << "Daugther1: " << DMParticleChi1::Definition()->GetPDGEncoding() << " energy= " << v1.e() / GeV << "GeV " << std::endl;
        std::cout << "Daugther2: " << DMParticleChi2::Definition()->GetPDGEncoding() << " energy= " << v2.e() / GeV << "GeV " << std::endl;
        std::cout << "cos(angle) CM wrt e+: "<<DMeta_CM<<std::endl;

        return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);

        break;
      }
      }
    } //end DMtype==1
    else {
      G4cout << "DMProcessAnnihilation::PostStepDoIt DMType== " << myDarkMatterAnnihilation->GetDMType() << G4endl;
      G4cout << "The decay to final state LDM particles is not yet implemented " << G4endl;
      //revert to Decay()==0a
      G4ThreeVector DMDirection = incidentDir;
      G4double DMM = myDarkMatterAnnihilation->GetMA() * GeV;
      G4double DMKinE = incidentE - DMM;
      G4DynamicParticle *movingDM = new G4DynamicParticle(theDMParticlePtr, DMDirection, DMKinE);
      aParticleChange.Initialize(aTrack);

      // Set DM:
      aParticleChange.SetNumberOfSecondaries(1);
      aParticleChange.AddSecondary(movingDM);
      // Kill projectile:
      aParticleChange.ProposeEnergy(0.);
      aParticleChange.ProposeTrackStatus(fStopAndKill);

      std::cout << "DM PDG ID = " << theDMParticlePtr->GetPDGEncoding() << " emitted by " << aTrack.GetDefinition()->GetParticleName() << " with energy = "
          << incidentE / GeV << "GeV, DM total energy = " << incidentE / GeV << std::endl;

      return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }
  }
  return 0;
}
