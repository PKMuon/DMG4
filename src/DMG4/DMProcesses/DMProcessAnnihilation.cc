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
#include "G4ParticleTypes.hh"
#include "G4Positron.hh"
#include "G4ProcessManager.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4Element.hh"
#include "Randomize.hh"
#include "G4ElementVector.hh"

#include "DarkMatterParametersRegistry.hh"

#include "AnnihilationStepLimiter.hh"

#define EDEP_ALONG_STEP
#define ATOMIC_EFFECTS


DMProcessAnnihilation::DMProcessAnnihilation(DarkMatterAnnihilation *DarkMatterPointerIn, G4ParticleDefinition *theDMParticlePtrIn, G4double BiasSigmaFactorIn,AnnihilationStepLimiter *LimiterIn) :
    G4VDiscreteProcess("DMProcessAnnihilation", fUserDefined),  // fElectromagnetic
    myDarkMatterAnnihilation(DarkMatterPointerIn), theDMParticlePtr(theDMParticlePtrIn), BiasSigmaFactor(BiasSigmaFactorIn), DMpar(0), iBranchingType(0), mChi(0), mChi1(0), mChi2(
        0),m_limiter(LimiterIn) {
  SetProcessSubType(1); //fBremsstrahlung? // TODO: verify this

  DMpar = DarkMatterParametersRegistry::GetInstance();
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
  G4double deltaE = 0.005; //according to B. Banto studies (5 MeV), documentation before
  G4double Emax=(myDarkMatterAnnihilation->GetMA()*myDarkMatterAnnihilation->GetMA()-2*Mel*Mel)/(2*Mel); //this is resonant energy in GeV

  xi=Emax/(Emax+deltaE);

  G4cout<<"DMProcessAnnihilation, init xi: "<<xi<<G4endl;


  //Atomic parts
  maxShellElectronEnergy=0;


}
DMProcessAnnihilation::~DMProcessAnnihilation(){
  
}
G4bool DMProcessAnnihilation::IsApplicable(const G4ParticleDefinition &pDef) {
  return ("e+" == pDef.GetParticleName());
}



G4double DMProcessAnnihilation::GetMeanFreePath(const G4Track &aTrack, G4double, /*previousStepSize*/
G4ForceCondition* /*condition*/) {


  G4double DensityMat = aTrack.GetMaterial()->GetDensity() / (g / cm3);
  G4double etot = aTrack.GetTotalEnergy() / GeV; //this is the total energy of the positron at the beginning of the step


  //this method only checks if the material is ok, and if the energy is above threshold
  if (myDarkMatterAnnihilation->EmissionAllowed(etot, DensityMat)) {
    /*
     * This part is related to atomic motion effects
     * See: https://gitlab.cern.ch/P348/DMG4/-/issues/14?work_item_iid=21
     *
     *1- Select a random element for this material, by considering for simplicity a per-atom cross-section proportional to Z^2
     *   FOR THE MOMENT, THIS ONLY SUPPORTS ELEMENTS with one-material only.
     */

    auto elms=aTrack.GetMaterial()->GetElementVector();

    if ((*elms).size()>1){
      G4cout<<"ERROR, DMProcessAnnihilation with multi-material elements not yet supported"<<G4endl;
      exit(1);
    }
    const G4Element* elm=(*elms)[0];
    G4int Z=(G4int)(elm->GetZ());

    /*2- If necessary, compute the electron energies for this material*/
    if (shellElectronEnergies.find(Z)==shellElectronEnergies.end()){

        for (G4int is=0;is<elm->GetNbOfAtomicShells();is++){
          auto v=this->SimulateElectronEnergies(elm,is);
          shellElectronZ[Z][is]=elm->GetNbOfShellElectrons(is);
          shellElectronEnergies[Z][is]=v;
          auto eMaxLoop=*max_element(v.begin(),v.end());
          if (eMaxLoop>maxShellElectronEnergy){
            maxShellElectronEnergy=eMaxLoop;
          }
        }
    }

    /*3- Now check if the positron energy is above threshold for e+e- --> A+B (where A,B are the particles in the final state)*/
    double smin=myDarkMatterAnnihilation->sMin();
    double sMandelstam=2*Mel*Mel+2*Mel*etot; //GeV^2

#ifdef ATOMIC_EFFECTS
    /*If atomic effects are included, we should write sMandelstam=2Mel*Mel+2*Mel(E_-zP_-).
      Largest sMandelstam value: z=-1, sMandelstam=2Mel*Mel+2*E+*(E_- + P_-).
      Take the largest positron energy for the most energetic shell*/
    double EeleMAX=maxShellElectronEnergy+Mel;
    double PeleMAX=sqrt(EeleMAX*EeleMAX-Mel*Mel);
    sMandelstam=2*Mel*Mel+2*etot*(EeleMAX+PeleMAX);
#endif


    if (sMandelstam<smin) return DBL_MAX;

    /*Finally compute the cross section*/
#ifdef ATOMIC_EFFECTS
    G4double CrossSection=myDarkMatterAnnihilation->GetSigmaTotAtomicEffects(etot, shellElectronZ[Z], shellElectronEnergies[Z]);
#else
    G4double CrossSection =myDarkMatterAnnihilation->GetSigmaTot(etot); //keep this line here
#endif


    G4double Emax=(myDarkMatterAnnihilation->GetMA()*myDarkMatterAnnihilation->GetMA()-2*Mel*Mel)/(2*Mel); //this is in GeV


#ifdef EDEP_ALONG_STEP
    /*
     * This part takes into account the energy dependence of the cross section along the step.
     * See: https://gitlab.cern.ch/P348/DMG4/-/issues/14
     *
     * The xi parameter is set to Er/(Er+deltaE), where deltaE is the max energy loss across this step
     * This is obtained from the step limiter if it is available
     */
    if (m_limiter){
      G4double dEmax=m_limiter->GetMaxEloss(etot*GeV)/GeV;
      xi=Emax/(Emax+dEmax);
      if (xi<.8) //G4 suggestion
        xi=.8;
      //G4cout<<"DMProcessGetMeanFreePath: E_f="<<aTrack.GetTotalEnergy()/GeV<<"[GeV]; max energy loss: "<<dEmax<<"[GeV]; xi="<<xi<<" Emax/xi="<< Emax/xi<< "[GeV]"<< G4endl;
    }

    //First case, the energy at the beginning of the step is smaller than the resonant energy
    if (etot < Emax){
      this->CrossSectionStepE=etot;
#ifdef ATOMIC_EFFECTS
      this->CrossSectionStepVal=myDarkMatterAnnihilation->GetSigmaTotAtomicEffects(this->CrossSectionStepE, shellElectronZ[Z], shellElectronEnergies[Z]);
#else
      this->CrossSectionStepVal=myDarkMatterAnnihilation->GetSigmaTot(this->CrossSectionStepE);
#endif

    }
    //Second case, the energy at the beginning of the step is larger than the resonant energy and also than Emax/xi
    //take the value xi*E (0<xi<1)
    else if (etot > (Emax/xi)){
      this->CrossSectionStepE=xi*etot;
#ifdef ATOMIC_EFFECTS
      this->CrossSectionStepVal=myDarkMatterAnnihilation->GetSigmaTotAtomicEffects(this->CrossSectionStepE, shellElectronZ[Z], shellElectronEnergies[Z]);
#else
      this->CrossSectionStepVal=myDarkMatterAnnihilation->GetSigmaTot(this->CrossSectionStepE);
#endif
    }
    //Third case, the energy at the beginning of the step is larger than the resonant energy but smaller than Emax/xi
    else{
      this->CrossSectionStepE=Emax;
#ifdef ATOMIC_EFFECTS
      this->CrossSectionStepVal=myDarkMatterAnnihilation->GetSigmaTotAtomicEffects(this->CrossSectionStepE, shellElectronZ[Z], shellElectronEnergies[Z]);
#else
      this->CrossSectionStepVal=myDarkMatterAnnihilation->GetTotalCrossSectionMax();
#endif
    }
    //avoid numerical manipulations of this->CrossSectionStepVal
    CrossSection = this->CrossSectionStepVal;
#endif

    CrossSection *= picobarn;
    //The DarkMatterAnnihilation classes compute the cross section for eps = epsilBench. Here, we revert back to epsilon
    //We do not do this for Z' annihilation, since in that case the eps passed by user is used.
    if (myDarkMatterAnnihilation->GetDMType()!=11){
      CrossSection *= (myDarkMatterAnnihilation->Getepsil() * myDarkMatterAnnihilation->Getepsil()) / (myDarkMatterAnnihilation->GetepsilBench() * myDarkMatterAnnihilation->GetepsilBench());
    }
    CrossSection /= myDarkMatterAnnihilation->GetSigmaNorm();

    if (CrossSection<=0.){
      return DBL_MAX;
    }

    G4double n = aTrack.GetMaterial()->GetTotNbOfAtomsPerVolume(); //The annihilation cross section already contains a multiplicative factor "Z".
    G4double XMeanFreePath = 1. / (n * CrossSection);
    XMeanFreePath /= BiasSigmaFactor;

    return XMeanFreePath;

  }
  return DBL_MAX;
}

//This method is called if this process is selected as that responsible for the step limitation
G4VParticleChange* DMProcessAnnihilation::PostStepDoIt(const G4Track &aTrack, const G4Step &aStep) {
  const G4double incidentE = aTrack.GetTotalEnergy(); //this is the energy at the end of the step
  const G4double initialE =  aStep.GetPreStepPoint()->GetTotalEnergy(); //this is the energy at the beginning of the step
  G4double DMMass = myDarkMatterAnnihilation->GetMA()*GeV; // in MeV

  auto elms=aTrack.GetMaterial()->GetElementVector();
  if ((*elms).size()>1){
    G4cout<<"ERROR, DMProcessAnnihilation with multi-material elements not yet supported"<<G4endl;
    exit(1);
  }
  const G4Element* elm=GetRandomElement(elms);
  G4int Z=(G4int)(elm->GetZ());




#ifdef ATOMIC_EFFECTS
  const G4double initialCrossSection=myDarkMatterAnnihilation->GetSigmaTotAtomicEffects(initialE/GeV,shellElectronZ[Z], shellElectronEnergies[Z]);
  const G4double finalCrossSection=myDarkMatterAnnihilation->GetSigmaTotAtomicEffects(incidentE/GeV,shellElectronZ[Z], shellElectronEnergies[Z]);
#else
  const G4double initialCrossSection=myDarkMatterAnnihilation->GetSigmaTot(initialE/GeV); //this is the cross section at the beginning of the step
  const G4double finalCrossSection=myDarkMatterAnnihilation->GetSigmaTot(incidentE/GeV); //this is the cross section at the end of the step
#endif



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

  if (std::isinf(this->CrossSectionStepVal) || std::isinf(finalCrossSection)) {
    G4cout<<"ERROR, DMProcessAnnihilation with inf CS! Skipping event!"<<G4endl;
    aStep.GetTrack()->SetKineticEnergy(0.);
    aStep.GetTrack()->SetTrackStatus(fStopAndKill);
    aStep.GetPostStepPoint()->SetProcessDefinedStep(0); //important for the G4SteppingAction
    aParticleChange.Initialize(aTrack);
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }


#ifdef EDEP_ALONG_STEP
  //G4cout<<"PostStepDoIt1 "<<aStep.GetPreStepPoint()->GetTotalEnergy()/GeV<<" "<<aStep.GetPostStepPoint()->GetTotalEnergy()/GeV<<" "<<aTrack.GetTotalEnergy()/GeV<<G4endl;
  //G4cout<<"PostStepDoIt1a "<<aTrack.GetTrackID()<<" "<<aTrack.GetParentID()<<" "<<diffE<<G4endl;
  //G4cout<<"PostStepDoIt2 "<<initialCrossSection<<" "<<finalCrossSection<<G4endl;
  //G4cout<<"PostStepDoIt3 "<<this->CrossSectionStepE<<" "<<this->CrossSectionStepVal<<" "<<prob<<G4endl;

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
#ifdef ATOMIC_EFFECTS
    double maxV=myDarkMatterAnnihilation->GetTotalCrossSectionMaxAtomicEffects(shellElectronZ[Z],shellElectronEnergies[Z]);
#else
    double maxV=myDarkMatterAnnihilation->GetTotalCrossSectionMax();
#endif
    std::cout << "DM PDG ID = " << theDMParticlePtr->GetPDGEncoding() << " emitted by " << aTrack.GetDefinition()->GetParticleName() << " with energy = "
        << incidentE / GeV << " GeV, DM energy = " << incidentE / GeV << " GeV [event n.: " <<G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID()<<"]"<<std::endl;
    std::cout << "DM cross section=" <<initialCrossSection<<" [MAX VALUE: "<<maxV<<"]"<<std::endl;

    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  } else { //simulate the decay e+e- -->A' -->ff

    if (myDarkMatterAnnihilation->GetDMType() == 1) { //dark photon

      //1: Get the cosine of the final state f in the CM frame.
      G4double DMeta_CM = myDarkMatterAnnihilation->SimulateEmissionResonant(incidentE);
      //2: Get the phi angle of the final state f in the CM frame
      G4double DMphi_CM = G4UniformRand() * 2 * CLHEP::pi;

      G4LorentzVector p4_posi(aTrack.GetMomentum(),aTrack.GetTotalEnergy());
#ifdef ATOMIC_EFFECTS
      G4double zEle=2*G4UniformRand()-1;
      G4double phiEle=2*G4UniformRand()*CLHEP::pi;
      G4double kinEle=0; //TODO
      G4double momEle=0;

      G4LorentzVector p4_ele(G4ThreeVector(momEle*sqrt(1-zEle*zEle)*cos(phiEle),momEle*sqrt(1-zEle*zEle)*sin(phiEle),momEle*zEle),kinEle+CLHEP::electron_mass_c2);
#else
      G4LorentzVector p4_ele(G4ThreeVector(0,0,0),CLHEP::electron_mass_c2);
#endif
      //3: Define the Lorentz Vector of the CM (e+ + e-)
      G4LorentzVector vCM=p4_posi+p4_ele;
      G4double ss = vCM.m2();

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



const G4Element* DMProcessAnnihilation::GetRandomElement(const G4ElementVector *elms){
  G4int idx=0;
  if (elms->size()==1){
    idx=0;
  }
  else{

    G4cerr<<"ERROR! DMProcessAnnihilation::GetRandomElement called for a multi-element material, not yet implemented"<<G4endl;
    exit(1);

    /*If there'sMandelstam more than one element in the material,
      consider a probability proportional to Z^2
    */
    std::vector<G4double> cumulativeProd;
    G4double totProd=0;
    for (auto elem : (*elms)){
      totProd+=elem->GetZ()*elem->GetZ();
      cumulativeProd.push_back(totProd);
    }
    G4double xprob =  totProd*G4UniformRand();
    // select material
    for(idx=0; idx<cumulativeProd.size(); ++idx) {
      if(xprob <= cumulativeProd[idx]) {
        break;
      }
    }
  }
  return (*elms)[idx];
}

//A.C. use an exponential model, as in G4 Compton model
//An exponential model
G4double DMProcessAnnihilation::GetOneRandomEleEnergy(G4double B){
  G4double u=G4UniformRand();
  G4double x=-B*log(u);
  return x;
}

std::vector<G4double> DMProcessAnnihilation::SimulateElectronEnergies(const G4Element *elm,G4int is){

  std::vector<G4double> v;
  G4double B=elm->GetAtomicShell(is); //binding energy, positive, in g4 units

  const G4int nEne=100;
  for (G4int ii=0;ii<nEne;ii++){
    G4double x=this->GetOneRandomEleEnergy(B); //in G4 units
    v.push_back(x/GeV); //convert it in GeV
  }
  return v;
}
