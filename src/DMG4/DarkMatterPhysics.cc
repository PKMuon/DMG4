#include "DarkMatterPhysics.hh"
#include "DarkMatterParametersRegistry.hh"

#include "DarkMatter.hh"
#include "DarkMatterAnnihilation.hh"
#include "DarkPhotons.hh"
#include "DarkZ.hh"
#include "DarkZAnnihilation.hh"
#include "DarkMuPhilicScalars.hh"
#include "DarkMuPhilicPseudoScalars.hh"
#include "ALP.hh"
#include "DarkPhotonsAnnihilation.hh"
#include "DarkScalarsAnnihilation.hh"
#include "DarkPseudoScalarsAnnihilation.hh"
#include "DarkAxialsAnnihilation.hh"
#include "DarkScalars.hh"
#include "DarkPseudoScalars.hh"
#include "DarkMassSpin2.hh"
#include "DarkMassSpin2Annihilation.hh"
#include "DarkAxials.hh"

#include "DMProcessDMBrem.hh"
#include "DMProcessPrimakoffALP.hh"
#include "DMProcessAnnihilation.hh"
#include "AnnihilationStepLimiter.hh"

#include "DMParticleAPrime.hh"
#include "DMParticleZPrime.hh"
#include "DMParticleALP.hh"
#include "DMParticleScalar.hh"
#include "DMParticlePseudoScalar.hh"
#include "DMParticleAxial.hh"

#include "DMParticleChi.hh"
#include "DMParticleChiScalar.hh"
#include "DMParticleChi1.hh"
#include "DMParticleChi2.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"

#include "G4BuilderType.hh"
#include "G4SystemOfUnits.hh"
#include "G4ProcessManager.hh"


DarkMatterPhysics::DarkMatterPhysics() 
: G4VPhysicsConstructor("DarkMatterPhysics")
{
  SetPhysicsType(bUnknown);
  //fMessenger = new DarkMatterPhysicsMessenger();

  if(!DarkMatterPhysicsConfigure()) {
    G4cout << "Dark Matter physics is not properly configured, exiting" << G4endl;
    exit(1);
  }
  this->Init();
 }


DarkMatterPhysics::DarkMatterPhysics(void *ptr)
: G4VPhysicsConstructor("DarkMatterPhysics")
{
  SetPhysicsType(bUnknown);
  //fMessenger = new DarkMatterPhysicsMessenger();

  if(!DarkMatterPhysicsConfigure(ptr)) {
    G4cout << "Dark Matter physics is not properly configured, exiting" << G4endl;
    exit(1);
  }
  this->Init();
}


DarkMatterPhysics::~DarkMatterPhysics()
{
  if(myDarkMatter) delete myDarkMatter;
}

void DarkMatterPhysics::Init(){
  //call an instance of the class
  DarkMatterParametersRegistry* DMpar = DarkMatterParametersRegistry::GetInstance();

  G4double EThresh = DMpar->GetRegisteredParam("EThresh");
  G4int DMProcessType = DMpar->GetRegisteredParam("DMProcessType");
  double DMMass   = DMpar->GetRegisteredParam("DMMass");
  double Epsilon  = DMpar->GetRegisteredParam("Epsilon");
  G4double ANucl     = DMpar->GetRegisteredParam("ANucl");
  G4double ZNucl     = DMpar->GetRegisteredParam("ZNucl");
  G4double Density   = DMpar->GetRegisteredParam("Density");
  G4int DecayType = DMpar->GetRegisteredParam("DecayType");

  G4double RDM =  DMpar->GetRegisteredParam("RDM",1./3);
  G4double fFactor = DMpar->GetRegisteredParam("Ffactor",0.1);
  //G4double IDMTheta = DMpar->GetRegisteredParam("IDMTheta",1.e-3); // is it to be used somewhere as argument?
  G4int BranchingType = DMpar->GetRegisteredParam("BranchingType",0);

/*
 * A.C. all quantities obtained from DMPar have intrinsic G4 units
 * In the following, we pass them to DarkMatter classes, that use following convention:
 *
 * Masses, energies => GeV
 * Density => g/cm3
 *
 * We convert them here
 */

  EThresh/=GeV;
  DMMass/=GeV;
  Density/=(g/cm3);



  switch(DMProcessType)
    {
    case 1:
      G4cout << "Initialize DarkPhotons\n";
      myDarkMatter = new DarkPhotons(DMMass, EThresh, 1., ANucl, ZNucl, Density, Epsilon, DecayType);
      break;
    case 2:
      G4cout << "Initialize DarkScalars\n";
      myDarkMatter = new DarkScalars(DMMass, EThresh, 1., ANucl, ZNucl, Density, Epsilon, DecayType);
      break;
    case 3:
      G4cout << "Initialize DarkAxials\n";
      myDarkMatter = new DarkAxials(DMMass, EThresh, 1., ANucl, ZNucl, Density,  Epsilon, DecayType);
      break;
    case 4:
      G4cout << "Initialize DarkPseudoScalars\n";
      myDarkMatter = new DarkPseudoScalars(DMMass, EThresh, 1., ANucl, ZNucl, Density,  Epsilon, DecayType);
      break;
    case 5:
      G4cout << "Initialize DarkMassSpin2\n";
      if(DecayType) { // Temporary plug
        G4cout << G4endl << "DarkMassSpin2 with decays is not yet implemented, exiting" << G4endl << G4endl;
        exit(1);
      }
      myDarkMatter = new DarkMassSpin2(DMMass, EThresh, 1., ANucl, ZNucl, Density,  Epsilon, DecayType);
      break;
    case 21:
      G4cout << "Initialize ALP\n";
      myDarkMatter = new ALP(DMMass, EThresh, 1., ANucl, ZNucl, Density,  Epsilon, DecayType);
      break;
    case 31:
      G4cout << "Initialize DarkZ\n";
      myDarkMatter = new DarkZ(DMMass, EThresh, 1., ANucl, ZNucl, Density,  Epsilon, DecayType);
      break;
    case 32:
      G4cout << "Initialize DarkMuPhilicScalars\n";
      if(DecayType) { // Temporary plug
        G4cout << G4endl << "DarkMuPhilicScalar with decays is not yet implemented, exiting" << G4endl << G4endl;
        exit(1);
      }
      myDarkMatter = new DarkMuPhilicScalars(DMMass, EThresh, 1., ANucl, ZNucl, Density,  Epsilon, DecayType);
      break;
    case 34:
      G4cout << "Initialize DarkMuPhilicPseudoScalars\n";
      if(DecayType) { // Temporary plug
        G4cout << G4endl << "DarkMuPhilicPseudoScalar with decays is not yet implemented, exiting" << G4endl << G4endl;
        exit(1);
      }
      myDarkMatter = new DarkMuPhilicPseudoScalars(DMMass, EThresh, 1., ANucl, ZNucl, Density,  Epsilon, DecayType);
      break;
    case 11:
      G4cout << "Initialize DarkPhotonsAnnihilation\n";
      myDarkMatter = new DarkPhotonsAnnihilation(DMMass, EThresh, 1., ANucl, ZNucl, Density, Epsilon, DecayType, RDM,
                                                 DMpar->GetRegisteredParam("AlphaD", 0.5), BranchingType, fFactor);
      break;
    case 12:
      G4cout << "Initialize DarkScalarsAnnihilation\n";
      myDarkMatter = new DarkScalarsAnnihilation(DMMass, EThresh, 1., ANucl, ZNucl, Density, Epsilon, DecayType, RDM,
              DMpar->GetRegisteredParam("AlphaD", 0.5), BranchingType, fFactor );
      break;
    case 13:
      G4cout << "Initialize DarkAxialsAnnihilation\n";
      myDarkMatter = new DarkAxialsAnnihilation(DMMass, EThresh, 1., ANucl, ZNucl, Density, Epsilon, DecayType, RDM,
              DMpar->GetRegisteredParam("AlphaD", 0.5), BranchingType, fFactor );
      break;
    case 14:
      G4cout << "Initialize DarkPseudoScalarsAnnihilation\n";
      myDarkMatter = new DarkPseudoScalarsAnnihilation(DMMass, EThresh, 1., ANucl, ZNucl, Density, Epsilon, DecayType, RDM,
              DMpar->GetRegisteredParam("AlphaD", 0.5), BranchingType, fFactor );
       break;
    case 15:
      G4cout << "Initialize DarkMassSpin2Annihilation\n";
      myDarkMatter = new DarkMassSpin2Annihilation(DMMass, EThresh, 1., ANucl, ZNucl, Density, Epsilon, DecayType, RDM,
              DMpar->GetRegisteredParam("AlphaD", 0.5), BranchingType, fFactor );
       break;
    case 16:
      G4cout << "Initialize DarkZAnnihilation\n";
      myDarkMatter = new DarkZAnnihilation(DMMass, EThresh, 1., ANucl, ZNucl, Density, Epsilon, DecayType, RDM,
          DMpar->GetRegisteredParam("AlphaD", 0.5), BranchingType,fFactor);
      break;
     default:
       G4cout << G4endl << "Wrong DM process type specified: " << DMProcessType << " , exiting" << G4endl << G4endl;
       exit(1);
     }

   BiasSigmaFactor = DMpar->GetRegisteredParam("BiasSigmaFactor0") * (myDarkMatter->GetepsilBench()*myDarkMatter->GetepsilBench()) / (myDarkMatter->Getepsil()*myDarkMatter->Getepsil());

   //For the e+ e- --> Z' --> ff process, we compute the cross section using epsil, so the code above has to be changed
   if (DMProcessType==16){
     BiasSigmaFactor = DMpar->GetRegisteredParam("BiasSigmaFactor0");
   }
}


void DarkMatterPhysics::ConstructParticle()
{
  /*A.C.
   * The following lines are necessary to construct the particles that will be propagated for annihilation
   *
   */
  DarkMatterParametersRegistry* DMpar = DarkMatterParametersRegistry::GetInstance();

  G4int DMProcessType = (G4int)(DMpar->GetRegisteredParam("DMProcessType"));
  G4int DecayType = (G4int)(DMpar->GetRegisteredParam("DecayType"));
  G4int BranchingType = (G4int)(DMpar->GetRegisteredParam("BranchingType", 0.));

  switch(DMProcessType)
    {
    case 1:  //dark-photon bremmstrahlung
      DMParticleAPrime::Definition();
      if (DecayType == 0) { //Only invisible, do nothing
      }
      else {  //Require final state particles
        if (BranchingType >= 2){
          DMParticleChi1::Definition();
          DMParticleChi2::Definition();
        }
      }
      break;
    case 2:
      DMParticleScalar::Definition();
      break;
    case 3:
      DMParticleAxial::Definition();
      break;
    case 4:
      DMParticlePseudoScalar::Definition();
      break;
    case 5:
      DMParticleAPrime::Definition();
      break;
    case 11: //annihilation processes
      DMParticleAPrime::Definition();
      break;
    case 12:
      DMParticleScalar::Definition();
      break;
    case 13:
      DMParticleAxial::Definition();
      break;
    case 14:
      DMParticlePseudoScalar::Definition();
      break;
    case 15:
      DMParticleAPrime::Definition(); // A' for the moment, the spin 2 particle not yet implemented
      if (DecayType == 0) { //Only invisible, do nothing
        break;
      }
      else {  //Require final state particles
        if ((BranchingType == 0) || (BranchingType == 1)) {
          DMParticleChi::Definition();
        } else {
          if (BranchingType == 3) {
            G4cout << G4endl << "Several decay channels according to BranchingType are not allowed in annihilation, exiting " << G4endl << G4endl;
            exit(1);
          }
          DMParticleChi1::Definition();
          DMParticleChi2::Definition();
        }
        break;
      }
    case 16:
      DMParticleZPrime::Definition();
      if (DecayType == 0) { //Only invisible, do nothing
      }
      else {  //Require final state particles
        if ((BranchingType == 0) || (BranchingType == 10)) { // neutrinos final state
        } else if ((BranchingType == 1) || (BranchingType == 11)) { // DM final state
          DMParticleChiScalar::Definition();
        } else {
          G4cout << G4endl << "BranchingType not implemented, exiting " << G4endl << G4endl;
          exit(1);
        }
      }
      break;
    case 21:
      DMParticleALP::Definition();
      break;
    case 31:
      DMParticleZPrime::Definition();
      break;
    case 32:
      DMParticleZPrime::Definition();
      break;
    case 34:
      DMParticleZPrime::Definition();
      break;
    default:
      break;
    }
}


void DarkMatterPhysics::ConstructProcess()
{
  // Which DM particle?
  G4ParticleDefinition* theDMParticlePtr = 0;
  if(myDarkMatter->GetParentPDGID() == 11) {
    if(myDarkMatter->GetDMType() == 1) {
      theDMParticlePtr = DMParticleAPrime::Definition();
    }
    if(myDarkMatter->GetDMType() == 2) {
      theDMParticlePtr = DMParticleScalar::Definition();
    }
    if(myDarkMatter->GetDMType() == 3) {
      theDMParticlePtr = DMParticleAxial::Definition();
    }
    if(myDarkMatter->GetDMType() == 4) {
      theDMParticlePtr = DMParticlePseudoScalar::Definition();
    }
    if(myDarkMatter->GetDMType() == 5) {
      theDMParticlePtr = DMParticleAPrime::Definition(); // A' for the moment, the spin 2 particle not yet implemented
    }
  }
  if(myDarkMatter->GetParentPDGID() == -11) { // Annihilation
    if(myDarkMatter->GetDMType() == 1) {
      theDMParticlePtr = DMParticleAPrime::Definition();
    }
    if(myDarkMatter->GetDMType() == 2) {
      theDMParticlePtr = DMParticleScalar::Definition();
    }
    if(myDarkMatter->GetDMType() == 3) {
      theDMParticlePtr = DMParticleAxial::Definition();
    }
    if(myDarkMatter->GetDMType() == 4) {
      theDMParticlePtr = DMParticlePseudoScalar::Definition();
    }
    if(myDarkMatter->GetDMType() == 5) {                    // Annihilation through spin 2 DM, A' for the moment
      theDMParticlePtr = DMParticleAPrime::Definition();
    }
    if(myDarkMatter->GetDMType() == 11) {                    // Annihilation through Z' (Lmu-Ltau or B-L models)
      theDMParticlePtr = DMParticleZPrime::Definition();
    }
  }
  if(myDarkMatter->GetParentPDGID() == 13) {
    theDMParticlePtr = DMParticleZPrime::Definition(); // Always Z' for the moment, scalar etc. particles from muons not yet implemented
  }
  if(myDarkMatter->GetParentPDGID() == 22) {
    theDMParticlePtr = DMParticleALP::Definition();
  }

  if(!theDMParticlePtr) {G4cout << "DarkMatterPhysics::ConstructProcess: did not manage to determine the DM particle type, exiting" << G4endl; exit(1);}

  myDarkMatter->SetMA(theDMParticlePtr->GetPDGMass()/GeV);
  myDarkMatter->SetDMPDGID(theDMParticlePtr->GetPDGEncoding());
  myDarkMatter->PrepareTable();

  G4PhysicsListHelper * phLHelper = G4PhysicsListHelper::GetPhysicsListHelper();

  phLHelper->DumpOrdingParameterTable();

  // if one need to (re-)associate certain process with a particle, note
  // the following snippet
  //G4ProcessManager * pMgr = Mocktron::Definition()->GetProcessManager();
  //pmanager->RemoveProcess(idxt);
  //pmanager->AddProcess(new G4MonopoleTransportation(fMpl),-1, 0, 0);

  // ... here one can set up the model parameters from external config
  //     sources, internal attributes previously set by messengers, etc

  // ... here the processes asociated with new physics should be registered
  //     as follows
  
  if(myDarkMatter->GetParentPDGID() == 11) {
    DMProcessDMBrem* DMBremPointer = new DMProcessDMBrem(myDarkMatter, theDMParticlePtr, BiasSigmaFactor);

    // Instead of using ordtable:
    G4ProcessManager* processManager = (G4Electron::ElectronDefinition())->GetProcessManager();
    processManager->AddDiscreteProcess(DMBremPointer);
    processManager = (G4Positron::PositronDefinition())->GetProcessManager();
    processManager->AddDiscreteProcess(DMBremPointer);

//    phLHelper->RegisterProcess( DMBremPointer, G4Electron::ElectronDefinition() );
//    phLHelper->RegisterProcess( DMBremPointer, G4Positron::PositronDefinition() );
  }
  if(myDarkMatter->GetParentPDGID() == -11) {
    DarkMatterParametersRegistry* DMpar = DarkMatterParametersRegistry::GetInstance();
    G4double annihilationStepLimiterFactor=DMpar->GetRegisteredParam("AnnihilationStepLimiterFactor",5.);
    DarkMatterAnnihilation *dmAnnihil=dynamic_cast<DarkMatterAnnihilation*>(myDarkMatter);
    AnnihilationStepLimiter *dmLimiterProc=new AnnihilationStepLimiter(dmAnnihil,"StepLimiterAnnihilation");
    dmLimiterProc->SetFactor(annihilationStepLimiterFactor);

    DMProcessAnnihilation* dmAnnihilProc = new DMProcessAnnihilation(dmAnnihil, theDMParticlePtr, BiasSigmaFactor, dmLimiterProc);

    /* Add the step limiter to the list of discrete processes for the e+
     * before adding the annihilation process
     * This will trigger the calculation of the max step length, including the maximum energy loss across the new step
     * So these quantities can be used by this process for this step
     */
    G4ParticleDefinition* posi=G4Positron::Definition();
    G4ProcessManager* processManager=posi->GetProcessManager();
    processManager->AddDiscreteProcess(dmLimiterProc);

    //Now add the annihilation process
    //phLHelper->RegisterProcess(dmAnnihilProc, G4Positron::PositronDefinition()); // old way, required ordtable
    processManager->AddDiscreteProcess(dmAnnihilProc);

  }
  if(myDarkMatter->GetParentPDGID() == 13) {
    DMProcessDMBrem* DMBremPointer = new DMProcessDMBrem(myDarkMatter, theDMParticlePtr, BiasSigmaFactor);

    // Instead of using ordtable:
    G4ProcessManager* processManager = (G4MuonMinus::MuonMinusDefinition())->GetProcessManager();
    processManager->AddDiscreteProcess(DMBremPointer);
    processManager = (G4MuonPlus::MuonPlusDefinition())->GetProcessManager();
    processManager->AddDiscreteProcess(DMBremPointer);

    //phLHelper->RegisterProcess( DMBremPointer, G4MuonMinus::MuonMinusDefinition() );
    //phLHelper->RegisterProcess( DMBremPointer, G4MuonPlus::MuonPlusDefinition() );
  }
  if(myDarkMatter->GetParentPDGID() == 22) {
    DMProcessPrimakoffALP* DMPrimakoffALPPointer = new DMProcessPrimakoffALP(myDarkMatter, theDMParticlePtr, BiasSigmaFactor);

    // Instead of using ordtable:
    G4ProcessManager* processManager = (G4Gamma::GammaDefinition())->GetProcessManager();
    processManager->AddDiscreteProcess(DMPrimakoffALPPointer);

    //phLHelper->RegisterProcess( DMPrimakoffALPPointer, G4Gamma::GammaDefinition() );
  }
}
