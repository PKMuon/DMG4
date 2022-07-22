#include "DarkMatterPhysics.hh"
#include "DarkMatterParametersFactory.hh"

#include "DarkMatter.hh"
#include "DarkPhotons.hh"
#include "DarkZ.hh"
#include "ALP.hh"
#include "DarkPhotonsAnnihilation.hh"
#include "DarkScalarsAnnihilation.hh"
#include "DarkPseudoScalarsAnnihilation.hh"
#include "DarkAxialsAnnihilation.hh"
#include "DarkScalars.hh"
#include "DarkPseudoScalars.hh"
#include "DarkAxials.hh"

#include "DMProcessDMBrem.hh"
#include "DMProcessPrimakoffALP.hh"
#include "DMProcessAnnihilation.hh"

#include "DMParticleAPrime.hh"
#include "DMParticleZPrime.hh"
#include "DMParticleALP.hh"
#include "DMParticleScalar.hh"
#include "DMParticlePseudoScalar.hh"
#include "DMParticleAxial.hh"
#include "DMParticleChi1.hh"
#include "DMParticleChi2.hh"

#include "DMParticleChi.hh"
#include "DMParticleChi1.hh"
#include "DMParticleChi2.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"

#include "G4BuilderType.hh"
#include "G4SystemOfUnits.hh"


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
  DarkMatterParametersFactory* DMpar = DarkMatterParametersFactory::GetInstance();

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
  G4double Theta = DMpar->GetRegisteredParam("Theta",0.1);
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
    case 21:
      G4cout << "Initialize ALP\n";
      myDarkMatter = new ALP(DMMass, EThresh, 1., ANucl, ZNucl, Density,  Epsilon, DecayType);
      break;
    case 31:
      G4cout << "Initialize DarkZ\n";
      myDarkMatter = new DarkZ(DMMass, EThresh, 1., ANucl, ZNucl, Density,  Epsilon, DecayType);
      break;
    case 11:
      G4cout << "Initialize DarkPhotonsAnnihilation\n";
      myDarkMatter = new DarkPhotonsAnnihilation(DMMass, EThresh, 1., ANucl, ZNucl, Density, Epsilon, DecayType,RDM,
                                                 DMpar->GetRegisteredParam("AlphaD", 0.5),BranchingType,fFactor );
      break;
    case 12:
      G4cout << "Initialize DarkScalarsAnnihilation\n";
      myDarkMatter = new DarkScalarsAnnihilation(DMMass, EThresh, 1., ANucl, ZNucl, Density, Epsilon, DecayType,RDM,
              DMpar->GetRegisteredParam("AlphaD", 0.5),BranchingType,fFactor );
      break;
    case 13:
      G4cout << "Initialize DarkAxialsAnnihilation\n";
      myDarkMatter = new DarkAxialsAnnihilation(DMMass, EThresh, 1., ANucl, ZNucl, Density, Epsilon, DecayType,RDM,
              DMpar->GetRegisteredParam("AlphaD", 0.5),BranchingType,fFactor );
      break;
    case 14:
      G4cout << "Initialize DarkPseudoScalarsAnnihilation\n";
      myDarkMatter = new DarkPseudoScalarsAnnihilation(DMMass, EThresh, 1., ANucl, ZNucl, Density, Epsilon, DecayType,RDM,
              DMpar->GetRegisteredParam("AlphaD", 0.5),BranchingType,fFactor );
       break;
     default:
       G4cout << G4endl << "Wrong DM process type specified: " << DMProcessType << " , exiting" << G4endl << G4endl;
       exit(1);
     }

   BiasSigmaFactor = DMpar->GetRegisteredParam("BiasSigmaFactor0") * 0.0001 * 0.0001 / (myDarkMatter->Getepsil()*myDarkMatter->Getepsil());
}


void DarkMatterPhysics::ConstructParticle()
{
  // This call to particle definition must be first or at least go before
  // Physics::ConstructProcess()
  DMParticleAPrime::Definition();
  DMParticleZPrime::Definition();
  DMParticleALP::Definition();
  DMParticleScalar::Definition();
  DMParticlePseudoScalar::Definition();
  DMParticleAxial::Definition();

  /*A.C.
   * The following lines are necessary to construct the particles that will be propagated for annihilation
   *
   */
  DarkMatterParametersFactory* DMpar = DarkMatterParametersFactory::GetInstance();

  G4int DMProcessType = (G4int)(DMpar->GetRegisteredParam("DMProcessType"));
  G4int DecayType = (G4int)(DMpar->GetRegisteredParam("DecayType"));
  G4int BranchingType = (G4int)(DMpar->GetRegisteredParam("BranchingType", 0.));

  switch(DMProcessType)
    {
    case 1:  //dark-photon bremmstrahlung
      if (DecayType == 0) { //Only invisible, do nothing
      }
      else {  //Require final state particles
        if (BranchingType >= 2){
          DMParticleChi1::Definition();
          DMParticleChi2::Definition();
        }
      }
      break;

    case 11: //annihilation processes
    case 12:
    case 13:
    case 14:
      if (DecayType == 0) { //Only invisible, do nothing
      }
      else {  //Require final state particles
        if ((BranchingType == 0) || (BranchingType == 1)) {
          DMParticleChi::Definition();
        } else {
          if ((BranchingType == 3)) {
            G4cout << G4endl << "Several decay channels according to BranchingType are not allowed in annihilation, exiting " << G4endl << G4endl;
            exit(1);
          }
          DMParticleChi1::Definition();
          DMParticleChi2::Definition();
        }
        break;
      }
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
  }
  if(myDarkMatter->GetParentPDGID() == 13) {
    theDMParticlePtr = DMParticleZPrime::Definition(); // Decay to SM particles
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
    phLHelper->RegisterProcess( new DMProcessDMBrem(myDarkMatter, theDMParticlePtr, BiasSigmaFactor),
                                G4Electron::ElectronDefinition() );
    phLHelper->RegisterProcess( new DMProcessDMBrem(myDarkMatter, theDMParticlePtr, BiasSigmaFactor),
                                G4Positron::PositronDefinition() );
  }
  if(myDarkMatter->GetParentPDGID() == -11) {
    phLHelper->RegisterProcess( new DMProcessAnnihilation(myDarkMatter, theDMParticlePtr, BiasSigmaFactor),
                                G4Positron::PositronDefinition() );
  }
  if(myDarkMatter->GetParentPDGID() == 13) {
    phLHelper->RegisterProcess( new DMProcessDMBrem(myDarkMatter, theDMParticlePtr, BiasSigmaFactor),
                                G4MuonMinus::MuonMinusDefinition() );
    phLHelper->RegisterProcess( new DMProcessDMBrem(myDarkMatter, theDMParticlePtr, BiasSigmaFactor),
                                G4MuonPlus::MuonPlusDefinition() );
  }
  if(myDarkMatter->GetParentPDGID() == 22) {
    phLHelper->RegisterProcess( new DMProcessPrimakoffALP(myDarkMatter, theDMParticlePtr, BiasSigmaFactor),
                                G4Gamma::GammaDefinition() );
  }
}
