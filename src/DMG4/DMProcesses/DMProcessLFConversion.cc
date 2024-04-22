#include "DMProcessLFConversion.hh"

#include "DarkMatter.hh"
#include "DarkLFCScalars.hh"

#include "DMParticleLFCScalar.hh"

#include "G4ProcessType.hh"
#include "G4EmProcessSubType.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4Electron.hh"
#include "G4TauMinus.hh"
#include "G4TauPlus.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"

DMProcessLFConversion::DMProcessLFConversion(DarkMatter* DarkMatterPointerIn, G4ParticleDefinition* theDMParticlePtrIn,
    G4double BiasSigmaFactorIn)
  : G4VDiscreteProcess("DMProcessLFConversion", fUserDefined),
  myDarkMatter(DarkMatterPointerIn),
  theDMParticlePtr(theDMParticlePtrIn),
  BiasSigmaFactor(BiasSigmaFactorIn)
{
  SetProcessSubType(1);
}

G4bool DMProcessLFConversion::IsApplicable(const G4ParticleDefinition& pDef)
{
  if(myDarkMatter->GetParentPDGID() == 15)
    return ("tau-" == pDef.GetParticleName() || "tau+" == pDef.GetParticleName());
  else if(myDarkMatter->GetParentPDGID() == 13)
    return ("mu-" == pDef.GetParticleName() || "mu+" == pDef.GetParticleName());
  else if(myDarkMatter->GetParentPDGID() == 11)
    return ("e-" == pDef.GetParticleName() || "e+" == pDef.GetParticleName());

  return false;
}

G4double DMProcessLFConversion::GetMeanFreePath( const G4Track& aTrack,
    G4double, /*previousStepSize*/
    G4ForceCondition* /*condition*/ )
{
  G4double DensityMat = aTrack.GetMaterial()->GetDensity()/(g/cm3);
  G4double ekin = aTrack.GetKineticEnergy()/GeV;

  if( myDarkMatter->EmissionAllowed(ekin, DensityMat) ) {

    G4double CrossSection = myDarkMatter->GetSigmaTot(ekin); //A.C. by DarkMatter definition, this is in picobarn
    CrossSection *= picobarn;

    //The DarkMatter classes compute the cross section for eps = epsilBench. Here, we revert back to epsilon
    CrossSection *= (myDarkMatter->Getepsil()* myDarkMatter->Getepsil())/(myDarkMatter->GetepsilBench()* myDarkMatter->GetepsilBench());
    CrossSection /= myDarkMatter->GetSigmaNorm();


    G4double n = aTrack.GetMaterial()->GetTotNbOfAtomsPerVolume();
    G4double XMeanFreePath = 1./(n*CrossSection);

    XMeanFreePath /= BiasSigmaFactor;


    return XMeanFreePath;

  }
  return DBL_MAX;
}

G4VParticleChange* DMProcessLFConversion::PostStepDoIt(const G4Track& aTrack,
    const G4Step & aStep)
{
  aParticleChange.Initialize(aTrack);

  // information about the incident particle
  const G4double incidentE = aTrack.GetTotalEnergy();
  G4ThreeVector incidentDir = aTrack.GetMomentumDirection();
  const G4double incidentKinE = aTrack.GetKineticEnergy();

  // compute process information (change in energy, momentum direction, ...)
  G4double XAcc=0., angles[2];

  // TODO: CHECK THIS
  //if(myDarkMatter->GetParentPDGID() == 11 || myDarkMatter->GetParentPDGID() == 13) {
    //XAcc = myDarkMatter->SimulateEmissionByMuonLFClog10(incidentE/GeV, angles); // 2-dim log-uniform sampling, angles are for the recoil lepton
  if(myDarkMatter->GetParentPDGID() == 15 || myDarkMatter->GetParentPDGID() == 13 || myDarkMatter->GetParentPDGID() == 11) {
    XAcc = myDarkMatter->SimulateEmissionVector(incidentE/GeV, angles); // 2-dim log-uniform sampling, angles are for the recoil lepton
  }

  // Check if it failed? In this case XAcc = 0
  if(XAcc > 0.) myDarkMatter->EmissionSimulated();

  G4double recoilE;
  // correct for MinitialLepton != MfinalLepton in LFC
  if(myDarkMatter->GetDMType() == 41) {
    G4double MRecoilLepton = 0.;
    // mass of recoil lepton
    if(myDarkMatter->GetDaughterPDGID() == 11) MRecoilLepton = Mel; // electron mode
    if(myDarkMatter->GetDaughterPDGID() == 13) MRecoilLepton = Mmu; // muon mode
    if(myDarkMatter->GetDaughterPDGID() == 15) MRecoilLepton = Mtau; // tau mode
    recoilE = incidentE * (1.0 - XAcc) - MRecoilLepton;           // KinE = TotE - Mass
  }
  else recoilE = incidentKinE - incidentE * XAcc;

  G4double recoilTheta = 0.,
           recoilPhi = 0.;
  G4double DMTheta = angles[0], DMPhi = angles[1];

  if(myDarkMatter->GetParentPDGID() == 11 || myDarkMatter->GetParentPDGID() == 13) {
    recoilTheta = angles[0];
    recoilPhi = angles[1];
    // TODO: modify for visible study (not for now)
    DMTheta = 0.;
    DMPhi = 0.;
  }
  G4double DME = incidentE * XAcc;
  G4double DMM = myDarkMatter->GetMA()*GeV;
  G4double DMKinE = DME - DMM;
  if(DMKinE < 0.) DMKinE = 0.;

  // Initialize DM direction vector:
  G4ThreeVector DMDirection(0., 0., .1);
  {
    DMDirection.setMag(1.);
    DMDirection.setTheta( DMTheta );
    DMDirection.setPhi( DMPhi );
    DMDirection.rotateUz(incidentDir);
  }
  // Initialize new projectile particle direction vector:
  G4ThreeVector projDirection(0., 0., 1.);
  projDirection.setMag(1.);
  projDirection.setTheta( recoilTheta );
  projDirection.setPhi( recoilPhi );
  projDirection.rotateUz(incidentDir);

  aParticleChange.SetNumberOfSecondaries(2);
  // create G4DynamicParticle object for the particle1
  G4DynamicParticle* aParticle1 = nullptr;

  // Select final state lepton accordingly
  auto charge = aTrack.GetDefinition()->GetPDGCharge();
  switch (myDarkMatter->GetDaughterPDGID()) {
    case 11:
      // Conversion to LFC scalar and Electron particle
      if(charge < 0) {
        aParticle1 = new G4DynamicParticle(G4Electron::Electron(),
            projDirection,
            recoilE);
      }
      else if(charge > 0) {
        aParticle1 = new G4DynamicParticle(G4Positron::Positron(),
            projDirection,
            recoilE);
      }
      break;
    case 13:
      // Conversion to LFC scalar and Muon particle
      if(charge < 0) {
        aParticle1 = new G4DynamicParticle(G4MuonMinus::MuonMinus(),
            projDirection,
            recoilE);
      }
      else if(charge > 0) {
        aParticle1 = new G4DynamicParticle(G4MuonPlus::MuonPlus(),
            projDirection,
            recoilE);
      }
      break;
    case 15:
      // Conversion to LFC scalar and Tau particle
      if(charge < 0) {
        aParticle1 = new G4DynamicParticle(G4TauMinus::TauMinus(),
            projDirection,
            recoilE);
      }
      else if(charge > 0) {
        aParticle1 = new G4DynamicParticle(G4TauPlus::TauPlus(),
            projDirection,
            recoilE);
      }
      break;
    default:
      break;
  }

  if (!aParticle1) {
    throw std::runtime_error("DMProcessLFConversion: ERROR: Daughter particle in LFConversion process is neither an electron, muon nor tau particle, exiting");
    exit(1);
  }

  aParticleChange.AddSecondary(aParticle1);
  // create G4DynamicParticle object for the particle2
  G4DynamicParticle* aParticle2 = new G4DynamicParticle(theDMParticlePtr,
      DMDirection,
      DMKinE);
  aParticleChange.AddSecondary(aParticle2);
  //
  // Kill the incident particle
  //
  aParticleChange.ProposeMomentumDirection(0., 0., 0.);
  aParticleChange.ProposeEnergy(0.);
  aParticleChange.ProposeTrackStatus(fStopAndKill);

  std::cout << "DM PDG ID = " << theDMParticlePtr->GetPDGEncoding() 
    << " emitted by " << aTrack.GetDefinition()->GetParticleName()
    << " with energy = " << incidentE/GeV << " DM kinetic energy = " << DMKinE/GeV << std::endl;

  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}
