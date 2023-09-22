#include "DMProcessLFConversion.hh"

#include "DarkMatter.hh"
#include "DarkLFCScalars.hh"

#include "DMParticleLFCScalar.hh"

#include "G4ProcessType.hh"
#include "G4EmProcessSubType.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"
#include "G4TauMinus.hh"
#include "G4TauPlus.hh"
#include "G4MuonMinus.hh"

DMProcessLFConversion::DMProcessLFConversion(DarkMatter* DarkMatterPointerIn, G4ParticleDefinition* theDMParticlePtrIn,
    G4double BiasSigmaFactorIn)
  : G4VDiscreteProcess("DMProcessLFConversion", fUserDefined),
  myDarkMatter(DarkMatterPointerIn),
  theDMParticlePtr(theDMParticlePtrIn),
  BiasSigmaFactor(BiasSigmaFactorIn)
{
  SetProcessSubType(2);
}

G4bool DMProcessLFConversion::IsApplicable(const G4ParticleDefinition& pDef)
{
  return ( "tau-" == pDef.GetParticleName() || "mu-" == pDef.GetParticleName() );
}

G4double DMProcessLFConversion::GetMeanFreePath( const G4Track& aTrack,
    G4double, /*previousStepSize*/
    G4ForceCondition* /*condition*/ )
{
  G4double DensityMat = aTrack.GetMaterial()->GetDensity()/(g/cm3);
  G4double ekin = aTrack.GetKineticEnergy()/GeV;
  if( myDarkMatter->EmissionAllowed(ekin, DensityMat) ) {

    G4double XMeanFreePath = myDarkMatter->GetMeanFreePathFactor()/myDarkMatter->GetSigmaTot(ekin);
    XMeanFreePath /= BiasSigmaFactor;

    //std::cout << "DMMeanFreePath = " << XMeanFreePath << std::endl;

    return XMeanFreePath;

  }
  return DBL_MAX;
}

G4VParticleChange* DMProcessLFConversion::PostStepDoIt(const G4Track& aTrack,
    const G4Step & aStep)
{
  aParticleChange.Initialize(aTrack);
  G4Material* aMaterial = aTrack.GetMaterial();

  // information about the incident particle
  const G4double incidentE = aTrack.GetTotalEnergy();
  G4ThreeVector incidentDir = aTrack.GetMomentumDirection();
  const G4double incidentKinE = aTrack.GetKineticEnergy();

  // compute process information (change in energy, momentum direction, ...)
  G4double XAcc=0., angles[2];

  if(myDarkMatter->GetParentPDGID() == 15 || myDarkMatter->GetParentPDGID() == 13) {
    // TODO: apply correct sampling
    XAcc = myDarkMatter->SimulateEmissionByMuonLFClog10(incidentE/GeV, angles); // 2-dim log-uniform sampling, angles are for the recoil lepton

  }

  // Check if it failed? In this case XAcc = 0

  if(XAcc > 0.) myDarkMatter->EmissionSimulated();

  G4double recoilE;
  // correct for MinitialLepton != MfinalLepton in LFC
  if(myDarkMatter->GetDMType() == 99) {
    G4double MRecoilLepton = 0.;
    // mass of recoil lepton
    if(myDarkMatter->GetParentPDGID() == 13) MRecoilLepton = Mtau; // muon mode
    if(myDarkMatter->GetParentPDGID() == 15) MRecoilLepton = Mmu; // tau mode
    recoilE = incidentE * (1.0 - XAcc) - MRecoilLepton;           // KinE = TotE - Mass
  }
  else recoilE = incidentKinE - incidentE * XAcc;

  G4double recoilTheta = 0.,
           recoilPhi = 0.;
  G4double DMTheta = angles[0], DMPhi = angles[1];
  if(myDarkMatter->GetParentPDGID() == 15 || myDarkMatter->GetParentPDGID() == 13) {
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
  G4DynamicParticle* aParticle1;
  // muon mode
  if(myDarkMatter->GetParentPDGID() == 13) {
    aParticle1 = new G4DynamicParticle(G4TauMinus::TauMinus(), 
        projDirection,
        recoilE);
  } 
  // tau mode
  if(myDarkMatter->GetParentPDGID() == 15) {
    aParticle1 = new G4DynamicParticle(G4MuonMinus::MuonMinus(), 
        projDirection,
        recoilE);
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
