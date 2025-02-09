#include "globals.hh"
#include "G4ios.hh"

#include "SteppingActionDMG4.hh"
//#include "RunAction.hh"
//#include "RunActionDMG4.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "DarkMatter.hh"
#include "DarkPhotons.hh"
#include "DarkScalars.hh"
#include "ALP.hh"

#include "DMParticles/DMParticle.hh"

#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"

#include "G4ProcessManager.hh"
#include "G4ProcessType.hh"
#include "Randomize.hh"

#include "G4ParticleTypes.hh"
#include "G4DynamicParticle.hh"
#include "G4EventManager.hh"
#include "G4TrackVector.hh"
#include "G4SystemOfUnits.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"

#include "G4ParticleTypes.hh"
#include "G4DynamicParticle.hh"
#include "G4EventManager.hh"
#include "G4TrackVector.hh"


SteppingActionDMG4::SteppingActionDMG4(DetectorConstruction* myDC, EventAction* myEA)
: myDetector(myDC), eventAction(myEA)
{
  eventAction->SetSteppingAction(this);
}


void SteppingActionDMG4::UserSteppingAction(const G4Step* aStep)
{
  if(myDetector->GetKillEvent()) {
    aStep->GetTrack()->SetKineticEnergy(0.);
    aStep->GetTrack()->SetTrackStatus(fStopAndKill);
  }

  G4StepPoint* SPointPreStep = aStep->GetPreStepPoint();
  G4StepPoint* SPointPostStep = aStep->GetPostStepPoint();

  if(SPointPostStep->GetProcessDefinedStep() != 0) {
    if(aStep->GetTrack()->GetCurrentStepNumber() == 1 && aStep->GetTrack()->GetCreatorProcess() != 0) {
      G4String creatorProcess = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
      if(creatorProcess.find(string("DMProcess")) != string::npos) {

        myDetector->SetAEmission(1);
        eventAction->CountEmission();
        eventAction->GetEventInfoPointer()->AddParentInfo(aStep);

        G4cout << "Dark Matter production with E = " << SPointPostStep->GetKineticEnergy()/GeV << " at (" << SPointPostStep->GetPosition()[0] << ", " << SPointPostStep->GetPosition()[1] << ", " << SPointPostStep->GetPosition()[2] << ")" << G4endl;
      }
    }
  }

  if(SPointPostStep->GetProcessDefinedStep() != 0) {
    if((SPointPostStep->GetProcessDefinedStep()->GetProcessName()).find(string("Decay")) != string::npos) {
      if((aStep->GetTrack()->GetParticleDefinition()->GetParticleName()).find(string("DM")) != string::npos) {

        eventAction->GetEventInfoPointer()->AddDaughterInfo(aStep);

        G4cout << "Dark Matter decay with E = " << SPointPostStep->GetKineticEnergy()/GeV << " at (" << SPointPostStep->GetPosition()[0] << ", " << SPointPostStep->GetPosition()[1] << ", " << SPointPostStep->GetPosition()[2] << ")" << G4endl;
      }
    }
  }

}


void SteppingActionDMG4::Reset()
{
  eventAction->GetDarkMatterPointer()->ResetNEmissions();
}


void SteppingActionDMG4::Finalize()
{
}
