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
: eventAction(myEA)
{
  eventAction->SetSteppingAction(this);
}


void SteppingActionDMG4::UserSteppingAction(const G4Step* aStep)
{
  G4StepPoint* SPointPreStep = aStep->GetPreStepPoint();
  //G4StepPoint* SPointPostStep = aStep->GetPostStepPoint();
  //G4double ekin = aStep->GetTrack()->GetKineticEnergy()/GeV;
  //G4double ekinprestep = SPointPreStep->GetKineticEnergy()/GeV;
  //G4ParticleDefinition* theParticleDefinition = aStep->GetTrack()->GetDefinition();

  if(aStep->GetPostStepPoint()->GetProcessDefinedStep() != 0) {
    if((aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()).find(string("DMProcess")) != string::npos) {

      eventAction->CountEmission();

      G4cout << "Dark Matter production at E = " << SPointPreStep->GetKineticEnergy()/GeV << G4endl;
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
