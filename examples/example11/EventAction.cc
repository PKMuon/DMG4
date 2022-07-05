#include "globals.hh"

#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "SteppingActionDMG4.hh"

#include "DarkMatter.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4UImanager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"


EventAction::EventAction(DetectorConstruction* myDC, DarkMatter* DMPointer)
: myDetector(myDC), myDarkMatter(DMPointer)
{;}


EventAction::~EventAction()
{}


void EventAction::BeginOfEventAction(const G4Event* event)
{
  theSteppingAction->Reset();

  myDetector->SetAEmission(0);
}


void EventAction::EndOfEventAction(const G4Event* evt)
{
  theSteppingAction->Finalize();

  G4cout << ">>> End of event " << evt->GetEventID() << endl;  
}
