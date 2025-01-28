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

#include "EventInfo.hh"


EventAction::EventAction(DetectorConstruction* myDC, DarkMatter* DMPointer)
: myDetector(myDC), myDarkMatter(DMPointer), NEmissions(0)
{
  myEventInfo = new EventInfo();
  myDetector->ResetKillEvent();
  ofstream outFile("Report.txt", std::ofstream::out);
  outFile.close();
}


EventAction::~EventAction()
{
  G4cout << "Total number of DM emissions = " << NEmissions << G4endl;
  ofstream outFile("Report.txt", std::ofstream::app);
  if(myDetector->GetAEmission()) outFile << "Total number of DM emissions " << NEmissions << G4endl;
  outFile.close();
}


void EventAction::BeginOfEventAction(const G4Event* event)
{
  theSteppingAction->Reset();
  myEventInfo->Reset();

  myDetector->SetAEmission(0);
}


void EventAction::EndOfEventAction(const G4Event* evt)
{
  theSteppingAction->Finalize();

  // Write event to file
  ofstream outFile("Report.txt", std::ofstream::app);
  outFile << "EVENT" << std::endl;
  outFile << evt->GetEventID() << " " << NEmissions << std::endl;
  myEventInfo->DumpInfo(outFile);
  outFile << "ENDEVENT" << std::endl;
  outFile.close();

  G4cout << ">>> End of event " << evt->GetEventID() << std::endl;  
}
