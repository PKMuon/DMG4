#pragma once

#include "G4UserStackingAction.hh"
#include "G4ClassificationOfNewTrack.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction;
class EventAction;

class StackingAction : public G4UserStackingAction
{
  public:
    StackingAction(DetectorConstruction* myDC);
    virtual ~StackingAction();

    virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*);
    virtual void NewStage();
    virtual void PrepareNewEvent();

  private:

    DetectorConstruction* myDetector;
    G4int stage;



};
