#pragma once

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class DetectorConstruction;
class EventAction;


class SteppingActionDMG4 : public G4UserSteppingAction
{
  public:
    SteppingActionDMG4(DetectorConstruction* myDC, EventAction* myEA);
    virtual ~SteppingActionDMG4(){};

    virtual void UserSteppingAction(const G4Step*);
    
    virtual void Reset();

    virtual void Finalize();

  private:

    EventAction* eventAction;
};
