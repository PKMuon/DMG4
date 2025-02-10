#include "globals.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "SteppingActionDMG4.hh"
#include "StackingAction.hh"

#include "QGSP_BERT.hh"
#include "FTFP_BERT.hh"

#include "DarkMatterPhysics.hh"
#include "DarkMatter.hh"

#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#include "G4RunManager.hh"
#include "G4PhysListFactory.hh"
#include "G4ios.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "Configure.hh"
#include <iostream>


int main(int argc,char** argv) {

  char * wCardFile = nullptr;
  int parserOutFlag = parseArguments(argc, argv, wCardFile);
  // Exit but before print usage
  if(parserOutFlag < 0 ) {
    print_usage( std::cerr, argv[0] );
    G4cerr << "Exit due to previous error(s)." << std::endl;
  }
  // Exit with failure
  if (parserOutFlag < 1) return EXIT_FAILURE;

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes
  DetectorConstruction* mkexp = new DetectorConstruction;
  runManager->SetUserInitialization(mkexp);

  // ___ Here the "extension" part starts ___
  G4PhysListFactory factory;
  G4VModularPhysicsList * phys = factory.GetReferencePhysList("FTFP_BERT");
  // ^^^ most of the standard physics lists are available by this interface

//  G4PhysicsListHelper * phLHelper = G4PhysicsListHelper::GetPhysicsListHelper();
//  phLHelper->DumpOrdingParameterTable();

  DarkMatterPhysics* myPhysics = new DarkMatterPhysics();
  phys->RegisterPhysics(myPhysics);
  // ^^^ Here the "extension" part ends ^^^

  runManager->SetUserInitialization(phys);  // init phys


#ifdef G4VIS_USE
  // Visualization, if you choose to have it!
  G4VisManager* visManager = new G4VisExecutive;

  visManager->Initialize();
#endif
   
  // UserAction classes
  runManager->SetUserAction(new PrimaryGeneratorAction(mkexp));

  // Define UserEventAction to store and write event information to output
  EventAction* myEA = new EventAction(mkexp, myPhysics->GetDarkMatterPointer());
  runManager->SetUserAction(myEA);
  // Define UserSteppingAction to save DM information
  runManager->SetUserAction(new SteppingActionDMG4(mkexp, myEA));
  // Define UserStackingAction to skip non-DM events
  runManager->SetUserAction(new StackingAction(mkexp));

  // User interactions
  G4UImanager * UI = G4UImanager::GetUIpointer();  

  if(argc==1)
  // Define (G)UI terminal for interactive mode  
  { 
    // G4UIterminal is a (dumb) terminal.
    G4UIsession * session = new G4UIterminal;
    UI->ApplyCommand("/control/execute prerun.g4mac");    
    session->SessionStart();
    delete session;
  }
  else
  // Batch mode
  { 
    G4String command = "/control/execute ";
    G4String fileName = wCardFile;
    //G4String fileName = argv[1];
    UI->ApplyCommand(command+fileName);
  }

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}
