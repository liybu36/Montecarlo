//---------------------------------------------------------------------------//
//                                                                           //
//                                                                           //
//                         G4DS Simulation                                   //
//                                                                           //
// --------------------------------------------------------------------------//

#include "DSLogger.hh"
#include "DSManager.hh"
#include "DSStorage.hh"
#include "DSRunAction.hh"
#include "DSDetectorConstruction.hh"
#include "DSPrimaryGeneratorAction.hh"
#include "DSPhysicsList.hh"
#include "DSEventAction.hh"
#include "DSTrackingAction.hh"
#include "DSSteppingAction.hh"
#include "DSStackingAction.hh"


#include "DSIO.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4UIsession.hh"

#ifdef G4VIS_USE
#include "DSVisManager.hh"
#include "G4TrajectoryDrawByParticleID.hh"
#endif

#ifdef G4UI_USE_ROOT
#include "G4UIRoot.hh"
#endif


#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#ifdef G4LIB_USE_GDML
#include "G4TransportationManager.hh"
#include "G4GDMLParser.hh"
#endif

#include <iostream>
#include <fstream>
#include <string>
#include "G4ios.hh"
#include <stdlib.h>



//Global variable to count the amount of photons
DSManager*                      runManager;
DSDetectorConstruction*         theDetector;
DSPhysicsList*                  thePhysicsList;
// Functions called from main().
void PrintHeader();
void PrintUsage(void);

using namespace std;


int main(int argc,char** argv) {
  
  PrintHeader();
  G4bool IsStack = false ;
  
  //---------------------------------------------------------------------------//
  // Definition of promary parameters 
  DSLogger::SetSeverity(DSLogger::routine);
  DSIO::Get()->CheckFileName("output");
  string filename ;
  if (argc>1) {
    ifstream ifs(argv[1]);
    string s ;
    while (getline (ifs, s)) {
       if (!s.find("/DSlog")) {
          DSLogger::SetSeverity(DSLogger::toEnum(s.substr(s.find(" ")+1)));
      } else if (!s.find("/run/filename")) {
          filename = s.substr(s.find(" ")+1);
	  //DSIO::Get()->CheckFileName(s.substr(s.find(" ")+1));
      }  else if (!s.find("/ds/stack")) {
          IsStack = true ;
      }
    } ifs.close();
  }
  if(argc==3)  filename = argv[2];
  DSIO::Get()->CheckFileName(filename);
  
  
  DSLog(routine) << "Output file name: " << DSIO::Get()->GetFileName() <<   endlog;
  // Copy fo the input file in the log file
  DSIO::Get()->OpenLogFiles(); 
  DSIO::Get()->GetStreamLogFile() << endl;
  DSIO::Get()->GetStreamLogFile() <<"------------------------------------------" << endl;
  DSIO::Get()->GetStreamLogFile() <<"                   G4DS                   " << endl;
  DSIO::Get()->GetStreamLogFile() <<                                                endl;
  DSIO::Get()->GetStreamLogFile() <<"     The Dark Side Geant4 Simulator        " << endl;
  DSIO::Get()->GetStreamLogFile() <<"------------------------------------------" << endl;
  DSIO::Get()->GetStreamLogFile() <<  endl;
  DSIO::Get()->GetStreamLogFile() <<  endl;
  DSIO::Get()->GetStreamLogFile() <<"------------------------------------------" << endl;
  DSIO::Get()->GetStreamLogFile() <<"                   Input                  " << endl;
  DSIO::Get()->GetStreamLogFile() <<"------------------------------------------" << endl;
  DSIO::Get()->GetStreamLogFile() <<  endl;
  if (argc>1) {
    ifstream ifs(argv[1]);
    string s ;
    while (getline (ifs, s)) 
      DSIO::Get()->GetStreamLogFile() << s << endl  ;
    DSIO::Get()->GetStreamLogFile() <<  endl;
    DSIO::Get()->GetStreamLogFile() <<  endl;
    DSIO::Get()->GetStreamLogFile() <<"------------------------------------------" << endl;
    DSIO::Get()->GetStreamLogFile() <<"                   Output                 " << endl;
    DSIO::Get()->GetStreamLogFile() <<"------------------------------------------" << endl;
    DSIO::Get()->GetStreamLogFile() <<  endl;
    
    ifs.close();
  }  
  
  //---------------------------------------------------------------------------//
  
  DSLog(trace) << "Creating G4 Run Manager" << endlog;
  runManager   = DSManager::Get();     

  // Register detector geometry and materials.
  DSLog(trace) << "Creating and registering G4 geometry" << endlog;
  DSDetectorConstruction *myDetector =  new DSDetectorConstruction;
  runManager->SetUserInitialization(myDetector);
 
  // Register Geant4 physics processes
  DSLog(trace) << "Creating and registering G4  physics processes" << endlog;
  DSPhysicsList *myPhysicsList       = new DSPhysicsList;
  runManager->SetUserInitialization(myPhysicsList);
  
  
  // Register event generator.
  DSLog(trace) << "Creating and registering event generator" << endlog;
  DSPrimaryGeneratorAction *myGenerator =  new DSPrimaryGeneratorAction;
  runManager->SetUserAction(myGenerator);
   
  // Register run action. What to do at beginning and end of each run.
  DSLog(trace) << "Registering G4 run action." << endlog;
  DSRunAction *myRunAction = new DSRunAction;
  runManager->SetUserAction(myRunAction);

  // Register event action, ie. what to save/compute for each event.
  DSLog(trace) << "Registering G4 event action." << endlog;
  DSEventAction *myEventAction = new DSEventAction;
  runManager->SetUserAction(myEventAction);
 
  //Register tracking action
  DSLog(trace) << "Registering G4 tracking action." << endlog;
  DSTrackingAction *myTrackAction = new DSTrackingAction;
  runManager->SetUserAction(myTrackAction);

  // Register stepping action, ie. what to save/compute for each step.
  DSLog(trace) << "Registering G4 stepping action." << endlog;
  DSSteppingAction *myStepAction =  new DSSteppingAction;
  runManager->SetUserAction(myStepAction);

  if(IsStack) {
    // Register stacking action, ie. what to save/compute for each step.
    DSLog(trace) << "Registering G4 stacking action." << endlog;
    DSStackingAction *myStackAction =  new DSStackingAction;
    runManager->SetUserAction(myStackAction);
  }
 
  G4UIsession* session=0;

  if (argc==1) {  // Define UI session for interactive mode.
#if defined (G4UI_USE_ROOT)
    // G4URoot is a ROOT based GUI.
    session = new G4UIRoot(argc,argv);
#elif defined (G4UI_USE_XM)
    // G4UIXm is an Xm based GUI.
    session = new G4UIXm(argc,argv);
#else
    // G4UIterminal is a (dumb) terminal.
#ifdef G4UI_USE_TCSH
    session = new G4UIterminal(new G4UItcsh);
#else
    session = new G4UIterminal();
#endif
#endif
  } 

#ifdef G4VIS_USE
   // visualization manager
  G4VisManager* visManager = new DSVisManager;
  visManager->Initialize();
  G4TrajectoryDrawByParticleID* model = new G4TrajectoryDrawByParticleID;

  model->SetDefault("cyan");
  model->Set("gamma", "red");
  model->Set("neutron", "blue");
  model->Set("opticalphoton", "green");
  model->Set("e+", "magenta");
  model->Set("e-", G4Colour(0.3, 0.3, 0.3));

  visManager->RegisterModel(model);

  visManager->SelectTrajectoryModel(model->Name());
#endif



  // get the pointer to the User Interface manager 
  G4UImanager* myUI = G4UImanager::GetUIpointer();  

  if (session)  { // Define UI session for interactive mode.
    // G4UIterminal is a (dumb) terminal.
    myUI->ApplyCommand("/control/execute vis.mac");
#if defined (G4UI_USE_ROOT) || defined (G4UI_USE_XM)
    // Customize the G4UIXm menubar with a macro file :
    myUI->ApplyCommand("/control/execute visTutor/gui.mac");
#endif
    session->SessionStart();
    delete session;
  } else  {         // Batch mode
    if(argc>1 && strcmp(argv[1], "-h")) {
      G4cout << "Entering batch mode...\n";
      G4cout << "Executing script file from command line: " << argv[1] << '\n';
      G4UImanager* myUI0 = G4UImanager::GetUIpointer();
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      myUI0->ApplyCommand(command + fileName);
    } else {
      PrintUsage();
    }
  }

  
  //commands to xport the geometry in GDML
#ifdef G4LIB_USE_GDML

  if((DSStorage::Get())->GetExportGDML())
    {           
      G4VPhysicalVolume* g4wv = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume(); 

      G4GDMLParser parser;
      ifstream _ff("geo.gdml");
      if (_ff.is_open()) { 
        _ff.close();
        remove("geo.gdml");
      }
      //parser.Write("geo.gdml",g4wv,true,"http://service-spi.web.cern.ch/service-spi/app/releases/GDML/GDML_3_0_0/schema/gdml.xsd");
      parser.SetOverlapCheck (true);
      parser.Write("geo.gdml",g4wv,true);
    }
#endif
 
   // job termination
}


//---------------------------------------------------------------------------//

void PrintHeader(void)
{
  G4cout << "Dark Side Monte Carlo Simulation" << G4endl;
  G4cout << "------------------------------------------" << G4endl;
  G4cout << "Version: 1.00" << G4endl;
  G4cout << "Last Update: 21-02-2013" << G4endl;
}

//---------------------------------------------------------------------------//

void PrintUsage(void)
{
  G4cout << "Usage:" << G4endl;
  G4cout << "g4DS -h : Displays this message" << G4endl;
  G4cout << "g4DS <filename> : Executes script <filename>" << G4endl;
  G4cout << "g4DS : Executes G4DS interactively" << G4endl << G4endl;
}

//---------------------------------------------------------------------------//

/*
 * $Log: g4ds.cc,v $
 * Revision 1.2  2014/12/17 12:05:09  dfranco
 * added the possibility to set the output filename from command line (./g4ds file.mac output)
 *
 * Revision 1.1  2014/05/07 12:20:35  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.6  2013/11/04 13:20:10  dfranco
 * added colors in the visualization
 *
 * Revision 1.5  2013/07/07 09:42:22  dfranco
 * Added a default stacking for killing particles from long living nucleus decays
 *
 *
 */
