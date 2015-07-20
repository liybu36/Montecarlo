#ifndef DSStackingActionMessenger_h
#define DSStackingActionMessenger_h


#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

#include "DSStackingAction.hh"

class G4UserStackingAction;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class DSStackingAction;
class G4StackManager;


class DSStackingActionMessenger: public G4UImessenger {

  public:
    DSStackingActionMessenger(DSStackingAction*);
   ~DSStackingActionMessenger();

  void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    G4UIdirectory*              fDirectory;
    G4UIcmdWithAString*         fSelectCmd;  
    DSStackingAction*           fStacking;
};

#endif


/*
 * $Log: DSStackingActionMessenger.hh,v $
 * Revision 1.1  2014/05/07 12:20:55  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.4  2013/07/07 09:42:20  dfranco
 * Added a default stacking for killing particles from long living nucleus decays
 *
 *
 */
