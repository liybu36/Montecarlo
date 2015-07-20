#ifndef DSPrimaryGeneratorActionMessenger_h
#define DSPrimaryGeneratorActionMessenger_h


#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

class DSPrimaryGeneratorAction;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;


class DSPrimaryGeneratorActionMessenger: public G4UImessenger {

  public:
    DSPrimaryGeneratorActionMessenger(DSPrimaryGeneratorAction*);
   ~DSPrimaryGeneratorActionMessenger();

  void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    DSPrimaryGeneratorAction*   fGeneratorPrimary;   
    G4UIdirectory*              fDirectory;
    G4UIcmdWithAString*         fSelectCmd;  
    
};

#endif
/*
 * $Log: DSPrimaryGeneratorActionMessenger.hh,v $
 * Revision 1.1  2014/05/07 12:20:54  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
