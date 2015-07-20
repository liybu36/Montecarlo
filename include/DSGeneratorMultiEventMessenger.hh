#ifndef DSGeneratorMultiEventMessenger_HH
#define DSGeneratorMultiEventMessenger_HH 1

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

#include "DSIO.hh"
#include "DSGeneratorMultiEvent.hh"


class DSGeneratorMultiEvent;


class DSGeneratorMultiEventMessenger: public G4UImessenger {

  public:
    DSGeneratorMultiEventMessenger(DSGeneratorMultiEvent*  );
   ~DSGeneratorMultiEventMessenger();

    void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    DSGeneratorMultiEvent*               generator;   
    G4UIdirectory*                       fDirectory;
    G4UIcmdWithAString*                  fPDGCmd;
     
    //MultiEvent  			 fMultiParticle;
    //MultiEvent& 			 ConvertFromString(string);
    
};

#endif
/*
 * Revision 1.1  2013/03/20 09:54:28  dfranco
 * New version of g4ds
 *
 *
 */
