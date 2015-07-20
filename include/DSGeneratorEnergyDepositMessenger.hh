#ifndef DSGeneratorEnergyDepositMessenger_HH
#define DSGeneratorEnergyDepositMessenger_HH 1

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

#include "DSIO.hh"
#include "DSGeneratorEnergyDeposit.hh"


class DSGeneratorEnergyDeposit;


class DSGeneratorEnergyDepositMessenger: public G4UImessenger {

  public:
    DSGeneratorEnergyDepositMessenger(DSGeneratorEnergyDeposit*  );
   ~DSGeneratorEnergyDepositMessenger();

    void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    DSGeneratorEnergyDeposit*            generator;   
    G4UIdirectory*                       fDirectory;
    G4UIcmdWithAString*                  fFileCmd;
     
     
};

#endif
/*
 * Revision 1.1  2013/03/20 09:54:28  dfranco
 * New version of g4ds
 *
 *
 */
