#ifndef DSStackingRDMMessenger_h
#define DSStackingRDMMessenger_h 1

#include "DSIO.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4ParticleGun.hh"
#include "G4UIcmdWithAnInteger.hh"

class G4UIcmdWithADoubleAndUnit;
class DSStackingRDM;

using namespace std;

class DSStackingRDMMessenger: public G4UImessenger {

  public:
    DSStackingRDMMessenger(DSStackingRDM*  );
   ~DSStackingRDMMessenger();

  void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    DSStackingRDM*                       stacking;   
    G4UIdirectory*                       fDirectory;
    G4UIcmdWithAnInteger*                fKillParticleCmd;
    G4UIcommand*                         fKillLEParticleCmd; 
   
    
};

#endif
