
#ifndef DSStackingRDMChainMessenger_h
#define DSStackingRDMChainMessenger_h 1
#include "DSIO.hh"

#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4ParticleGun.hh"


class G4UIcmdWithADoubleAndUnit;
class DSStackingRDMChain;

using namespace std;

class DSStackingRDMChainMessenger: public G4UImessenger {

  public:
    DSStackingRDMChainMessenger(DSStackingRDMChain*  );
   ~DSStackingRDMChainMessenger();

  void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    DSStackingRDMChain*                  stacking;   
    G4UIdirectory*                       fDirectory;
    G4UIcmdWithADoubleAndUnit*           fLifeTimeCmd;
};

#endif
