
#ifndef DSGeneratorRDMDecayGunMessenger_h
#define DSGeneratorRDMDecayGunMessenger_h 1
#include "DSIO.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4ParticleGun.hh"
#include "DSGeneratorRDMUIcmdWithNucleusAndUnit.hh"

class DSPrimaryGeneratorAction;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class DSGeneratorRDMDecayGun;

using namespace std;

class DSGeneratorRDMDecayGunMessenger: public G4UImessenger {

  public:
    DSGeneratorRDMDecayGunMessenger(DSGeneratorRDMDecayGun*  );
   ~DSGeneratorRDMDecayGunMessenger();

  void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    DSGeneratorRDMDecayGun*              generator;   
    G4UIdirectory*                       fDirectory;
    DSGeneratorRDMUIcmdWithNucleusAndUnit*      fIonCmd;
};

#endif


/*
 * Revision 1.1  2013/03/22 13:50:09  dfranco
 * Added radioactive decay and radioactive decay chain generators. They require the correspondent stacking actions. Two mac files are included as examples
 *
 */
