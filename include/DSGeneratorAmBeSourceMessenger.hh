#ifndef DSGeneratorAmBeSourceMessenger_h
#define DSGeneratorAmBeSourceMessenger_h 1

#include "DSIO.hh"

#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

class DSPrimaryGeneratorAction;
class G4UIcmdWithAString;
class DSGeneratorAmBeSource;

using namespace std;

class DSGeneratorAmBeSourceMessenger: public G4UImessenger {

public:
  
  //default constructor
  DSGeneratorAmBeSourceMessenger(DSGeneratorAmBeSource* );

  //destructor
  virtual ~DSGeneratorAmBeSourceMessenger();

  void SetNewValue(G4UIcommand*, G4String);

  
private:

  DSGeneratorAmBeSource*	generator;
  G4UIdirectory*		fDirectory;
  G4UIcmdWithAString*		fNeutrinoCmd;
  G4UIcmdWithAString*		fDisableParticleType;

};

#endif
