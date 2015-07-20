//---------------------------------------------------------------------------//

#ifndef DSGeneratorSCSMessenger_HH
#define DSGeneratorSCSMessenger_HH

//---------------------------------------------------------------------------//


#include <G4UIdirectory.hh>
#include <G4UImessenger.hh>
#include <G4UIcmdWithAString.hh>

#include "DSGeneratorSCS.hh"


//---------------------------------------------------------------------------//


class DSGeneratorSCS;


class DSGeneratorSCSMessenger: public G4UImessenger {


  public:
    DSGeneratorSCSMessenger( DSGeneratorSCS* );
   ~DSGeneratorSCSMessenger();

    void SetNewValue(G4UIcommand*, G4String);


  private:
    DSGeneratorSCS*       fGenerator;
    G4UIdirectory*        fDirectory;
    G4UIcmdWithAString*   fIsotopeCmd;    


};

#endif
