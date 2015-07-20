//---------------------------------------------------------------------------//

#ifndef _DSGENERATORSPS_HH
#define _DSGENERATORSPS_HH

//---------------------------------------------------------------------------//


#include <G4Event.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include <G4ThreeVector.hh>

#include "DSGeneratorSCSMessenger.hh"
#include "DSVGenerator.hh"


//---------------------------------------------------------------------------//

class DSGeneratorSCSMessenger;

class DSGeneratorSCS: public DSVGenerator {


  public:

    DSGeneratorSCS();
   ~DSGeneratorSCS();

    virtual void DSGeneratePrimaries( G4Event* ); 

    void SetIsotope( G4String val ) { fIsotope = val; }

  private:

    void          LoadAr39CrossSection();
    G4double      PickAnEnergy();


    DSGeneratorSCSMessenger*  fMessenger;

    G4String                  fIsotope;
    G4bool                    fIsFirstEvent;
    G4ParticleDefinition*     fParticle;  
    G4ParticleTable*          fParticleTable;

    G4ThreeVector             fDirection;
    G4ThreeVector             fPosition;

    G4double                  fSpectrumEne[1000];
    G4double                  fSpectrumCumul[1000];
    G4int                     fSamples;

};

#endif
