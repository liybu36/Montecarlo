// --------------------------------------------------------------------------//
/** 
 * AUTHOR: D. Franco
 * CONTACT: dfranco@in2p3.fr
*/
// --------------------------------------------------------------------------//

#ifndef _DSGENERATORMultiEvent_HH
#define _DSGENERATORMultiEvent_HH

//---------------------------------------------------------------------------//

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "DSVGenerator.hh"
#include "G4Event.hh"

using namespace std ;

struct MultiEvent {
  int    Id;
  float  BRTOT;
  int    PDG;
  float  Energy;
  float  BR;
}; 

//---------------------------------------------------------------------------//
class DSGeneratorMultiEventMessenger;
class DSGeneratorMultiEvent : public DSVGenerator {
  public:

    //default constructor
    DSGeneratorMultiEvent();

    //copy constructor
    //DSGeneratorMultiEvent(const DSGeneratorMultiEvent &);

    //destructor
    ~DSGeneratorMultiEvent();

    //public interface
    virtual void DSGeneratePrimaries(G4Event *event);



    G4int GetNumberOfParticles() const { return fPDG.size(); }
    void SetfParticles( G4ThreeVector fPart);   

    void SetListOfParticles (string val)   { fListOfStringParticles.push_back(val) ; }

  private:

    DSGeneratorMultiEventMessenger*   fMessenger;

    MultiEvent                   ConvertFromString(string);
    
    G4ParticleTable*             fParticleTable;
    G4ParticleDefinition*        fParticle;
    G4ThreeVector                fPosition;
    G4ThreeVector                fDirection;


    vector<G4int>                fPDG;
    vector<G4double>             fPDGEnergy;
    vector<G4double>             fPDGBR;

    G4bool                       fRead ;

    vector<MultiEvent>           fListOfParticles;
    vector<G4String>             fListOfStringParticles;
    G4int                        fNVertex;
    vector<G4double>             fBR;

};
#endif
/*
 * Revision 1.1  2013/03/20 09:54:28  dfranco
 * New version of g4ds
 *
 *
 */
