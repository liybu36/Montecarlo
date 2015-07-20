// --------------------------------------------------------------------------//
/** 
 * AUTHOR: A. Meregaglia
*/
// --------------------------------------------------------------------------//

#ifndef _DSGENERATOREnergyDeposit_HH
#define _DSGENERATOREnergyDeposit_HH

//---------------------------------------------------------------------------//

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "DSVGenerator.hh"
#include "G4Event.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSPosDistribution.hh"

using namespace std ;


//---------------------------------------------------------------------------//
class DSGeneratorEnergyDepositMessenger;
class DSGeneratorEnergyDeposit : public DSVGenerator {
  public:

    //default constructor
    DSGeneratorEnergyDeposit();


    //destructor
    ~DSGeneratorEnergyDeposit();

    //public interface
    virtual void DSGeneratePrimaries(G4Event *event);



  void SetFile (string val)   { fFile=val ; }

  private:




  DSGeneratorEnergyDepositMessenger*   fMessenger;
  G4SPSAngDistribution*                fSPSAng;
  G4SPSPosDistribution*                fSPSPos;
  
    
  G4ParticleTable*             fParticleTable;
  G4ParticleDefinition*        fParticle;
  G4ThreeVector                fPosition;
  G4ThreeVector                fDirection;

  G4bool                       fRead ;
  G4int                        fnumev;
  G4String                     fFile; 
};
#endif
/*
 * Revision 1.1  2013/03/20 09:54:28  dfranco
 * New version of g4ds
 *
 *
 */
