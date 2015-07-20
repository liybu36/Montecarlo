// --------------------------------------------------------------------------//
/** 
 * AUTHOR: D. Franco
 * CONTACT: davide.franco@mi.infn.it
*/
// --------------------------------------------------------------------------//

#ifndef _DSGENERATORG4GUN_HH
#define _DSGENERATORG4GUN_HH

//---------------------------------------------------------------------------//


#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include "DSVGenerator.hh"
#include "DSGeneratorG4GunMessenger.hh"
#include "G4Event.hh"
#include "G4SPSPosDistribution.hh"
#include "G4SPSAngDistribution.hh"
#include "G4SPSEneDistribution.hh"

//---------------------------------------------------------------------------//

class DSGeneratorG4Gun : public DSVGenerator {
public:

  //default constructor
  DSGeneratorG4Gun();

  //copy constructor
  //DSGeneratorG4Gun(const DSGeneratorG4Gun &);

  //destructor
  virtual ~DSGeneratorG4Gun();

  //public interface
  virtual void DSGeneratePrimaries(G4Event *event);



  //void SetParticlePosition(G4ThreeVector pos){ fParticleGun->SetParticlePosition(pos); }
  //void SetParticleDefinition(G4ParticleDefinition* part) { fParticleGun->SetParticleDefinition(part); }
  //void SetParticleMomentumDirection(G4ParticleMomentum mom){ fParticleGun->SetParticleMomentumDirection(mom); }
  //void SetParticleEnergy(G4double ene){ fParticleGun->SetParticleEnergy(ene); }
  

private:

  DSGeneratorG4GunMessenger*   fTheMessenger;
  G4ParticleGun*               fParticleGun;
};
#endif
/*
 * $Log: DSGeneratorG4Gun.hh,v $
 * Revision 1.1  2014/05/07 12:20:52  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
