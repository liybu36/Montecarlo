// --------------------------------------------------------------------------//
/** 
 * AUTHOR: D. Franco
 * CONTACT: davide.franco@mi.infn.it
*/
// --------------------------------------------------------------------------//

#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SPSRandomGenerator.hh"

#include "DSLogger.hh"
#include "DSGeneratorG4Gun.hh"
#include "DSGeneratorG4GunMessenger.hh"

#include "DSEventHandler.hh"
//---------------------------------------------------------------------------//


DSGeneratorG4Gun::DSGeneratorG4Gun(): DSVGenerator("DSGeneratorG4Gun") {

  fParticleGun  = new G4ParticleGun ;
  

  fTheMessenger = new DSGeneratorG4GunMessenger(this);
  if(!fParticleGun) {
    DSLog(error) << "Could not allocate G4ParticleGun! Out of memory?"<<endlog;
    DSLog(fatal) << endlog;
  }
  DSLog(routine) << "G4ParticleGun Constructed." << endlog;
}
//---------------------------------------------------------------------------//

//DSGeneratorG4Gun::DSGeneratorG4Gun(const DSGeneratorG4Gun & other)
//{;}

//---------------------------------------------------------------------------//

DSGeneratorG4Gun::~DSGeneratorG4Gun()
{
  delete fTheMessenger;
  delete fParticleGun;
  
}

//---------------------------------------------------------------------------//


void DSGeneratorG4Gun::DSGeneratePrimaries(G4Event *event) {
  
  //fParticleGun->SetNumberOfParticles(1000);
  
  
  
  
  G4ThreeVector position  =  GetVParticlePosition() ;
  for(int i=0; i<GetVNumberOfParticles(); ++i) {

    G4ThreeVector direction =  GetVParticleDirection() ;
    G4double      energy    =  GetVParticleEnergy(GetVParticleDefinition());


    fParticleGun->SetParticleDefinition(GetVParticleDefinition());
    fParticleGun->SetParticlePosition(position);
    fParticleGun->SetParticleMomentumDirection(direction);     
    fParticleGun->SetParticleEnergy(energy);     
    if(GetVParticleCharge() != 0.) fParticleGun->SetParticleCharge(GetVParticleCharge());

    DSLog(development)<< fParticleGun->GetParticleDefinition()->GetParticleName()  << " " 
                      << "Energy: "    << fParticleGun->GetParticleEnergy()/MeV    << " MeV; "
                      << "Position: "  << fParticleGun->GetParticlePosition()/cm   << " cm;  "
                      << "Direction: " << fParticleGun->GetParticleMomentumDirection()   << " "
                      << endlog;

    fParticleGun->GeneratePrimaryVertex(event);
    
    if(GetVNumberOfParticles() == 1) {
      DSEventHandler::Get()->SetPDG(GetVParticleDefinition()->GetPDGEncoding());
      DSEventHandler::Get()->SetPosition(position);
      DSEventHandler::Get()->SetDirection(direction);
      DSEventHandler::Get()->SetEnergy(energy);
    
    } else {
      if(i ==0 ) {
	DSEventHandler::Get()->SetPDG(GetVParticleDefinition()->GetPDGEncoding());
	DSEventHandler::Get()->SetPosition(position);
	DSEventHandler::Get()->SetDirection(direction);
	DSEventHandler::Get()->SetEnergy(energy);
      
      } else {
	DSEventHandler::Get()->SetDPDG(GetVParticleDefinition()->GetPDGEncoding());
	DSEventHandler::Get()->SetDPosition(position);
	DSEventHandler::Get()->SetDDirection(direction);
	DSEventHandler::Get()->SetDEnergy(energy);
	DSEventHandler::Get()->SetDaughters();
      
      
      }
    }
  }

}


/*
 * $Log: DSGeneratorG4Gun.cc,v $
 * Revision 1.3  2014/07/18 13:54:49  dfranco
 * Added a new position generator in the Gas Pocket; calibrated the numnber of photons generated per electron in S2; added a new command to generate multiple particles in the same position
 *
 * Revision 1.2  2014/05/21 10:28:08  dfranco
 * added the possibility to shoot ions
 *
 * Revision 1.1  2014/05/07 12:21:02  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.3  2014/03/06 11:11:34  dfranco
 * fixed a bug when saving the primary particle information
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
