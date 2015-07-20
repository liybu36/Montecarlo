// --------------------------------------------------------------------------//
/** 
 * AUTHOR: D. Franco
 * CONTACT: davide.franco@mi.infn.it
*/
// --------------------------------------------------------------------------//

#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ParticleTable.hh"
#include "G4SPSRandomGenerator.hh"

#include "DSLogger.hh"
#include "DSEventHandler.hh"
#include "DSGeneratorMultiEvent.hh"
#include "DSGeneratorMultiEventMessenger.hh"
//---------------------------------------------------------------------------//


DSGeneratorMultiEvent::DSGeneratorMultiEvent(): DSVGenerator("DSGeneratorMultiEvent") {

  fRead  = false ;
  fNVertex = 0;
  
  fPosition = G4ThreeVector(0.,0.,0.) ;
  
  fParticleTable = G4ParticleTable::GetParticleTable();  ;
  fMessenger = new DSGeneratorMultiEventMessenger(this);
  
  DSLog(routine) << "MultiEvent Generatore Built" << endlog;
}
//---------------------------------------------------------------------------//

//DSGeneratorMultiEvent::DSGeneratorMultiEvent(const DSGeneratorMultiEvent & other)
//{;}

//---------------------------------------------------------------------------//

DSGeneratorMultiEvent::~DSGeneratorMultiEvent()
{
  delete fMessenger;
}

//---------------------------------------------------------------------------//


void DSGeneratorMultiEvent::DSGeneratePrimaries(G4Event *event) {
  if(!fRead) {
  
    
    fRead = true ;


    cout << fListOfStringParticles.size() << endl ;
    // Convert strings into particle properties
    for(int i=0;i<(G4int) fListOfStringParticles.size();i++) 
      fListOfParticles.push_back(ConvertFromString(fListOfStringParticles[i]));
    
    DSLog(routine)<< "Generated particles in each vertex:" << endlog;
    
    if(!int(fListOfParticles.size())) {
      DSLog(error) << "Number of selected particles equals to 0"<<endlog;
      DSLog(fatal) << endlog;  
    }

    // Branching ratio normalizer
    G4int Id = 0;
    G4double BRNorm = 0;
    for(G4int i=0; i< (G4int) fListOfParticles.size(); i++) { 
      if(fListOfParticles[i].Id > Id) {
        Id = fListOfParticles[i].Id;
        BRNorm += fListOfParticles[i].BRTOT;
      }   
    }
        
    // Build the cumulative probability ditribution
    Id = 0;    
    G4double BRCum = 0;
    for(G4int i=0; i< (G4int) fListOfParticles.size(); i++) {
      G4int pdg =  fListOfParticles[i].PDG;
      DSLog(routine) << i << " " 
                     << fParticleTable->FindParticle(pdg)->GetParticleName() << " "
                     << fListOfParticles[i].Energy << " MeV"
                     << endlog;
      if(fListOfParticles[i].Id > Id) {
        fNVertex++; 
        Id = fListOfParticles[i].Id;
        BRCum += fListOfParticles[i].BRTOT;
        fBR.push_back(BRCum/BRNorm);
      }
    }

  }


  // Choose the primary event
  G4int source = 0;
  G4double rand = G4UniformRand() ;  
  for(G4int i=0;i<fNVertex;i++) {
    if(rand < fBR[i]) {
      source = i ;
      break;
    }
  }
  
  // Get fixed or random position from user commands
  fPosition = GetVParticlePosition();
 
  for(G4int i = 0; i <G4int(fListOfParticles.size()); i++) {
    if(fListOfParticles[i].Id != source + 1) continue ;
    if(G4UniformRand() > fListOfParticles[i].BR) continue;  
 
    const G4int pdg = fListOfParticles[i].PDG ;
    fParticle = fParticleTable->FindParticle(pdg);     

    const G4double mass = fParticle->GetPDGMass();
    const G4double energy = fListOfParticles[i].Energy*MeV  + mass;

    fDirection = GetVParticleDirection();
    
    G4double pmom = std::sqrt(energy*energy - mass*mass);
    G4double px = pmom*fDirection.x();
    G4double py = pmom*fDirection.y();
    G4double pz = pmom*fDirection.z();

    G4PrimaryVertex*   vertex   = new G4PrimaryVertex(fPosition,0);
    G4PrimaryParticle* particle = new G4PrimaryParticle(fParticle,px,py,pz);

    vertex->SetPrimary( particle );
    event->AddPrimaryVertex( vertex );    
    
    DSEventHandler::Get()->SetEnergy(DSEventHandler::Get()->GetEnergy() + fListOfParticles[i].Energy*MeV);
    DSEventHandler::Get()->SetPosition(fPosition);
   
    DSEventHandler::Get()->SetDId(source+1);     
    //    DSEventHandler::Get()->SetDPID(0);
    DSEventHandler::Get()->SetDParentTrackID(0);
    DSEventHandler::Get()->SetDPosition(fPosition);  
    DSEventHandler::Get()->SetDDirection(fDirection); 
    DSEventHandler::Get()->SetDPDG(pdg);
    DSEventHandler::Get()->SetDTime(0.);
    DSEventHandler::Get()->SetDEnergy(fListOfParticles[i].Energy*MeV);
    DSEventHandler::Get()->SetDaughters();  
  }
  if(!event->GetNumberOfPrimaryVertex ()) {
    const G4int pdg =  12;
    fParticle = fParticleTable->FindParticle(pdg);     

    const G4double mass = fParticle->GetPDGMass();
    const G4double energy = mass;

    fDirection = GetVParticleDirection();
    
    G4double pmom = std::sqrt(energy*energy - mass*mass);
    G4double px = pmom*fDirection.x();
    G4double py = pmom*fDirection.y();
    G4double pz = pmom*fDirection.z();

    G4PrimaryVertex*   vertex   = new G4PrimaryVertex(fPosition,0);
    G4PrimaryParticle* particle = new G4PrimaryParticle(fParticle,px,py,pz);

    vertex->SetPrimary( particle );
    event->AddPrimaryVertex( vertex );    
  }
}

void DSGeneratorMultiEvent::SetfParticles( G4ThreeVector fPart) {
  fPDG.push_back(G4int(fPart.x()));       // PDG Code number 
  fPDGEnergy.push_back(fPart.y()*MeV);    // Energy in MeV
  fPDGBR.push_back(fPart.z());            // Branching Ratio
}   

MultiEvent DSGeneratorMultiEvent::ConvertFromString(string text) {  
  MultiEvent myMultiParticle;
  sscanf (text.c_str(),"%d %f %d %f %f",
          &myMultiParticle.Id,
          &myMultiParticle.BRTOT,
          &myMultiParticle.PDG,
          &myMultiParticle.Energy,
          &myMultiParticle.BR
  );
  return myMultiParticle;

}


/*
 * Revision 1.1  2013/03/20 09:54:28  dfranco
 * New version of g4ds
 *
 *
 */
