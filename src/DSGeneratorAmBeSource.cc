#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SPSRandomGenerator.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "DSGeneratorAmBeSource.hh"
#include "DSGeneratorAmBeSourceMessenger.hh"
#include "DSLogger.hh"

#include "DSEventHandler.hh"

//---------------------------------------------------------------------------//

DSGeneratorAmBeSource::DSGeneratorAmBeSource(): DSVGenerator("DSGeneratorAmBeSource") {

  isFirstTime = true ;
  fNeutrinoType = -1;
  fDisableParticleType=0; //neither n nor gamma disabled   
  fPosition = G4ThreeVector(0.,0.,0.) ;
  fParticleTable = G4ParticleTable::GetParticleTable();    

  fTheMessenger = new DSGeneratorAmBeSourceMessenger(this);

  DSLog(routine) << "AmBe Source Generator Built" << endlog;
 
  InitNeutron0G();
  InitNeutron1G();
  InitNeutron2G();

}
//---------------------------------------------------------------------------//


DSGeneratorAmBeSource::~DSGeneratorAmBeSource() {
  delete fTheMessenger;  
}
//---------------------------------------------------------------------------//


void DSGeneratorAmBeSource::DSGeneratePrimaries(G4Event *event) {


  if(isFirstTime) {
    if(fNeutrinoType < 0) {
      DSLog(error) << "Error: Must set the AmBeSource source type!"<<endlog;
      DSLog(fatal) << endlog;
    }  

    DSLog(routine) << "AmBeSource  " << fNeutrinoType << endlog ;

    isFirstTime = false ;
  }
 
  fPosition = GetVParticlePosition();

  // We choose one of the channels to simulate 0 Gamma, 1 Gamma, 2 Gammas 

  G4int    Scheme=0;
  //0 -NoGamma
  //1 -1Gamma
  //2 -2Gamma
 
  G4double ScRnd = G4UniformRand();

  if (ScRnd < 0.36)  Scheme=0;

  else if(ScRnd < 0.97) Scheme=1;

  else Scheme=2;


  //Here I choose different branches of the alpha + Be-9 reactions
  //full simulation 0
  //neutrons without Gamma 
  //neutrons with 1 Gamma
  //neutrons with 2 Gammas

  if (fNeutrinoType == 1) Scheme=0;
  if (fNeutrinoType == 2) Scheme=1;
  if (fNeutrinoType == 3) Scheme=2;

  if(fDisableParticleType==1 && fNeutrinoType == 1){
    DSLog(error) << "the 'neutrons without Gamma' channel is incompatible with neutrons disabled, since no particles are being generated." << endlog;
    DSLog(fatal) << endlog;
  }

  //Here we generate the Energy of the neutron, but only if it is not disabled
  if(fDisableParticleType!=1){
    G4double KinE = 0.0;

    if (Scheme==0){ KinE = ShootEnergyNeutron0G()*MeV;}
    
    if (Scheme==1){ KinE = ShootEnergyNeutron1G()*MeV;}
    
    if (Scheme==2){ KinE = ShootEnergyNeutron2G()*MeV;}
    
    fParticle = fParticleTable->FindParticle(2112);     

    G4double mass = fParticle->GetPDGMass();
    G4double energy = KinE + mass;	
    
    fDirection = GetVParticleDirection();
    
    G4double pmom = std::sqrt(energy*energy - mass*mass);
    G4double px = pmom*fDirection.x();
    G4double py = pmom*fDirection.y();
    G4double pz = pmom*fDirection.z();
    
    G4PrimaryVertex*   vertex   = new G4PrimaryVertex(fPosition,0);
    G4PrimaryParticle* particle = new G4PrimaryParticle(fParticle,px,py,pz);
    
    vertex->SetPrimary( particle );
    event->AddPrimaryVertex( vertex );    
    
    DSEventHandler::Get()->SetDId(1);
    DSEventHandler::Get()->SetDPDG(2112);
    DSEventHandler::Get()->SetDTime(0.);
    DSEventHandler::Get()->SetDPosition(fPosition);
    DSEventHandler::Get()->SetDEnergy(KinE/MeV);
    DSEventHandler::Get()->SetDDirection(fDirection);
    DSEventHandler::Get()->SetDaughters();
  }
  
  //simulate gammas, but only if they are not disabled
  if(fDisableParticleType!=2){
    if (Scheme>0) { 
      
      fParticle = fParticleTable->FindParticle(22);     
      
      G4double mass = fParticle->GetPDGMass();
      G4double energy = 4.439*MeV;
      
      fDirection = GetVParticleDirection();
      
      G4double pmom = std::sqrt(energy*energy - mass*mass);
      G4double px = pmom*fDirection.x();
      G4double py = pmom*fDirection.y();
      G4double pz = pmom*fDirection.z();
      
      G4PrimaryVertex*   vertex   = new G4PrimaryVertex(fPosition,0);
      G4PrimaryParticle* particle = new G4PrimaryParticle(fParticle,px,py,pz); 
      
      vertex->SetPrimary( particle );
      event->AddPrimaryVertex( vertex );    
      
      DSEventHandler::Get()->SetDId(2);
      DSEventHandler::Get()->SetDPosition(fPosition);
      DSEventHandler::Get()->SetDDirection(fDirection);
      DSEventHandler::Get()->SetDPDG(22);
      DSEventHandler::Get()->SetDTime(0.);
      DSEventHandler::Get()->SetDEnergy(energy/MeV);
      DSEventHandler::Get()->SetDaughters(); 
    }
    
    if (Scheme == 2) {
      
      fParticle = fParticleTable->FindParticle(22);     
      
      G4double mass = fParticle->GetPDGMass();
      G4double energy = 3.215*MeV;
      
      fDirection = GetVParticleDirection();
      
      G4double pmom = std::sqrt(energy*energy - mass*mass);
      G4double px = pmom*fDirection.x();
      G4double py = pmom*fDirection.y();
      G4double pz = pmom*fDirection.z();
      
      G4PrimaryVertex*   vertex   = new G4PrimaryVertex(fPosition,0);
      G4PrimaryParticle* particle = new G4PrimaryParticle(fParticle,px,py,pz);
      
      vertex->SetPrimary( particle );
      event->AddPrimaryVertex( vertex );    
      
      DSEventHandler::Get()->SetDId(3);
      DSEventHandler::Get()->SetDPosition(fPosition);
      DSEventHandler::Get()->SetDDirection(fDirection);
      DSEventHandler::Get()->SetDPDG(22);
      DSEventHandler::Get()->SetDTime(0.);
      DSEventHandler::Get()->SetDEnergy(energy/MeV);
      DSEventHandler::Get()->SetDaughters(); 
    }
  } //end fDisableParticleType-condition
}



//-------------------------------------------------------------------------
//    Initialize channel 'neutron without Gamma'
//-------------------------------------------------------------------------

void DSGeneratorAmBeSource::InitNeutron0G() {

   G4double sum = 0;
   G4double norma =0; 

   G4double NeutronEnergySp [58] =    
{
  0.0,
  0.0,      
  6.8	,
  5.0	,
  2.5	,
  1.5	,
  0.76	,//End 1 MeV
  0.4	,
  0.25	,
  0.1   ,//End 1.6 MeV
  0.	,
  0.	,//End 2
  0.	,
  0.	,
  0.	,
  0.	,
  0.	,
  0.	,
  0.	,
  0.	,
  0.	,
  0.	,//End 4
  0.	,
  0.	,
  0.	,
  0.	,
  0.	,
  0.	,
  0.    ,
  0.50	,//5.4-5.6VeV
  0.75  ,//End 6
  2.    ,       
  4.8   ,
  4.8	,
  3.8	,
  3.5	,
  3.75	,
  4.28	,
  4.5	,
  4.5	,
  4.25	,
  3.75	,
  3.	,
  2.2	,
  1.5	,
  1.05	,
  1.	,
  1.25	,
  2.	,
  2.25	,
  2.05	,//End 10
  1.82	,
  1.5	,
  1	,
  0.6	,
  0.2	,//End 11
  0.1	,
  0.0	
};  

  norma= 0;

  for(G4int i = 0; i< 58 ; i++) {
      norma += NeutronEnergySp[i];
   }

  for(G4int i = 0; i< 58 ; i++) {
      sum += NeutronEnergySp[i]/norma;
      Neutron0G_Probability[i]=sum;
      Neutron0G_EnergyBin[i]=G4double(i)*200;
   }


}

//-------------------------------------------------------------------------
//    Energy Shooter for Am-Be neutron without Gamma
//-------------------------------------------------------------------------

G4double DSGeneratorAmBeSource::ShootEnergyNeutron0G() {
  G4double deltaX,x,y;
  G4double val = G4UniformRand()  ;


  for(G4int i=0;i<58;i++) {
    if(Neutron0G_Probability[i] >= val) {

      if(i == 0) {  
        return Neutron0G_EnergyBin[0]/1000;
      }

      deltaX = val - Neutron0G_Probability[i] ;
      y = Neutron0G_EnergyBin[i] - Neutron0G_EnergyBin[i-1] ;
      x = Neutron0G_Probability[i] - Neutron0G_Probability[i-1] ;

      return ((deltaX*y/x + Neutron0G_EnergyBin[i])/1000) ;

    }
  }
  return 0;
}

//-------------------------------------------------------------------------
//    Initialize channel 'neutron with 1 Gamma'
//-------------------------------------------------------------------------

void DSGeneratorAmBeSource::InitNeutron1G() {

   G4double sum = 0;
   G4double norma =0; 
   G4double NeutronEnergySp [58] =    
{
  0.0	,
  0.0	,
  0.0	,
  0.0	,
  0.0	,//End 1 MeV
  0.0	, 
  0.4 	,
  1.2 	,            
  2.4	,
  4.0	,//End 2 Mev 
  5.6	,
  6.6	,
  7.1	,
  7.8	,
  8.8	,
  10	,
  11.25	,
  10.75	,
  10	,
  9.0	,
  8.4	,
  8.0  	,
  8.0	,
  8.7	,
  9.0	,
  9.3	,
  8.7	,
  8.1	,
  7.	,
  6.2	,
  4.5	,
  2.7	,
  0.2	,
  0.0,0.0,0.0,0.0,0.0,
  0.0,0.0,0.0,0.0,0.0,
  0.0,0.0,0.0,0.0,0.0,
  0.0,0.0,0.0,0.0,0.0
};  

  norma= 0;

  for(G4int i = 0; i< 58 ; i++) {
     norma += NeutronEnergySp[i];
   }

   for(G4int i = 0; i< 58 ; i++) {
     sum += NeutronEnergySp[i]/norma;
     Neutron1G_Probability[i]=sum;
     Neutron1G_EnergyBin[i]=G4double(i)*200;
   }

}
//-------------------------------------------------------------------------
//    Energy Shooter for AmBe neutron with 2 Gammas
//-------------------------------------------------------------------------

G4double DSGeneratorAmBeSource::ShootEnergyNeutron1G() {

  G4double deltaX,x,y;
  G4double val = G4UniformRand()  ;

  for(G4int i=0;i<58;i++) {
    if(Neutron1G_Probability[i] >= val) {

      if(i == 0) {
        return Neutron1G_EnergyBin[0]/1000;
      }

      deltaX = val - Neutron1G_Probability[i] ;
      y = Neutron1G_EnergyBin[i] - Neutron1G_EnergyBin[i-1] ;
      x = Neutron1G_Probability[i] - Neutron1G_Probability[i-1] ;

      return ((deltaX*y/x + Neutron1G_EnergyBin[i])/1000) ;

    }
  }
  return 0;
}

//-------------------------------------------------------------------------
//    Init channel 'neutron with 2 Gammas'
//-------------------------------------------------------------------------

void DSGeneratorAmBeSource::InitNeutron2G() {

   G4double sum = 0;
   G4double norma =0; 

   G4double NeutronEnergySp [58] =    
{
  0.0	,
  0.0	,
  0.0	,
  0.3	,
  1.3	,
  2.2	,
  2.0	,
  2.06	,
  2.25	,
  1.9	,
  1.6	,
  1.3	,
  1.25	,
  1.0	,
  0.5	,
  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
  0.0,0.0,0.0,0.0
};  

  norma= 0;

  for(G4int i = 0; i< 58 ; i++) {
     norma += NeutronEnergySp[i];
  }

  for(G4int i = 0; i< 58 ; i++) {
     sum += NeutronEnergySp[i]/norma;
     Neutron2G_Probability[i]=sum;
     Neutron2G_EnergyBin[i]=G4double(i)*200;
  }


}

//-------------------------------------------------------------------------
//    Energy Shooter for AmBe neutron with 2 Gammas
//-------------------------------------------------------------------------

G4double DSGeneratorAmBeSource::ShootEnergyNeutron2G() {

  G4double deltaX,x,y;
  G4double val = G4UniformRand()  ;

  for(G4int i=0;i<58;i++) {
    if(Neutron2G_Probability[i] >= val) {

      if(i == 0) {
        return Neutron2G_EnergyBin[0]/1000;
      }

      deltaX = val - Neutron2G_Probability[i] ;
      y = Neutron2G_EnergyBin[i] - Neutron2G_EnergyBin[i-1] ;
      x = Neutron2G_Probability[i] - Neutron2G_Probability[i-1] ;

      return ((deltaX*y/x + Neutron2G_EnergyBin[i])/1000) ;

    }
  }
  return 0;
}




