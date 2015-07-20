//---------------------------------------------------------------------------//
/**                                                            
 *      
 * CLASS DECLARATION:  DSVGenerator.cc
 *
 *---------------------------------------------------------------------------//
 *
 * DESCRIPTION: 
 *
 */ 
// Begin description of class here
/**
 *
 * Pure virtual base class for DS generators. 
 * 
 */
// End class description
//
/**  
 * SPECIAL NOTES:
 *
 */
// 
// --------------------------------------------------------------------------//
/** 
 * AUTHOR: davide.franco@mi.infn.it
 */
// --------------------------------------------------------------------------//

#include "DSVGenerator.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "DSLogger.hh"
#include "DSVGeneratorMessenger.hh"
#include "DSStorage.hh"
#include "G4SPSRandomGenerator.hh"
#include "G4TransportationManager.hh"
#include "G4Material.hh"
#include "Randomize.hh"
#include <fstream>
//---------------------------------------------------------------------------//

DSVGenerator::DSVGenerator(const G4String &myname):
  fGeneratorName(myname), fG4Messenger(0), fReportingFrequency(1000) {
  fVolumeFlag = false ;
  fEnergyFlag = false ;
  fNumberOfHits = 0;
  G4ThreeVector zero(0., 0., 0.) ;
  fRndGen = new G4SPSRandomGenerator;
  IsMyEnergyDistribution = false ;
  fCharge = 0.0;
  fNumberOfParticles = 1 ;
  gNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  fPositionFlag = 1;

  IsG2 = false;
  //SetRadius0(0.);
  //SetCentreCoords(zero);    

  IsVolumetric = false ;
  fMessenger  = new DSVGeneratorMessenger(this);
}

//---------------------------------------------------------------------------//

//DSVGenerator::DSVGenerator(const DSVGenerator & other)
//{;}

//---------------------------------------------------------------------------//

DSVGenerator::~DSVGenerator(){  
  delete fSPSPos;
  delete fSPSAng;
  delete fSPSEne;
}

G4ThreeVector DSVGenerator::GetVParticleDirection() { 

  if( fPositionFlag != 1 ) {
//  if(fVolumeFlag) {
    fSPSAng->SetAngDistType("iso"); 
    fSPSAng->SetPosDistribution(fSPSPos); 
    return fSPSAng->GenerateOne(); 
  }
  
  return fDirection;
}

G4ThreeVector DSVGenerator::GetVParticlePosition() { 

  G4ThreeVector thePos( 0., 0., 0.);

  if( fPositionFlag == 1 ) return fPosition;
  else if( fPositionFlag == 2 ) return fSPSPos->GenerateOne(); 
  else if( fPositionFlag == 3 ) thePos = GenerateInCryostats();
  else if( fPositionFlag == 4 ) thePos = GenerateInTeflon();
  else if( fPositionFlag == 5 ) thePos = GenerateInFusedSilica();
  else if( fPositionFlag == 6 ) thePos = GenerateInPMTPhotocathode();
  else if( fPositionFlag == 7 ) thePos = GenerateInPMTStem();
  else if( fPositionFlag == 8 ) thePos = GenerateInLiquidArgon();
  else if( fPositionFlag == 9 ) thePos = GenerateInHolderSource();
  return thePos ; 
}


G4double DSVGenerator::GetVParticleEnergy(G4ParticleDefinition *pDef) { 
  if(fEnergyFlag) 
    return fSPSEne->GenerateOne(pDef); 
  if(IsMyEnergyDistribution) {
    double val = G4UniformRand();
    for(int i=1;i< int(fMyProb.size());++i) {
      if(fMyProb[i] > val) {
        double x2 = fMyProb[i];
	double x1 = fMyProb[i-1];
        double y2 = fMyEne[i];
	double y1 = fMyEne[i-1];
	if(x2 - x1 == 0) continue ;
	double _m = (y2 - y1)/(x2 - x1);
	double q = y2 - x2*_m;
	return _m*val + q ;
	
      }
    }
  
  
  }
  
  return fEnergy;
}

    
void DSVGenerator::SetIsVolumeDistribution(G4bool val) {
  fVolumeFlag = val;
  if(fVolumeFlag && !IsVolumetric) {
    DSLog(routine) << " Random spatial distribution activated " << endlog ;
    fSPSAng = new G4SPSAngDistribution() ;
    fSPSPos = new G4SPSPosDistribution() ;
    fSPSPos->SetBiasRndm(fRndGen);
    fSPSAng->SetBiasRndm(fRndGen);  
    IsVolumetric = true ;
  }
}

void DSVGenerator::SetIsEnergyDistribution(G4bool val) {
  fEnergyFlag = val;
  if(fEnergyFlag) {
    fSPSEne = new G4SPSEneDistribution ;
    fSPSEne->SetBiasRndm(fRndGen);
  }
}

void DSVGenerator::SetEnergyFileName(string filename) {
  IsMyEnergyDistribution = true ;
  int count = 0;
  double x,y;
  double prob = 0 ;
  ifstream ff(filename.c_str());
  while(!ff.eof()) {
    ff >> x >> y ;
    if(ff.eof()) break ;
    prob += y ;
    fMyEne.push_back(x*keV);
    fMyProb.push_back(prob);
    count++;
  }
  ff.close();
  for(int i=0;i<count;++i) fMyProb[i] /= prob ;

}

void DSVGenerator::SetUniformTPC() { 
  fPositionFlag = 2;
  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");
  SetCentreCoords(G4ThreeVector(0, 0, -36.5*mm)) ;
  SetRadius(200.*mm);
  SetHalfZ(200.*mm); 
  
  if (DSStorage::Get()->Get5KGeometry()) {
    SetCentreCoords(G4ThreeVector(0, 0,0*mm)) ;
    SetRadius(840.*mm);
    SetHalfZ(908.*mm); 
  }
  if (DSStorage::Get()->Get20KGeometry()) {
    SetCentreCoords(G4ThreeVector(0, 0,0*mm)) ;
    SetRadius(160.*cm);
    SetHalfZ(1200.*mm); 
  }
  
  ConfineSourceToVolume("ActiveLAr");
}

void DSVGenerator::SetUniformGasPocket() { 
  fPositionFlag = 2;
  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");
  SetCentreCoords(G4ThreeVector(0, 0, +150.0*mm)) ;
  SetRadius(200.*mm);
  SetHalfZ(6.*mm); 
  if (DSStorage::Get()->Get5KGeometry()) {
    SetCentreCoords(G4ThreeVector(0, 0,905*mm)) ;
    SetRadius(840.*mm);
    SetHalfZ(10.*mm); 
  }
  ConfineSourceToVolume("GasPocket");
}

void DSVGenerator::SetTPCCenter() { 
  SetVParticlePosition(G4ThreeVector(0, 0, -36.54*mm));
  if(fVolumeFlag) SetCentreCoords(G4ThreeVector(0, 0, -36.54*mm)) ;
  if (DSStorage::Get()->Get5KGeometry()) SetCentreCoords(G4ThreeVector(0, 0, 0)) ;
}



G4ThreeVector DSVGenerator::GenerateInCryostats(){

  G4double center_z = IsG2 ? 0.*mm   : 161.5*mm;
  G4double height   = IsG2 ? 2231*mm : 1400.*mm;
  G4double radius   = IsG2 ? 910.*mm : 356.*mm;  

  if (DSStorage::Get()->Get5KGeometry()) {
    center_z = 130.*mm;
    height   = 2774*mm ; 
    radius   = 996*mm ; 
  }

  if (DSStorage::Get()->Get20KGeometry()) {
    center_z = 0.*mm;
    height   = 4200*mm ; 
    radius   = 1830*mm ; 
  }
  
  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");
  SetCentreCoords(G4ThreeVector(0,0,center_z));
  SetHalfZ( height/2. );
  SetRadius( radius );

  G4ThreeVector myPos(0.,0.,0.);
  while( !CheckMaterial( myPos, "StainlessSteel" ) ) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInLiquidArgon(){

  G4double center_z = IsG2 ? 0.*mm   : 161.5*mm;
  G4double height   = IsG2 ? 2231*mm : 1400.*mm;
  G4double radius   = IsG2 ? 910.*mm : 356.*mm;  

  if (DSStorage::Get()->Get5KGeometry()) {
    center_z = 0*mm;
    height   = 1820*mm;
    radius   = 990*mm;
  }
  if (DSStorage::Get()->Get20KGeometry()) {
    center_z = 0*mm;
    height   = 2500*mm;
    radius   = 1600*mm;
  }

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");
  SetCentreCoords(G4ThreeVector(0,0,center_z));
  SetHalfZ( height/2. );
  SetRadius( radius );

  G4ThreeVector myPos(0.,0.,0.);
  myPos = fSPSPos->GenerateOne();
  while( !CheckMaterial( myPos, "LiquidArgon" ) ) myPos = fSPSPos->GenerateOne();
  
  //std::cout<<myPos.x()<<" "<<myPos.y()<<" "<<myPos.z()<<std::endl;

  return myPos;
}



G4ThreeVector DSVGenerator::GenerateInTeflon(){

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");
  SetCentreCoords(G4ThreeVector(0,0,160.*mm));
  SetHalfZ(860./2.*mm);
  SetRadius(235.*mm);
  if (DSStorage::Get()->Get5KGeometry()) {
    SetCentreCoords(G4ThreeVector(0,0,0));
    SetHalfZ(905*mm);
    SetRadius(867.*mm);
  }
  if (DSStorage::Get()->Get20KGeometry()) {
    SetCentreCoords(G4ThreeVector(0,0,0));
    SetHalfZ(1200*mm);
    SetRadius(1500.*mm);
  }


  G4ThreeVector myPos(0.,0.,0.); 
  myPos = fSPSPos->GenerateOne();
  while( !CheckMaterial( myPos, "Teflon" ) ) myPos = fSPSPos->GenerateOne();
 
  return myPos;
}


G4ThreeVector DSVGenerator::GenerateInFusedSilica(){

  G4double pos_z  = 160.5 * mm ;
  G4double half_z = 10./2. * mm ;
  G4double radius = 213.*mm ; 
  if (DSStorage::Get()->Get5KGeometry()) {
    pos_z  = 913*mm ; 
    half_z = 13. * mm ;
    radius = 840.*mm ; 
  }
  if (DSStorage::Get()->Get20KGeometry()) {
    pos_z  = 1200*mm ; 
    half_z = 40. * mm ;
    radius = 1500.*mm ; 
  }
  
  if( G4UniformRand() > 0.5 ){
    if (DSStorage::Get()->Get5KGeometry()) {
      pos_z  = -903*mm ; 
    } else if (DSStorage::Get()->Get20KGeometry()) {
      pos_z  = 1200*mm ; 
    } else {
    pos_z = -221. * mm;
    half_z = 10./2. * mm;
    }
  }
  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");

  SetCentreCoords( G4ThreeVector(0,0,pos_z) );
  SetHalfZ( half_z );
  SetRadius(radius);


  G4ThreeVector myPos(0.,0.,0.);
  myPos = fSPSPos->GenerateOne();
  while( !CheckMaterial( myPos, "FusedSilica" ) ) myPos = fSPSPos->GenerateOne();
 
  return myPos;
}


G4ThreeVector DSVGenerator::GenerateInPMTPhotocathode(){

  G4double pos_z  = 170.5 * mm ;
  G4double half_z = 12./2. * mm ;
  if( G4UniformRand() > 0.5 ){
    pos_z = -220.5 * mm;
    half_z = 12./2. * mm;
  }

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");
  SetCentreCoords( G4ThreeVector(0,0,pos_z) );
  SetHalfZ( half_z );
  SetRadius(235.*mm);

  G4ThreeVector myPos(0.,0.,0.);
  while( !CheckMaterial( myPos, "Bialkali" ) ) myPos = fSPSPos->GenerateOne();
 
  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInPMTStem(){

  G4double pos_z  = 287.6 * mm ;
  G4double half_z = 1. * mm ;
  if( G4UniformRand() > 0.5 ){
    pos_z = -344.86 * mm;
    half_z = 1. * mm;
  }

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Cylinder");
  SetCentreCoords( G4ThreeVector(0,0,pos_z) );
  SetHalfZ( half_z );
  SetRadius(235.*mm);

  G4ThreeVector myPos(0.,0.,0.);
  while( !CheckMaterial( myPos, "FusedSilica" ) ) myPos = fSPSPos->GenerateOne();
 
  return myPos;
}

G4ThreeVector DSVGenerator::GenerateInHolderSource(){

  G4double center_r   = 3.5*mm + (DSStorage::Get()->GetHolderRadius());
  G4double center_phi = DSStorage::Get()->GetHolderPhi();
  G4double center_z   = DSStorage::Get()->GetHolderZ();
  G4double radius     = 15.*mm;  

  SetIsVolumeDistribution(true);
  SetPosDisType("Volume");
  SetPosDisShape("Sphere");
  SetCentreCoords(G4ThreeVector( center_r*(std::cos(center_phi)), center_r*(std::sin(center_phi)), center_z));
  SetRadius( radius );

  G4ThreeVector myPos(0.,0.,0.);
  while( !CheckMaterial( myPos, "Teflon" ) ) myPos = fSPSPos->GenerateOne();

  return myPos;
}

G4bool DSVGenerator::CheckMaterial(G4ThreeVector pos, G4String MatName) {
  
  G4bool found = false ;
  fNumberOfHits++;
  DSStorage::Get()->SetNumberOfHits( fNumberOfHits );
 
  G4ThreeVector null(0.,0.,0.);
  G4ThreeVector *ptr;
  ptr = &null;

  G4VPhysicalVolume *theVolume;
  theVolume = gNavigator->LocateGlobalPointAndSetup(pos,ptr,true);
  G4Material *amat = theVolume->GetLogicalVolume()->GetMaterial();
  G4String theMatName = amat->GetName();
  if(theMatName == MatName)  found = true ;

  return found ;

}


//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

/*
 * $Log: DSVGenerator.cc,v $
 * Revision 1.12  2015/04/24 12:36:54  pagnes
 * generators in materials updated for DS20k
 *
 * Revision 1.11  2015/04/23 14:04:02  pagnes
 * DS20K geometry added (config 10)
 *
 * Revision 1.10  2015/01/14 16:58:37  dfranco
 * added source vial and correspondent commands by Laura and Bianca; manual updated
 *
 * Revision 1.9  2014/10/22 14:03:47  dfranco
 * add method to generate events in Liqud Argon
 *
 * Revision 1.9  2014/10/22 
 * Added the possibility of generating events in the Liquid Argon volume 

 * Revision 1.8  2014/07/18 13:54:49  dfranco
 * Added a new position generator in the Gas Pocket; calibrated the numnber of photons generated per electron in S2; added a new command to generate multiple particles in the same position
 *
 * Revision 1.7  2014/07/09 13:06:18  pagnes
 * Generators in materials fixed
 *
 * Revision 1.6  2014/07/01 12:27:43  cz2
 * In SetUniformTPC, reset the center to (0,0,-36.5mm) and enlarge both the half-z and radius to 200mm.
 *
 * Revision 1.5  2014/06/12 13:10:06  perassos
 * Update of the method SetTPCCenter() to the new TPC center coordinates
 *
 * Revision 1.4  2014/05/21 10:28:08  dfranco
 * added the possibility to shoot ions
 *
 * Revision 1.3  2014/05/08 11:00:39  pagnes
 * Added generator in PMTs stem
 *
 * Revision 1.2  2014/05/07 12:47:33  dfranco
 * bug fixed
 *
 * Revision 1.1  2014/05/07 12:21:05  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.8  2014/04/22 09:52:28  perassos
 * Added the option to simulate bgd from the G2 cryostats
 *
 * Revision 1.7  2014/04/11 10:20:39  perassos
 * Added generation in materials
 *
 * Revision 1.6  2014/02/13 12:03:26  dfranco
 * added new commands for spatial distributions and fixed the manual
 *
 * Revision 1.5  2013/10/15 16:39:39  dfranco
 * added the possibility to define the energy spectrum reading an ascii file
 *
 * Revision 1.4  2013/05/06 19:54:38  swesterd
 * test
 *
 * Revision 1.3  2013/05/06 10:25:25  dfranco
 * Fixed cylindrical spatial distribution
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
