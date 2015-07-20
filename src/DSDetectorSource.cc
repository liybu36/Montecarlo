#include <iostream>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"
#include "DSDetectorSource.hh"
#include "DSIO.hh"
#include "DSLogger.hh"
#include "DSMaterial.hh"
#include "DSStorage.hh"
#include "DSParameters.hh"
#include "G4VPhysicalVolume.hh"
#include "G4DisplacedSolid.hh"

using namespace std;

////////////////////////////////////////////////
//////        Detector Description      ////////
////////////////////////////////////////////////
/*

From the center (0,0,0) to the top
 
LiquidArgon     ActiveLAr          27.537 cm 
GaseousArgon    GasPocket          28.474 cm 
TPB             TPB                28.484 cm 
FusedSilica     BellTop            28.802 cm 
LiquidArgon     InnerLiqArgon      29.304 cm 
LiquidArgon     PMTLArDisk_0       29.404 cm 
Bialkali        TPMT_0             29.7627 cm

Back of the PMT:

Vacuum          PMTHeadVac_0       32.853 cm    
Vacuum          PMTJoinVac_0       32.908 cm    
Vacuum          PMTBodyVac_0       36.365 cm 
Vacuum          PMTTopVac_0        41.592 cm 
Kovar           PMTTop_0           41.647 cm 
LiquidArgon     PMTAssemblyTop     43.528 cm 
LiquidArgon     OuterLiquidArgon   44.4 cm 
GaseousArgon    GaseousArgon       62.141 cm 
StainlessSteel  InnerCryostat      62.6 cm  
Vacuum          VacuumCryostat     63.4 cm  
GaseousArgon    TrunkAr            123.4 cm  
Air             World              126.301 cm
StainlessSteel  Trunk6             126.584 cm
GaseousArgon    TrunkAr            131.406 cm
StainlessSteel  Trunk6             131.689 cm
Air             World              2000 cm  


From the center (0,0,0) to the left

LiquidArgon     ActiveLAr          17.77 cm  
TPB             TPB                17.78 cm  
Teflon          Reflector          20.32 cm  
MetalCopper     FieldRings         20.7367 cm  
LiquidArgon     InnerLiqArgon      21.59 cm  
Teflon          TeflonSupport      23.495 cm  
LiquidArgon     OuterLiquidArgon   25.197 cm  
StainlessSteel  InnerCryostat      25.65 cm  
Vacuum          VacuumCryostat     31.6983 cm  
StainlessSteel  OuterCryostat      32.1 cm  


*/
////////////////////////////////////////////////



DSDetectorSource::DSDetectorSource(G4VPhysicalVolume *myMotherVolume) {
  
  fMotherVolume = myMotherVolume;
  const double myTwoPi = 2*M_PI*rad;
  bool   myCheckOverlap   = DSStorage::Get()->GetCheckOverlap();
  
  fPosition = DSStorage::Get()->GetSourcePosition();


  DSLog(routine) << " Constructing Source Geometry" << endlog ;

  G4RotationMatrix* rotZ90 = new G4RotationMatrix;
  rotZ90->rotateZ( M_PI/2.*rad );

  G4Tubs *mySolidVial            = new G4Tubs ("Vial_Solid", 0, 14.5*mm, 5*mm, 0, myTwoPi);
  G4LogicalVolume *myLogicVial   = new G4LogicalVolume(mySolidVial, DSMaterial::Get()->GetStainlessSteel(), "Vial_Logic");
  G4VPhysicalVolume *myPhysicVial = new G4PVPlacement( rotZ90, fPosition, "Vial", myLogicVial, fMotherVolume, true, 0, myCheckOverlap );


  G4Tubs *mySolidSource             = new G4Tubs ("Source_Solid", 0, 12.5*mm, 3*mm, 0, myTwoPi);
  G4LogicalVolume *myLogicSource    = new G4LogicalVolume(mySolidSource, DSMaterial::Get()->GetTeflon(), "Source_Logic");
  new G4PVPlacement( 0, G4ThreeVector(0,0,0), "Source", myLogicSource, myPhysicVial, true, 1, myCheckOverlap );
 

  DefineSurfaces(); 

}

DSDetectorSource::~DSDetectorSource(){
  ; //delete fMessenger;
}

void  DSDetectorSource::DefineSurfaces() {
  ;
}  



/*
 * $Log: DSDetectorSource.cc,v $
 * Revision 1.1  2014/11/06 17:39:45  dfranco
 * added source and fixed a bug: removed gaseous nitrogen from the neutron veto
 *
 *
 */
