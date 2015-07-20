#include "DSDetectorTester.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "DSMaterial.hh"
#include "DSLogger.hh"
#include "DSStorage.hh"
#include <iostream>
#include <fstream>
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4UIcommand.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"

using namespace std;

DSDetectorTester::DSDetectorTester(G4VPhysicalVolume *myMotherVolume) {
  
  fMotherVolume = myMotherVolume;  
  
  DSLog(routine) << " Constructing Tester Geometry" << endlog ;
  
  G4double mySteelRadius    = 20.0*cm;
  G4double mySteelHeight    = 40.0*cm/2.;
  G4double mySteelThickness = 1.0*cm;
  G4double myLArRadius      = mySteelRadius   - mySteelThickness ;
  G4double myLArHeight      = mySteelHeight   - mySteelThickness ;
  G4double myTeflonRadius   = mySteelRadius   - mySteelThickness;
  G4double myTeflonHeight   = 1.0*cm/2.;
  G4double myTPBRadius      = mySteelRadius   - mySteelThickness ;
  G4double myTPBHeight      = 3.28*mm/2.;
  

  //BlackHole Tank
  fSolidSteelTank  = new G4Tubs("BlackHole_Solid", 0,mySteelRadius,mySteelHeight,0,twopi*rad);
  fLogicSteelTank  = new G4LogicalVolume(fSolidSteelTank, DSMaterial::Get()->GetBlackHole(), "BlackHole_Logic");
  fPhysicSteelTank = new G4PVPlacement(0,
				 G4ThreeVector(0,0,0),
				 "BlackHole",
				 fLogicSteelTank,
				 fMotherVolume,
				 false,
				 1,
				 DSStorage::Get()->GetCheckOverlap());

  //LAr
  fSolidLAr  = new G4Tubs("LAr_Solid", 0,myLArRadius,myLArHeight,0,twopi*rad);
  fLogicLAr  = new G4LogicalVolume(fSolidLAr, DSMaterial::Get()->GetLiquidArgon(), "LAr_Logic");
  fPhysicLAr = new G4PVPlacement(0,
				 G4ThreeVector(0,0,0),
				 "LAr",
				 fLogicLAr,
				 fPhysicSteelTank,
				 false,
				 2,
				 DSStorage::Get()->GetCheckOverlap());

  //Teflon volume
  fSolidTeflon  = new G4Tubs("Teflon_Solid", 0,myTeflonRadius,myTeflonHeight,0,twopi*rad);
  fLogicTeflon  = new G4LogicalVolume(fSolidTeflon, DSMaterial::Get()->GetTeflon(), "Teflon_Logic");
  fPhysicTeflon = new G4PVPlacement(0,
				 G4ThreeVector(0,0, myLArHeight - myTeflonHeight),
				 "Teflon",
				 fLogicTeflon,
				 fPhysicLAr,
				 false,
				 3,
				 DSStorage::Get()->GetCheckOverlap());

  //TPB volume
  fSolidTPB  = new G4Tubs("TPB_Solid", 0,myTPBRadius,myTPBHeight,0,twopi*rad);
  fLogicTPB  = new G4LogicalVolume(fSolidTPB, DSMaterial::Get()->GetTPB(), "TPB_Logic");
  fPhysicTPB = new G4PVPlacement(0,
				 G4ThreeVector(0,0, myLArHeight - 2*myTeflonHeight - myTPBHeight),
				 "TPB",
				 fLogicTPB,
				 fPhysicLAr,
				 false,
				 4,
				 DSStorage::Get()->GetCheckOverlap());

  DefineSurfaces();
  
}
DSDetectorTester::~DSDetectorTester(){
  ;
}
void DSDetectorTester::DefineSurfaces(){
  
  
  // BlackHole - LAr
  G4OpticalSurface *fOpSteelLArSurface     = new G4OpticalSurface("OpSteelLArSurface");
  new G4LogicalBorderSurface("OpSteelLArBorderSurface", fPhysicLAr, fPhysicSteelTank, fOpSteelLArSurface );
  fOpSteelLArSurface->SetType( dielectric_dielectric );
  fOpSteelLArSurface->SetModel( glisur );
  fOpSteelLArSurface->SetFinish( polished );               
  G4MaterialPropertiesTable *fMTBSteelLAr = new G4MaterialPropertiesTable();
  fMTBSteelLAr->AddConstProperty("REFLECTIVITY", 0.0);                     
  fMTBSteelLAr->AddConstProperty("EFFICIENCY", 0.0);                     
  fOpSteelLArSurface->SetMaterialPropertiesTable( fMTBSteelLAr );


  // BlackHole - TPB
  G4OpticalSurface *fOpSteelTPBSurface     = new G4OpticalSurface("OpSteelTPBSurface");
  new G4LogicalBorderSurface("OpSteelTPBBorderSurface", fPhysicTPB, fPhysicSteelTank, fOpSteelTPBSurface );
  fOpSteelTPBSurface->SetType( dielectric_dielectric );
  fOpSteelTPBSurface->SetModel( glisur );
  fOpSteelTPBSurface->SetFinish( ground );               
  G4MaterialPropertiesTable *fMTBSteelTPB = new G4MaterialPropertiesTable();
  fMTBSteelTPB->AddConstProperty("REFLECTIVITY", 0.0);                     
  fMTBSteelTPB->AddConstProperty("EFFICIENCY", 0.0);                     
  fOpSteelTPBSurface->SetMaterialPropertiesTable( fMTBSteelTPB );

  // BlackHole - Teflon
  G4OpticalSurface *fOpSteelTeflonSurface     = new G4OpticalSurface("OpSteelTeflonSurface");
  new G4LogicalBorderSurface("OpSteelTeflonBorderSurface", fPhysicSteelTank, fPhysicTeflon, fOpSteelTeflonSurface );
  fOpSteelTeflonSurface->SetType( dielectric_metal );
  fOpSteelTeflonSurface->SetFinish( ground );    
  fOpSteelTeflonSurface->SetModel( glisur );
  G4MaterialPropertiesTable *fMTBSteelTeflon = new G4MaterialPropertiesTable();
  fMTBSteelTeflon->AddConstProperty("REFLECTIVITY", 1.0);      
  fOpSteelTeflonSurface->SetMaterialPropertiesTable( fMTBSteelTeflon );

  // LAr - Teflon
  //G4OpticalSurface *fOpLArTeflonSurface     = new G4OpticalSurface("OpLArTeflonSurface");
  //new G4LogicalBorderSurface("OpLArTeflonBorderSurface", fPhysicTeflon, fPhysicLAr, fOpLArTeflonSurface );
  //fOpLArTeflonSurface->SetType( dielectric_metal );
  //fOpLArTeflonSurface->SetModel( unified );
  //fOpLArTeflonSurface->SetFinish( ground );                           
  //fOpLArTeflonSurface->SetSigmaAlpha(0.3);                          
  //G4MaterialPropertiesTable *fMTBLArTeflon = new G4MaterialPropertiesTable();
  //fMTBLArTeflon->AddConstProperty("EFFICIENCY", 0.0);                     
  //fMTBLArTeflon->AddConstProperty("REFLECTIVITY", 0.1);                     
  //fMTBLArTeflon->AddConstProperty("SPECULARLOBECONSTANT",  1.0);                     
  //fMTBLArTeflon->AddConstProperty("BACKSCATTERCONSTANT",   1.0);                     
  //fMTBLArTeflon->AddConstProperty("SPECULARSPIKECONSTANT", 1.0);                     
  //OpLArTeflonSurface->SetMaterialPropertiesTable( fMTBLArTeflon );

  // LAr - TPB
  G4OpticalSurface *fOpLArTPBSurface     = new G4OpticalSurface("OpLArTPBSurface");
  new G4LogicalBorderSurface("OpLArTPBBorderSurface", fPhysicTPB, fPhysicLAr, fOpLArTPBSurface );
  fOpLArTPBSurface->SetType( dielectric_dielectric );
  fOpLArTPBSurface->SetModel( unified );
  fOpLArTPBSurface->SetFinish( polished );                           
  fOpLArTPBSurface->SetSigmaAlpha(0.0);                          
  G4MaterialPropertiesTable *fMTBLArTPB = new G4MaterialPropertiesTable();
  //fMTBLArTPB->AddConstProperty("EFFICIENCY", 0.0);                     
  //fMTBLArTPB->AddConstProperty("REFLECTIVITY", 0.1);                     
  //fMTBLArTPB->AddConstProperty("SPECULARLOBECONSTANT",  1.0);                     
  //fMTBLArTPB->AddConstProperty("BACKSCATTERCONSTANT",   1.0);                     
  //fMTBLArTPB->AddConstProperty("SPECULARSPIKECONSTANT", 1.0);                     
  fOpLArTPBSurface->SetMaterialPropertiesTable( fMTBLArTPB );

  // Teflon - TPB
  G4OpticalSurface *fOpTeflonTPBSurface     = new G4OpticalSurface("OpTeflonTPBSurface");
  new G4LogicalBorderSurface("OpTeflonTPBBorderSurface",fPhysicTPB, fPhysicTeflon,  fOpTeflonTPBSurface);
  fOpTeflonTPBSurface->SetType( dielectric_metal );
  fOpTeflonTPBSurface->SetModel(unified);
  fOpTeflonTPBSurface->SetFinish(groundfrontpainted); 
  fOpTeflonTPBSurface->SetSigmaAlpha(0.1);                          
  G4MaterialPropertiesTable *fMTBTeflonTPB = new G4MaterialPropertiesTable();
  G4double TeflonTPBENE[2] = {0.1*eV, 20.0*eV};
  G4double TeflonTPBREF[2] = {0.90, 0.90};
  G4double TeflonTPBEFF[2] = {0.00, 0.00};
  fMTBTeflonTPB->AddProperty("REFLECTIVITY", TeflonTPBENE, TeflonTPBREF, 2);			 
  fMTBTeflonTPB->AddProperty("EFFICIENCY",   TeflonTPBENE, TeflonTPBEFF, 2);		       
  
  //fMTBTeflonTPB->AddConstProperty("REFLECTIVITY", 0.7);                     
  //fMTBTeflonTPB->AddConstProperty("EFFICIENCY", 0.0);                     
  //fMTBTeflonTPB->AddConstProperty("SPECULARLOBECONSTANT", 1.0);                     
  //fMTBTeflonTPB->AddConstProperty("BACKSCATTERCONSTANT",  1.0);                     
  //fMTBTeflonTPB->AddConstProperty("SPECULARSPIKECONSTANT", 1.0);                     
  fOpTeflonTPBSurface->SetMaterialPropertiesTable( fMTBTeflonTPB );



}




/*
 * $Log: DSDetectorTester.cc,v $
 * Revision 1.1  2014/05/07 12:21:02  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.5  2013/06/10 13:09:19  dfranco
 * improved optical effects on surfaces in the test detector
 *
 * Revision 1.4  2013/06/05 16:28:30  dfranco
 * Improved optics in the test detector
 *
 * Revision 1.3  2013/06/04 16:56:41  dfranco
 * Added fake BlackHole material to abosrb photons
 *
 * Revision 1.2  2013/05/31 13:28:26  dfranco
 * Added WLSEFFICIENCY to TPB
 *
 * Revision 1.1  2013/05/31 13:02:15  dfranco
 * Added a detector tester, with simpplified geometry (configuration number = 4) to test optical properties of the materials
 *
 *
 */
