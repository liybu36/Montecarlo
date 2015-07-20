#include "DSScintCell.hh"
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
#include "G4VisAttributes.hh"

// This is for testing the 3" cells Shawn is working with in Princeton

using namespace std;

DSScintCell::DSScintCell(G4VPhysicalVolume *myMotherVolume) {
  
  fMotherVolume = myMotherVolume;  
  
  DSLog(routine) << " Constructing Tester Geometry" << endlog ;
  
  // Numbers for cell #VIII
  G4double myCellRadius = (7.62/2.)*cm;
  G4double myCellHeight = 7.62*cm;
  G4double myStemLength = 6.8*cm;
  G4double myStemRadius = 0.5*cm;
  G4double myAirGap = 3.*mm;
  G4double myCellThickness = 0.5*mm;

  // Cell
  fSolidCellBody = new G4Tubs("CellBody_Solid", 0, myCellRadius, myCellHeight/2.,0,twopi*rad);
  fSolidCellStem = new G4Tubs("CellStem_Solid", 0, myStemRadius, myStemLength/2.,0,twopi*rad);

  fSolidCell = new G4UnionSolid("Cell_Solid", fSolidCellBody, fSolidCellStem, 0, G4ThreeVector(0,0,-myCellHeight/2.-myStemLength/2.));
  fLogicCell = new G4LogicalVolume(fSolidCell, DSMaterial::Get()->GetFusedSilica(), "Cell_Logic");
  fPhysicCell= new G4PVPlacement(0,
				 G4ThreeVector(0,0,0),
				 "Cell",
				 fLogicCell,
				 fMotherVolume,
				 false,
				 1,
				 DSStorage::Get()->GetCheckOverlap());
  
  // Boron Scintillator
  fSolidScintBody = new G4Tubs("ScintBody_Solid", 0, myCellRadius-myCellThickness, myCellHeight/2.-myCellThickness-myAirGap/2.,0,twopi*rad);
  fSolidScintStem = new G4Tubs("ScintStem_Solid", 0, myStemRadius-myCellThickness, myStemLength/2.+myCellThickness/2.,0,twopi*rad);

  fSolidScint = new G4UnionSolid("Scint_Solid", fSolidScintBody, fSolidScintStem, 0, G4ThreeVector(0,0,-myCellHeight/2.-myStemLength/2.+myCellThickness+myAirGap/2.));
  fLogicScint = new G4LogicalVolume(fSolidScint, DSMaterial::Get()->GetBoronScintillator(), "Scint_Logic");
  fPhysicScint= new G4PVPlacement(0,
				 G4ThreeVector(0,0,-myAirGap/2.),
				 "Scint",
				 fLogicScint,
				 fPhysicCell,
				 false,
				 1,
				 DSStorage::Get()->GetCheckOverlap());
  fLogicScint->SetVisAttributes(new G4VisAttributes(G4Colour(0,0,1)));
  
  // Nitrogen bubble
  fSolidBubble = new G4Tubs("Bubble_Solid", 0, myCellRadius-myCellThickness, myAirGap/2., 0, twopi*rad);
  fLogicBubble = new G4LogicalVolume(fSolidBubble, DSMaterial::Get()->GetGaseousNitrogen(), "Bubble_Logic");
  fPhysicBubble= new G4PVPlacement(0,
				   G4ThreeVector(0,0,myCellHeight/2.-myAirGap/2.-myCellThickness),
				   "Bubble",
				   fLogicBubble,
				   fPhysicCell,
				   false,
				   1,
				   DSStorage::Get()->GetCheckOverlap());
  fLogicBubble->SetVisAttributes(new G4VisAttributes(G4Colour(0,1,0)));
  
  // PMT
  G4double pmtHeight = 12.*cm;
  fSolidPMT = new G4Tubs("PMT_Solid", 0, myCellRadius, pmtHeight/2., 0, twopi*rad);
  fLogicPMT = new G4LogicalVolume(fSolidPMT, DSMaterial::Get()->GetBialkali(), "VPMTWindow_0");
  fPhysicPMT= new G4PVPlacement(0,
				//				G4ThreeVector(0,0,myCellHeight/2.-myCellThickness/2.),
				G4ThreeVector(0,0,myCellHeight/2.+pmtHeight/2.),
				"VPMT_0",
				fLogicPMT,
				fMotherVolume,
				false,
				1,
				DSStorage::Get()->GetCheckOverlap());
  fLogicPMT->SetVisAttributes(new G4VisAttributes(G4Colour(1,0,0)));
  DSStorage::Get()->SetVetoPMTMaterialIndex(fLogicPMT->GetMaterial()->GetIndex());
  //  DSStorage::Get()->SetPMTMaterialIndex(fLogicPMT->GetMaterial()->GetIndex());

  DefineSurfaces();
  
}
DSScintCell::~DSScintCell(){
  ;
}
void DSScintCell::DefineSurfaces(){
  fOpLumirrorSurface = new G4OpticalSurface("OpLumirrorSurface");
  fOpLumirrorSurface->SetType(dielectric_metal);
  fOpLumirrorSurface->SetModel(unified);
  fOpLumirrorSurface->SetFinish(groundfrontpainted);
  fOpLumirrorSurface->
    SetMaterialPropertiesTable(DSMaterial::Get()->GetLumirrorMPT());

  fLumirrorSurface = new G4LogicalBorderSurface("LumirrorSurface",
					       fPhysicCell,
					       fMotherVolume,
					       fOpLumirrorSurface);

  fOpPMTSurface = new G4OpticalSurface("OpPMTSurface");
  fOpPMTSurface->SetType(dielectric_metal);
  fOpPMTSurface->SetModel(glisur);
  fOpPMTSurface->SetFinish(polished);
  fOpPMTSurface->
    SetMaterialPropertiesTable(DSMaterial::Get()->GetVPhotocathodeMPT());
  
  fPMTSurface = new G4LogicalBorderSurface("PMTSurface", 
					   fPhysicCell, 
					   fPhysicPMT,
					   //				   fOpLumirrorSurface);
					   fOpPMTSurface);
  
  /*
  G4MaterialPropertiesTable *fPMTSurfProp = new G4MaterialPropertiesTable();
  fPMTSurfProp->AddConstProperty("REFLECTIVITY", 0.2);
  fPMTSurfProp->AddConstProperty("EFFICIENCY",   0.0);
  fPMTSurface = new G4LogicalBorderSurface("PMTSurface", fPhysicCell, fPhysicPMT, fOpPMTSurface );  

  fOpPMTSurface = new G4OpticalSurface("OpPMTSurface");
  fOpPMTSurface->SetType( dielectric_metal );
  fOpPMTSurface->SetModel( glisur );
  fOpPMTSurface->SetFinish( polished );  
  fOpPMTSurface->SetMaterialPropertiesTable( fPMTSurfProp );
  */
}
