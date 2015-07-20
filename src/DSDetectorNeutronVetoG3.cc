#include "DSDetectorNeutronVetoG3.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "DSDetectorPMTNeutronVetoG3.hh"
#include "DSMaterial.hh"
#include "DSLogger.hh"
#include "DSStorage.hh"
#include "DSParameters.hh"
#include <iostream>
#include <fstream>
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4VisAttributes.hh"

using namespace std;


DSDetectorNeutronVetoG3::DSDetectorNeutronVetoG3(G4VPhysicalVolume *myMotherVolume) {
  
  fMotherVolume = myMotherVolume;
  
  G4int volumeNumber = 10000;
    
  DSLog(routine) << " Constructing NeutronVetoG3 Geometry" << endlog ;

  G4bool cellTest = false;
  const double myTwoPi = 2*M_PI*rad;


  // Volumes for constructing the top and bottom flanges
  G4double myTopFlange_d = 382.5*cm;
  G4double myBottomFlange_d = 210.0*cm;

  G4double myTopCylinder_h = 60*cm;
  G4double myBottomCylinder_h = 20*cm;

  fSubtrCylinderTop = new G4Tubs("SubtrCylinderTop_Solid", 0, myTopFlange_d/2., myTopCylinder_h/2., 0, myTwoPi );
  fSubtrCylinderBottom = new G4Tubs("SubtrCylinderBottom_Solid", 0, myBottomFlange_d/2., myBottomCylinder_h/2., 0, myTwoPi );



  // Stainless Steel Vessel
  G4double myIDiameter = 736.6*cm;
  G4double myWallThickness = 2.54*cm;
 
  G4double myTopCylinderZshift = 0.5 * pow( myIDiameter*myIDiameter - myTopFlange_d*myTopFlange_d, 0.5 ) + myTopCylinder_h/2. ;
  G4double myBottomCylinderZshift = 0.5 * pow( myIDiameter*myIDiameter - myBottomFlange_d*myBottomFlange_d, 0.5 ) + myBottomCylinder_h/2. ;


  fSolidSteelVessel_0 = new G4Orb("SteelVessel0_Solid", myIDiameter/2. + myWallThickness );
  fSolidSteelVessel_1 = new G4SubtractionSolid( "SteelVessel1_Solid", fSolidSteelVessel_0, fSubtrCylinderTop,    0, G4ThreeVector( 0, 0,  myTopCylinderZshift ));
  fSolidSteelVessel   = new G4SubtractionSolid( "SteelVessel_Solid",  fSolidSteelVessel_1, fSubtrCylinderBottom, 0, G4ThreeVector( 0, 0, -myBottomCylinderZshift ));
  fLogicSteelVessel  = new G4LogicalVolume(fSolidSteelVessel, DSMaterial::Get()->GetStainlessSteel(), "SteelVessel_Logic");
  fPhysicSteelVessel = new G4PVPlacement(0,
				 G4ThreeVector(0,0,0),
				 "SteelVessel",
				 fLogicSteelVessel,
				 fMotherVolume,
				 false,
				 volumeNumber++,
				 DSStorage::Get()->GetCheckOverlap());

  DSLog(routine) << fPhysicSteelVessel->GetName() << " = " << fPhysicSteelVessel->GetCopyNo() << endlog;



  // Top and Bottom Flanges
  G4ThreeVector myTopFlangePos( 0, 0, myTopCylinderZshift - myTopCylinder_h/2. + myWallThickness/2.);
  G4ThreeVector myBottomFlangePos( 0, 0, -myBottomCylinderZshift + myBottomCylinder_h/2. - myWallThickness/2.);


  fSolidTopFlange  = new G4Tubs("TopFlange_Solid", 0, myTopFlange_d/2., myWallThickness/2., 0, myTwoPi);
  fLogicTopFlange  = new G4LogicalVolume( fSolidTopFlange, DSMaterial::Get()->GetStainlessSteel(), "TopFlange_Logic");
  fPhysicTopFlange = new G4PVPlacement( 0, myTopFlangePos, "TopFlange", fLogicTopFlange, fMotherVolume, false, volumeNumber++, DSStorage::Get()->GetCheckOverlap());

  fSolidBottomFlange  = new G4Tubs("BottomFlange_Solid", 0, myBottomFlange_d/2., myWallThickness/2., 0, myTwoPi);
  fLogicBottomFlange  = new G4LogicalVolume( fSolidBottomFlange, DSMaterial::Get()->GetStainlessSteel(), "BottomFlange_Logic");
  fPhysicBottomFlange = new G4PVPlacement( 0, myBottomFlangePos, "BottomFlange", fLogicBottomFlange, fMotherVolume, false, volumeNumber++, DSStorage::Get()->GetCheckOverlap());



  // Boron Scintillator
  fSolidBScintillator_0 = new G4Orb( "BoronScintillator0_Solid", myIDiameter/2. );
  fSolidBScintillator_1 = new G4SubtractionSolid( "BoronScintillator1_Solid", fSolidBScintillator_0, fSubtrCylinderTop, 0, G4ThreeVector( 0, 0, myTopCylinderZshift ));
  fSolidBScintillator   = new G4SubtractionSolid( "BoronScintillator_Solid", fSolidBScintillator_1, fSubtrCylinderBottom, 0, G4ThreeVector( 0, 0, -myBottomCylinderZshift ));
  fLogicBScintillator  = new G4LogicalVolume(fSolidBScintillator, DSMaterial::Get()->GetGdWater(), "BoronScintillator_Logic");
  DSLog(routine) << "Gd Water in the neutron veto" << endlog ;  
  
  
  /*
  if(cellTest) {
      fLogicBScintillator  = new G4LogicalVolume(fSolidBScintillator, DSMaterial::Get()->GetGaseousNitrogen(), "BoronScintillator_Logic");
  } else {
    if(!DSStorage::Get()->GetScintillator()) {
      fLogicBScintillator  = new G4LogicalVolume(fSolidBScintillator, DSMaterial::Get()->GetBoronScintillator(), "BoronScintillator_Logic");
      DSLog(routine) << "Borate Scintillator in the neutron veto" << endlog ;
    } else { 
      fLogicBScintillator  = new G4LogicalVolume(fSolidBScintillator, DSMaterial::Get()->GetGdScintillator(), "BoronScintillator_Logic");
      DSLog(routine) << "Gd Scintillator in the neutron veto" << endlog ;
    }
  }*/
  
  
  fPhysicBScintillator = new G4PVPlacement(0,
					   G4ThreeVector(0,0,0),
					   "BoronScintillatorVolume",
					   fLogicBScintillator,
					   fPhysicSteelVessel,
					   false,
					   volumeNumber++,
					   DSStorage::Get()->GetCheckOverlap());

  DSLog(routine) << fPhysicBScintillator->GetName() << " = " << fPhysicBScintillator->GetCopyNo() << endlog;


  if(cellTest)
    {  
      G4double cellDepth = 0*m;
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
				     G4ThreeVector(1*m,0,cellDepth),
				     "Cell",
				     fLogicCell,
				     fPhysicBScintillator,
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
    }

  //PMT
  new DSDetectorPMTNeutronVetoG3( fPhysicBScintillator );

  /*  
  //Lumirror Sheath outside flat part of cryostat
  G4double sheathGap = 1.*nm;
  G4double cryoSheathMaxZ = DSParameters::Get()->GetCryoSheathZ()[0];
  G4double cryoSheathMinZ = DSParameters::Get()->GetCryoSheathZ()[1];
  G4double cryoSheathR = DSParameters::Get()->GetCryoSheathR()[0]+sheathGap;
  G4double sheathThickness = 1.*nm;
  fSolidCryoSheath = new G4Tubs("CryoSheath_Solid", 
				cryoSheathR, 
				cryoSheathR+sheathThickness,
				(cryoSheathMaxZ - cryoSheathMinZ)/2.,
				0, 2*M_PI);
  fLogicCryoSheath = new G4LogicalVolume(fSolidCryoSheath, DSMaterial::Get()->GetBoronScintillator(), "CryoSheath_Logic");
  fPhysicCryoSheath = new G4PVPlacement(0,
					G4ThreeVector(0,0,(cryoSheathMaxZ+cryoSheathMinZ)/2.),
					"CryoSheath",
					fLogicCryoSheath,
					fPhysicBScintillator,
					false,
					volumeNumber++,
					DSStorage::Get()->GetCheckOverlap());
  DSLog(routine) << fPhysicCryoSheath->GetName() << " = " << fPhysicCryoSheath->GetCopyNo() << endlog;

  fLogicCryoSheath->SetVisAttributes(new G4VisAttributes(G4Colour(0,.5,.5)));
  */



  DefineSurfaces();
}

DSDetectorNeutronVetoG3::~DSDetectorNeutronVetoG3(){
  ; //delete fMessenger;
}
void DSDetectorNeutronVetoG3::DefineSurfaces(){
  fOpUntreatedStainlessSteelSurface = new G4OpticalSurface("OpUntreatedStainlessSteelSurface");
  fOpUntreatedStainlessSteelSurface->SetType(dielectric_metal);
  fOpUntreatedStainlessSteelSurface->SetModel(unified);
  fOpUntreatedStainlessSteelSurface->SetFinish(groundbackpainted);
  fOpUntreatedStainlessSteelSurface->
    SetMaterialPropertiesTable(DSMaterial::Get()->GetUntreatedStainlessSteelMPT());

  fOpElectropolishedStainlessSteelSurface = new G4OpticalSurface("OpElectropolishedStainlessSteelSurface");
  fOpElectropolishedStainlessSteelSurface->SetType(dielectric_metal);
  fOpElectropolishedStainlessSteelSurface->SetModel(unified);
  fOpElectropolishedStainlessSteelSurface->SetFinish(groundbackpainted);
  fOpElectropolishedStainlessSteelSurface->
    SetMaterialPropertiesTable(DSMaterial::Get()->GetElectropolishedStainlessSteelMPT());

  fOpAluminumFoilSurface = new G4OpticalSurface("OpAluminumFoilSurface");
  fOpAluminumFoilSurface->SetType(dielectric_metal);
  fOpAluminumFoilSurface->SetModel(unified);
  fOpAluminumFoilSurface->SetFinish(groundbackpainted);
  fOpAluminumFoilSurface->
    SetMaterialPropertiesTable(DSMaterial::Get()->GetAluminumFoilMPT());

  fOpLumirrorSurface = new G4OpticalSurface("OpLumirrorSurface");
  fOpLumirrorSurface->SetType(dielectric_metal);
  fOpLumirrorSurface->SetModel(unified);
  fOpLumirrorSurface->SetFinish(groundfrontpainted);
  fOpLumirrorSurface->
    SetMaterialPropertiesTable(DSMaterial::Get()->GetLumirrorMPT());
  
  fSteelInnerSurface = new G4LogicalBorderSurface("SSteelInnerSurface",
						  fPhysicBScintillator,
						  fPhysicSteelVessel,
						  fOpLumirrorSurface);

  fSteelInnerSurfaceFlip = new G4LogicalBorderSurface("SSteelInnerSurface",
						      fPhysicSteelVessel,
						      fPhysicBScintillator,
						      fOpLumirrorSurface);
					   
  fSteelOuterSurface = new G4LogicalBorderSurface("SSteelOuterSurface",
						  fPhysicSteelVessel,
						  fMotherVolume,
						  fOpUntreatedStainlessSteelSurface);
  /*
  fCryoSheathSurface = new G4LogicalBorderSurface("CryoSheathSurface",
						  fPhysicBScintillator,
						  fPhysicCryoSheath,
						  fOpLumirrorSurface);
  */

  fTopFlangeSurface = new G4LogicalBorderSurface("TopFlangeSurface",
						 fPhysicBScintillator,
						 fPhysicTopFlange,
						 fOpLumirrorSurface);
  //fOpUntreatedStainlessSteelSurface);
}
