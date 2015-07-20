#include <iostream>
#include <fstream>
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4VisAttributes.hh"
#include "DSMaterial.hh"
#include "DSDetectorSourceHolder.hh"
#include "G4SystemOfUnits.hh"
#include "DSLogger.hh"
#include "DSStorage.hh"
#include "G4UIcommand.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "DSParameters.hh"
#include "G4PhysicalConstants.hh"
using namespace std;

DSDetectorSourceHolder::DSDetectorSourceHolder(G4VPhysicalVolume *myMotherVolume)
{
  fMotherVolume = myMotherVolume;

  G4bool rotate = false;
  
  G4int volumeNumber = 30013;
  DSLog(routine) << "Constructing Source Holder Geometry " << endlog;
  G4double source_x = -33.85*cm; 
  G4double source_y = 3.35*cm;
  G4double source_z = -3.65*cm; 
  //Gap between source holder wall and outer cryostat is 2mm
  G4double tubelength = 6.3*cm; 
  G4double tubeinnerradius = 1.35*cm; //2.7cm
  G4double tubethickness = 0.2*cm;
  G4double tubeendcapradius = tubeinnerradius+tubethickness;
  G4double tubeendcaplength = 0.2*cm;
  G4double leadlayerthickness = 0.2*cm;  // 0.1*cm;
  G4double theta =M_PI/2; 
  G4double phi =M_PI/2.;
  G4RotationMatrix* TubeRotation = new G4RotationMatrix;
  if(rotate)
    { TubeRotation -> rotateZ(phi);
      source_x=-103.95*cm;
      source_y=-17*cm;
    }
  TubeRotation -> rotateX(-theta);  

  /*  G4double gap = 30*cm;
  G4double distance = tubelength/2.+tubeendcaplength+gap;
  source_x = source_x-distance;
  source_y = source_y+distance;
  */
  const G4ThreeVector TubePosition(source_x, source_y, source_z); 
  const G4ThreeVector EndCapPosition(source_x, source_y+tubelength/2.+tubeendcaplength/2., source_z); 
  const G4ThreeVector FrontCapPosition(source_x, source_y-tubelength/2.-tubeendcaplength/2., source_z); 

  //outer layer is Tube Body 
  //  fSolidTubeBody = new G4Tubs("TubeBody_Solid", tubeinnerradius, tubeinnerradius+tubethickness, tubelength/2.,0,twopi*rad);
  fSolidTubeBody = new G4Tubs("TubeBody_Solid",0, tubeinnerradius+tubethickness, tubelength/2.+tubethickness,0,twopi*rad);
  fLogicTubeBody = new G4LogicalVolume(fSolidTubeBody, DSMaterial::Get()->GetStainlessSteel(),"TubeBody_Logic");
  fPhysicTubeBody = new G4PVPlacement(TubeRotation, 
				      TubePosition, 
				      "TubeBody",
				      fLogicTubeBody,
				      fMotherVolume, 
				      false, 
				      volumeNumber++,
				      DSStorage::Get()->GetCheckOverlap());

  fLogicTubeBody->SetVisAttributes(new G4VisAttributes(G4Color(0,0,1)));
 
  DSLog(routine) << "Constructing Tube Body: " << volumeNumber <<endlog;
  DSLog(routine) << "TubeInnerRadius= " << fSolidTubeBody->GetInnerRadius() << "   " << "TubeOuterRadius= " <<fSolidTubeBody->GetOuterRadius()<<endlog;
  DSLog(routine) << "TubeLength= " << (fSolidTubeBody->GetZHalfLength())*2. <<endlog;

  //inner layer is filled with metal Lead
  fSolidTubeLead = new G4Tubs("TubeLead_Solid", 0, tubeinnerradius,tubelength/2.,0,twopi*rad);
  fLogicTubeLead = new G4LogicalVolume(fSolidTubeLead,DSMaterial::Get()->GetMetalLead(),"TubeLead_Logic");
  fPhysicTubeLead = new G4PVPlacement(0, //TubeRotation,
				      G4ThreeVector(0,0,0),// TubePosition,
				      "TubeLead",
				      fLogicTubeLead,
				      fPhysicTubeBody,// fMotherVolume, 
				      false,
				      volumeNumber++,
				      DSStorage::Get()->GetCheckOverlap());
  fLogicTubeLead->SetVisAttributes(new G4VisAttributes(G4Color(0,1.0,0)));
  DSLog(routine) << "Constructing Thin layer inside Tube: " << volumeNumber << endlog;
 
  //inner layer is filled woth the material( vacuum or Nitrogern)

  fSolidTubeScint = new G4Tubs("TubeScint_Solid", 0, tubeinnerradius-leadlayerthickness,tubelength/2.-leadlayerthickness,0,twopi*rad);
  //  fSolidTubeScint = new G4Tubs("TubeScint_Solid", 0, tubeinnerradius,tubelength/2.,0,twopi*rad);
  fLogicTubeScint = new G4LogicalVolume(fSolidTubeScint,DSMaterial::Get()->GetVacuum(),"TubeScint_Logic");
  //  fLogicTubeScint = new G4LogicalVolume(fSolidTubeScint,DSMaterial::Get()->GetBoronScintillator(),"TubeScint_Logic");
  fPhysicTubeScint = new G4PVPlacement(0,//TubeRotation,
				       G4ThreeVector(0,0,0),//TubePosition,
				       "TubeScint",
				       fLogicTubeScint,
				       fPhysicTubeLead,//  fMotherVolume, 
				       false, 
				       volumeNumber++,
				       DSStorage::Get()->GetCheckOverlap());
  //  fPhysicTubeScint = new G4PVPlacement(0,0,"TubeScint",fLogicTubeScint,fPhysicTubeBody, false, volumeNumber++, DSStorage::Get()->GetCheckOverlap());
  fLogicTubeScint->SetVisAttributes(new G4VisAttributes(G4Color(1,1,0)));
  DSLog(routine) << "Constructing Scint part(Vacuum or Nitrogen) of Tube: " << volumeNumber << endlog;

  /*
  //the end caps of the tube
  fSolidTubeEndCap = new G4Tubs("TubeEndCap_Solid", 0, tubeendcapradius, tubeendcaplength/2.,0,twopi*rad);
  fLogicTubeEndCap = new G4LogicalVolume(fSolidTubeEndCap, DSMaterial::Get()->GetStainlessSteel(),"TubeEndCap_Logic");

  fPhysicTubeEndCap = new G4PVPlacement(TubeRotation, 
					EndCapPosition,
					"TubeEndCap",
					fLogicTubeEndCap,
					fMotherVolume, 
					false,
					volumeNumber++,
					DSStorage::Get()->GetCheckOverlap());
  
  fLogicTubeEndCap->SetVisAttributes(new G4VisAttributes(G4Color(0,0,1)));
  DSLog(routine) << "Constructing End Cap of Tube: " << volumeNumber << endlog;

  fPhysicTubeFrontCap = new G4PVPlacement(TubeRotation, 
					  FrontCapPosition,
					  "TubeFrontCap",
					  fLogicTubeEndCap,
					  fMotherVolume, 
					  false,
					  volumeNumber++,
					  DSStorage::Get()->GetCheckOverlap());

  DSLog(routine) << "Constructing Front Cap of Tube: " << volumeNumber << endlog;
  */
  /*
  //Source Sphere position
  fSolidSphere = new G4Sphere("Source_Sphere",0,source_sphere_radius,0,twopi*rad,0,pi*rad);
  fLogicSphere = new G4LogicalVolume(fSolidSphere, DSMaterial::Get()->GetStainlessSteel(),"Sphere_Logic");
  fPhysicSphere = new G4PVPlacement(0, SourcePosition, "Source-Sphere",fLogicSphere,fMotherVolume, false, 1,DSStorage::Get()->GetCheckOverlap());
  fLogicSphere->SetVisAttributes(new G4VisAttributes(G4Color(0,0,0)));
  */


  //DefineSurfaces();
}
DSDetectorSourceHolder::~DSDetectorSourceHolder(){
  ;
}


/*void DSDetectorSourceHolder::DefineSurfaces()
{
  fOpUntreatedStainlessSteelSurface = new G4OpticalSurface("OpUntreatedStainlessSteelSurface");
  fOpUntreatedStainlessSteelSurface->SetType(dielectric_metal);
  fOpUntreatedStainlessSteelSurface->SetModel(unified);
  fOpUntreatedStainlessSteelSurface->SetFinish(groundbackpainted);
  fOpUntreatedStainlessSteelSurface->SetMaterialPropertiesTable(DSMaterial::Get()->GetUntreatedStainlessSteelMPT());
  
  fOpElectropolishedStainlessSteelSurface = new G4OpticalSurface("OpElectropolishedStainlessSteelSurface");
  fOpElectropolishedStainlessSteelSurface->SetType(dielectric_metal);
  fOpElectropolishedStainlessSteelSurface->SetModel(unified);
  fOpElectropolishedStainlessSteelSurface->SetFinish(groundbackpainted);
  fOpElectropolishedStainlessSteelSurface->SetMaterialPropertiesTable(DSMaterial::Get()->GetElectropolishedStainlessSteelMPT());
  
  fOpAluminumFoilSurface = new G4OpticalSurface("OpAluminumFoilSurface");
  fOpAluminumFoilSurface->SetType(dielectric_metal);
  fOpAluminumFoilSurface->SetModel(unified);
  fOpAluminumFoilSurface->SetFinish(groundbackpainted);
  fOpAluminumFoilSurface->SetMaterialPropertiesTable(DSMaterial::Get()->GetAluminumFoilMPT());
  
  fOpLumirrorSurface = new G4OpticalSurface("OpLumirrorSurface");
  fOpLumirrorSurface->SetType(dielectric_metal);
  fOpLumirrorSurface->SetModel(unified);
  fOpLumirrorSurface->SetFinish(groundfrontpainted);
  fOpLumirrorSurface->SetMaterialPropertiesTable(DSMaterial::Get()->GetLumirrorMPT());
  
  fSteelOuterSurface = new G4LogicalBorderSurface("SSteelOuterSurface",
                                                  fPhysicSteelVessel,
                                                  fMotherVolume,
                                                  fOpUntreatedStainlessSteelSurface);
						  
  fTopFlangeSurface = new G4LogicalBorderSurface("TopFlangeSurface",
                                                 fPhysicBScintillator,
                                                 fPhysicTopFlange,
                                                 fOpLumirrorSurface);
                                                 
}
*/
