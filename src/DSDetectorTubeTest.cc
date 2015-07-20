#include <iostream>
#include <fstream>
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4VisAttributes.hh"
#include "DSMaterial.hh"
#include "DSDetectorTubeTest.hh"
#include "G4SystemOfUnits.hh"
#include "DSLogger.hh"
#include "DSStorage.hh"
#include "G4UIcommand.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "DSParameters.hh"
#include "G4PhysicalConstants.hh"
using namespace std;

DSDetectorTubeTest::DSDetectorTubeTest(G4VPhysicalVolume *myMotherVolume)
{
  fMotherVolume = myMotherVolume;
  G4int volumeNumber = 30003;
  DSLog(routine) << "Constructing Tube Test Geometry " << endlog;
  G4double source_x = 47.78*cm; // 62.088*cm;  // 67*cm;
  G4double source_y = 0*cm;
  G4double source_z = 11.39*cm;  //14.18*cm;  // 25.43*cm; // -3.65*cm;
  // G4double source_z = -3.65*cm; // -3.65*cm;
  //  G4double source_sphere_radius = 1*cm;
  G4double tubelength = 37*cm;   // 40*cm;
  G4double tubeinnerradius = 11.1*mm; // 1*cm;
  G4double tubethickness = 3.2*mm;
  G4double tubeendcapradius = tubeinnerradius+tubethickness;
  //  G4double tubeendcaplength = 0.1*cm;
  //  G4double tuberotationangle = 0*rad;
  G4double theta =(180.-52)*M_PI/180.;  // 3*M_PI/4; //rotation angle along Y axis 
  //  G4double theta = M_PI/2; //rotation angle along Y axis 
  //  G4double phi = M_PI;  //rotation angle along Z axis
  // const G4ThreeVector SourcePosition(source_x+source_sphere_radius+tubeendcaplength,source_y,source_z);
  G4RotationMatrix* TubeRotation = new G4RotationMatrix;
  TubeRotation -> rotateY(-theta);
  //  TubeRotation -> rotateZ(phi);
  //  const G4ThreeVector TubePosition(-(source_x-tubelength/2.), source_y, source_z); 
  const G4ThreeVector TubePosition(-source_x, source_y, source_z); 

  //outer layer is Tube Body 
  fSolidTubeBody = new G4Tubs("TubeBody_Solid", tubeinnerradius, tubeinnerradius+tubethickness, tubelength/2.,0,twopi*rad);
  fLogicTubeBody = new G4LogicalVolume(fSolidTubeBody, DSMaterial::Get()->GetStainlessSteel(),"TubeBody_Logic");
  fPhysicTubeBody = new G4PVPlacement(TubeRotation, TubePosition, "TubeBody",fLogicTubeBody,fMotherVolume, false, volumeNumber++,DSStorage::Get()->GetCheckOverlap());
  fLogicTubeBody->SetVisAttributes(new G4VisAttributes(G4Color(0,0,0)));
  DSLog(routine) << "Constructing Tube Body: " << volumeNumber <<endlog;
  DSLog(routine) << "TubeInnerRadius= " << fSolidTubeBody->GetInnerRadius() << "   " << "TubeOuterRadius= " <<fSolidTubeBody->GetOuterRadius()<<endlog;
  DSLog(routine) << "TubeLength= " << (fSolidTubeBody->GetZHalfLength())*2. <<endlog;
  /*  //the end caps of the tube
  fSolidTubeEndCap = new G4Tubs("TubeEndCap_Solid", 0, tubeendcapradius, tubeendcaplength/2.,0,twopi*rad);
  fLogicTubeEndCap = new G4LogicalVolume(fSolidTubeEndCap, DSMaterial::Get()->GetStainlessSteel(),"TubeEndCap_Logic");
  fPhysicTubeEndCap = new G4PVPlacement(TubeRotation, G4ThreeVector(source_x+tubeendcaplength/2.,source_y,source_z),"TubeEndCap",fLogicTubeEndCap,fMotherVolume, false,volumeNumber++,DSStorage::Get()->GetCheckOverlap());
  fLogicTubeEndCap->SetVisAttributes(new G4VisAttributes(G4Color(0,0,1)));
  DSLog(routine) << "Constructing End Cap of Tube: " << volumeNumber << endlog;
  */
  /*  //Source Sphere position
  fSolidSphere = new G4Sphere("Source_Sphere",0,source_sphere_radius,0,twopi*rad,0,pi*rad);
  fLogicSphere = new G4LogicalVolume(fSolidSphere, DSMaterial::Get()->GetStainlessSteel(),"Sphere_Logic");
  fPhysicSphere = new G4PVPlacement(0, SourcePosition, "Source-Sphere",fLogicSphere,fMotherVolume, false, 1,DSStorage::Get()->GetCheckOverlap());
  fLogicSphere->SetVisAttributes(new G4VisAttributes(G4Color(0,0,0)));
  */
  //inner layer is filled woth the material( vacuum or Nitrogern)
  fSolidTubeScint = new G4Tubs("TubeScint_Solid", 0, tubeinnerradius,tubelength/2.,0,twopi*rad);
  fLogicTubeScint = new G4LogicalVolume(fSolidTubeScint,DSMaterial::Get()->GetVacuum(),"TubeScint_Logic");
  //  fLogicTubeScint = new G4LogicalVolume(fSolidTubeScint,DSMaterial::Get()->GetBoronScintillator(),"TubeScint_Logic");
  fPhysicTubeScint = new G4PVPlacement(TubeRotation,TubePosition,"TubeScint",fLogicTubeScint,fMotherVolume, false, volumeNumber++, DSStorage::Get()->GetCheckOverlap());
  //  fPhysicTubeScint = new G4PVPlacement(0,0,"TubeScint",fLogicTubeScint,fPhysicTubeBody, false, volumeNumber++, DSStorage::Get()->GetCheckOverlap());
  fLogicTubeScint->SetVisAttributes(new G4VisAttributes(G4Color(1,0,0)));
  DSLog(routine) << "Constructing Scint part(Vacuum or Nitrogen) of Tube: " << volumeNumber << endlog;

  //DefineSurfaces();
}
DSDetectorTubeTest::~DSDetectorTubeTest(){
  ;
}


/*void DSDetectorTubeTest::DefineSurfaces()
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
