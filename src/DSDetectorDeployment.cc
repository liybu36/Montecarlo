#include "DSDetectorDeployment.hh"
#include <iostream>
#include <fstream>
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4VisAttributes.hh"
#include "DSMaterial.hh"
#include "DSDetectorNeutronVeto.hh"
#include "DSDetectorWaterTank.hh"
#include "G4SystemOfUnits.hh"
#include "DSLogger.hh"
#include "DSStorage.hh"
#include "G4UIcommand.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "DSParameters.hh"
#include "G4PhysicalConstants.hh"

using namespace std;

DSDetectorDeployment::DSDetectorDeployment(G4VPhysicalVolume *myMotherVolume)
{
    fMotherVolume = myMotherVolume;
    G4int volumeNumber = 10010;
    DSLog(routine) << "Constructing Organic Pipe Geometry " << endlog;
    G4double pipe_x = 250*cm;
    G4double pipe_y = 0*cm;
    G4double pipe_z = 200*cm;
    G4double piperotation = 0;
    G4double pipelength = 400*cm;
    G4double piperadius = 4*cm;
    G4double pipethickness = 2*cm;
    const G4ThreeVector pipeBodyPosition(pipe_x, pipe_y, pipe_z);
    G4RotationMatrix* pipeRotation = new G4RotationMatrix;
    pipeRotation -> rotateX(piperotation);
    
    fSolidpipeBody = new G4Tubs("pipeBody_Solid", piperadius, piperadius+pipethickness, pipelength/2.,0,twopi*rad);
    fLogicpipe = new G4LogicalVolume(fSolidpipeBody, DSMaterial::Get()->GetStainlessSteel(),"pipe_Logic");
    fPhysicpipe = new G4PVPlacement(pipeRotation, pipeBodyPosition, "pipe",fLogicpipe,fMotherVolume, false, volumeNumber++,DSStorage::Get()->GetCheckOverlap());
    fLogicpipe->SetVisAttributes(new G4VisAttributes(G4Color(0,1,0)));
    
    fSolidpipeScintBody = new G4Tubs("pipeScintBody_Solid", 0, piperadius, pipelength/2. ,0,twopi*rad);
    fLogicpipeScint = new G4LogicalVolume(fSolidpipeScintBody,DSMaterial::Get()->GetBoronScintillator(),"pipeScint_logic");
    fPhysicpipeScint = new G4PVPlacement(0 ,G4ThreeVector(0,0,0), "pipeScint",fLogicpipeScint,fPhysicpipe, false, volumeNumber++, DSStorage::Get()->GetCheckOverlap());
    fLogicpipeScint->SetVisAttributes(new G4VisAttributes(G4Color(1,0,0)));
    
    //  DefineSurfaces();
    
}
/*
void DSDetectorDeployment::DefineSurfaces()
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
 
 fpipeOuterSurface = new G4LogicalBorderSurface("pipeSSteelOuterSurface",
         fPhysicpipe,
	 fMotherVolume,
	 fOpUntreatedStainlessSteelSurface);

 fpipeInnerSurface = new G4LogicalBorderSurface("pipeSSteelInnerSurface",
         fPhysicpipeScint,
	 fPhysicpipe,
	 fOpUntreatedStainlessSteelSurface);
  
}
*/
