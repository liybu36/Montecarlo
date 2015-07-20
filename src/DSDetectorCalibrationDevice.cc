#include "DSDetectorCalibrationDevice.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
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

DSDetectorCalibrationDevice::DSDetectorCalibrationDevice(G4VPhysicalVolume *myMotherVolume) {

  fMotherVolume = myMotherVolume;
  
  G4int volumeNumber = 20000;

  DSLog(routine) << " Constructing CalibrationDevice Geometry" << endlog;

  //a cylindrical hollow tube filled with vacuum (in MC) or N2 (in the detector)
  G4double colInnerRadius = 1.0*cm;
  G4double wallThickness = 1.0*mm;
  G4double colOuterRadius = colInnerRadius+wallThickness;
  G4double colHalfLength = 20.*cm;
  G4double colEndCapHalfLength = 0.5*mm;

  G4double GapSourceCollimator = 0.5*cm; //Source outside the collimator, to avoid that energy is deposited in the collimator (big effect for electrons).
  
  /*
  //configuration used by Princeton guys, especially Hao
  G4double source_x=63.0*cm; //end of collimator at 45.1*cm 
  G4double source_y=0*cm;
  G4double source_z=-3.65*cm;
  G4double colRotationAngle=0; //rotation angle in rho-z plane
  */

  G4double colRotationAngle=30*degree; //rotation angle in rho-z plane
  G4RotationMatrix* colRot=new G4RotationMatrix;
  colRot->rotateY(pi/2*rad+colRotationAngle); //rotateY is a rotation around axis Y, the pi/2*rad are from the rotation out of the z-axis


  //32.2 cm: outer cryostat wall
  //3 cm: gap cryostat & collimator
  G4double colCenter_x=32.2*cm+cos(colRotationAngle)*(3*cm+colHalfLength);
  G4double colCenter_y=0*cm;
  G4double colCenter_z=32.2*cm*tan(colRotationAngle)+sin(colRotationAngle)*(3*cm+colHalfLength);
  G4ThreeVector colCenter = G4ThreeVector(colCenter_x, colCenter_y, colCenter_z);
  
  G4ThreeVector colEndFrontCenter = G4ThreeVector(colCenter.x()-(colHalfLength+colEndCapHalfLength)*cos(colRotationAngle),colCenter.y(),colCenter.z()-(colHalfLength+colEndCapHalfLength)*sin(colRotationAngle));
  G4ThreeVector colEndBackCenter = G4ThreeVector(colCenter.x()+(colHalfLength+colEndCapHalfLength)*cos(colRotationAngle),colCenter.y(),colCenter.z()+(colHalfLength+colEndCapHalfLength)*sin(colRotationAngle));

  //this is a sanity check to be compared with the position in the macro. (It is not used further)
  G4double source_x=colCenter.x()+10*cm*cos(colRotationAngle); //end of collimator at 45.1*cm 
  G4double source_y=0*cm;
  G4double source_z=colCenter.z()+10*cm*sin(colRotationAngle);

   
   //ignoring thickness of the collimator for now:
   //G4ThreeVector colCenter = G4ThreeVector(source_x-(colHalfLength+colEndCapHalfLength*2+GapSourceCollimator)*cos(colRotationAngle*degree), source_y, source_z-(colHalfLength+colEndCapHalfLength*2+GapSourceCollimator)*sin(colRotationAngle*degree));


   DSLog(routine) << "collimator dimensions: inner radius: " << colInnerRadius/cm << ", length: " << colHalfLength*2/cm << ", wall thickness: " << wallThickness/cm << ") cm" << endlog;
   DSLog(routine) << "collimator center: (" << colCenter.x()/cm << ", " << colCenter.y()/cm << ", " << colCenter.z()/cm << ") cm" << endlog;
   DSLog(routine) << "source position (in " << __FILE__ << "): (" << source_x/cm << ", " << source_y/cm << ", " << source_z/cm << ") cm" << endlog;
   DSLog(routine) << "rotation angle: (in " << __FILE__ << "): " << colRotationAngle/degree << " deg" << endlog;

  //main collimator volume
  fSolidCollimator = new G4Tubs("Collimator_Solid", colInnerRadius, colOuterRadius, colHalfLength, 0, twopi*rad);
  fLogicCollimator = new G4LogicalVolume(fSolidCollimator, DSMaterial::Get()->GetStainlessSteel(), "Collimator_Logic");
  fPhysicCollimator = new G4PVPlacement(colRot,
			                colCenter,
				        "CollimatorVolume",
				        fLogicCollimator,
			                fMotherVolume,
				        false,
				        volumeNumber++,
				        DSStorage::Get()->GetCheckOverlap());

  DSLog(routine) << fPhysicCollimator->GetName() << " = " << fPhysicCollimator->GetCopyNo() << endlog;

  //collimator end caps
  fSolidColEnd = new G4Tubs("CollimatorEndCap_Solid", 0, colOuterRadius, colEndCapHalfLength, 0, twopi*rad);
  fLogicColEnd = new G4LogicalVolume(fSolidColEnd, DSMaterial::Get()->GetStainlessSteel(), "CollimatorEndCap_Logic");

  //rear end cap
  fPhysicColEndBack = new G4PVPlacement(colRot,
					colEndBackCenter,
				        "CollimatorBackEndCapVolume",
				        fLogicColEnd,
				        fMotherVolume,
			                false,
				        volumeNumber++,
				        DSStorage::Get()->GetCheckOverlap());

  DSLog(routine) << fPhysicColEndBack->GetName() << " = " << fPhysicColEndBack->GetCopyNo() << endlog;

//front end cap
fPhysicColEndFront = new G4PVPlacement(colRot,
				       colEndFrontCenter,
				       "CollimatorFrontEndCapVolume",
				       fLogicColEnd,
				       fMotherVolume,
				       false,
				       volumeNumber++,
				       DSStorage::Get()->GetCheckOverlap());

  DSLog(routine) << fPhysicColEndFront->GetName() << " = " << fPhysicColEndFront->GetCopyNo() << endlog;

//inside of collimator
fSolidColFill = new G4Tubs("CollimatorFill_Solid", 0, colInnerRadius, colHalfLength, 0, twopi*rad);
//fLogicColFill = new G4LogicalVolume(fSolidColFill, DSMaterial::Get()->GetVacuum(), "CollimatorFill_Logic");
fLogicColFill = new G4LogicalVolume(fSolidColFill, DSMaterial::Get()->GetGaseousNitrogen(), "CollimatorFill_Logic");
fPhysicColFill = new G4PVPlacement(colRot,
				   colCenter,
				   "CollimatorFillVolume",
			   	   fLogicColFill,
				   fMotherVolume,
				   false,
				   volumeNumber++,
				   DSStorage::Get()->GetCheckOverlap());

  DSLog(routine) << fPhysicColFill->GetName() << " = " << fPhysicColFill->GetCopyNo() << endlog;
}

