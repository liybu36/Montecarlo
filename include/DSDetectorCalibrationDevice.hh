#ifndef DSDetectorCalibrationDevice_H
#define DSDetectorCalibrationDevice_H

#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
//#include "G4OpticalSurface.hh"

class DSMaterial;

class DSDetectorCalibrationDevice {

  public:

    DSDetectorCalibrationDevice(G4VPhysicalVolume*);
    ~DSDetectorCalibrationDevice();

    G4VPhysicalVolume* GetDetectorComponent() { return fPhysicCollimator; }


  private:
    DSMaterial*			dsmaterial;

    G4VPhysicalVolume*		fMotherVolume;

    G4Tubs* 			fSolidCollimator;
    G4LogicalVolume*		fLogicCollimator;
    G4VPhysicalVolume*		fPhysicCollimator;

    G4Tubs*			fSolidColEnd;
    G4LogicalVolume*		fLogicColEnd;
    G4VPhysicalVolume*		fPhysicColEndBack;
    G4VPhysicalVolume*		fPhysicColEndFront;

    G4Tubs*		 	fSolidColFill;
    G4LogicalVolume*		fLogicColFill;
    G4VPhysicalVolume*		fPhysicColFill;


};

#endif
