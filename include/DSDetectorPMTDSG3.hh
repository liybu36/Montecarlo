#ifndef DSDetectorPMTDSG3_H
#define DSDetectorPMTDSG3_H

#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4MaterialPropertiesTable.hh"

class  DSDetectorPMTDSG3 {

  public:

    DSDetectorPMTDSG3( G4VPhysicalVolume* , G4VPhysicalVolume* );
   ~DSDetectorPMTDSG3();

    //G4VPhysicalVolume* GetDetectorComponent()   { return  fPhysicPMTtop; }   // to be changed


  private:
    void                        DefineSurfaces();

    G4int*                      fPMTnum;

    G4VPhysicalVolume*          fMotherVolume;
    G4VPhysicalVolume*          fCapDiskVolume;

    G4Tubs*                     fSolidPMTBody;
    G4LogicalVolume*            fLogicPMTBody;
    //G4LogicalVolume*            fLogicPMTBody[379];
    G4VPhysicalVolume*          fPhysicPMTBody[379];

    G4Tubs*                     fSolidPMTBodyVac;
    G4LogicalVolume*            fLogicPMTBodyVac;
    G4VPhysicalVolume*          fPhysicPMTBodyVac;
    //G4LogicalVolume*            fLogicPMTBodyVac[379];
    //G4VPhysicalVolume*          fPhysicPMTBodyVac[379];

    G4Sphere*                   fSolidPMTHead;
    //G4LogicalVolume*            fLogicPMTHead[379];
    G4LogicalVolume*            fLogicPMTHead;
    G4VPhysicalVolume*          fPhysicPMTHead[379];

    G4Sphere*                   fSolidPMTHeadVac;
    G4LogicalVolume*            fLogicPMTHeadVac;
    G4VPhysicalVolume*          fPhysicPMTHeadVac;
    //G4LogicalVolume*            fLogicPMTHeadVac[379];
    //G4VPhysicalVolume*          fPhysicPMTHeadVac[379];




    // Surface
    G4OpticalSurface*           fOpPMTVacuumSurface;
    //G4LogicalBorderSurface*     fPMTVacuumSurface[379];
    G4LogicalBorderSurface*     fPMTVacuumSurface;
    G4MaterialPropertiesTable*  fPMTVacuumSurfProp;

    G4OpticalSurface*           fOpPMTLArSurface;
    G4LogicalBorderSurface*     fPMTLArSurface[379];
    G4MaterialPropertiesTable*  fPMTLArSurfProp;

    G4OpticalSurface*           fOpTeflonLArSurface;
    //G4LogicalBorderSurface*     fTeflonLArSurface[379];
    G4LogicalBorderSurface*     fTeflonLArSurface;
    G4MaterialPropertiesTable*  fTeflonLArSurfProp;


};

#endif
