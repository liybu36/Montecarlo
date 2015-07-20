#ifndef DSDetectorPMTDSG2_H
#define DSDetectorPMTDSG2_H

#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4MaterialPropertiesTable.hh"

#define NPMTS 279

class  DSDetectorPMTDSG2 {

  public:

    DSDetectorPMTDSG2( G4VPhysicalVolume* , G4VPhysicalVolume* );
   ~DSDetectorPMTDSG2();

    //G4VPhysicalVolume* GetDetectorComponent()   { return  fPhysicPMTtop; }   // to be changed


  private:
    void                        DefineSurfaces();

    G4int*                      fPMTnum;

    G4VPhysicalVolume*          fMotherVolume;
    G4VPhysicalVolume*          fTeflonCapDiskVolume;
    G4VPhysicalVolume*          fInnerLArVolume;

    G4Tubs*                     fSolidPMTSurrounding;
    G4LogicalVolume*            fLogicPMTSurrounding[NPMTS];
  
    G4Tubs*                     fSolidPMTHead;
    G4LogicalVolume*            fLogicPMTHead[NPMTS];

    G4Tubs*                     fSolidPMTHeadVac;
    G4LogicalVolume*            fLogicPMTHeadVac[NPMTS];

    G4Tubs*                     fSolidPMTWindow;
    G4LogicalVolume*            fLogicPMTWindow[NPMTS];

    G4Tubs*                     fSolidPMTJoinVac;
    G4LogicalVolume*            fLogicPMTJoinVac[NPMTS];

    G4Tubs*                     fSolidPMTBody;
    G4LogicalVolume*            fLogicPMTBody[NPMTS];

    G4Tubs*                     fSolidPMTBodyVac;
    G4LogicalVolume*            fLogicPMTBodyVac[NPMTS];

    G4Tubs*                     fSolidPMTTop;
    G4LogicalVolume*            fLogicPMTTop[NPMTS];

    G4Tubs*                     fSolidPMTTopVac;
    G4LogicalVolume*            fLogicPMTTopVac[NPMTS];

    G4Tubs*                     fSolidPMTLArDisk;
    G4LogicalVolume*            fLogicPMTLArDisk[NPMTS];


    G4VPhysicalVolume*          fPhysicPMTSurrounding[NPMTS];
    G4VPhysicalVolume*          fPhysicPMTHead[NPMTS];
    G4VPhysicalVolume*          fPhysicPMTHeadVac[NPMTS];
    G4VPhysicalVolume*          fPhysicPMTWindow[NPMTS];
    G4VPhysicalVolume*          fPhysicPMTJoinVac[NPMTS];
    G4VPhysicalVolume*          fPhysicPMTBody[NPMTS];
    G4VPhysicalVolume*          fPhysicPMTBodyVac[NPMTS];
    G4VPhysicalVolume*          fPhysicPMTTop[NPMTS];
    G4VPhysicalVolume*          fPhysicPMTTopVac[NPMTS];
    G4VPhysicalVolume*          fPhysicPMTLArDisk[NPMTS];


    // Surface
    G4OpticalSurface*           fOpPMTVacuumSurface;
    G4LogicalBorderSurface*     fPMTVacuumSurface[NPMTS];
    G4MaterialPropertiesTable*  fPMTVacuumSurfProp;

    G4OpticalSurface*           fOpPMTLArSurface;
    G4LogicalBorderSurface*     fPMTLArSurface[NPMTS];
    G4MaterialPropertiesTable*  fPMTLArSurfProp;

    G4OpticalSurface*           fOpTeflonLArSurface;
    G4LogicalBorderSurface*     fTeflonLArSurface[NPMTS];
    G4MaterialPropertiesTable*  fTeflonLArSurfProp;


};

#endif
