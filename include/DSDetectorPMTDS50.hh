#ifndef DSDetectorPMTDS50_H
#define DSDetectorPMTDS50_H

#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4MaterialPropertiesTable.hh"

class  DSDetectorPMTDS50 {

  public:

    DSDetectorPMTDS50( G4VPhysicalVolume* , G4VPhysicalVolume* );
   ~DSDetectorPMTDS50();

    //G4VPhysicalVolume* GetDetectorComponent()   { return  fPhysicPMTtop; }   // to be changed


  private:
    void                        DefineSurfaces();

    G4int*                      fPMTnum;

    G4VPhysicalVolume*          fMotherVolume;
    G4VPhysicalVolume*          fTeflonCapDiskVolume;
    G4VPhysicalVolume*          fInnerLArVolume;

    G4Tubs*                     fSolidPMTSurrounding;
    G4LogicalVolume*            fLogicPMTSurrounding[19];
  
    G4Tubs*                     fSolidPMTHead;
    G4LogicalVolume*            fLogicPMTHead[19];

    G4Tubs*                     fSolidPMTHeadVac;
    G4LogicalVolume*            fLogicPMTHeadVac[19];

    G4Tubs*                     fSolidPMTWindow;
    G4LogicalVolume*            fLogicPMTWindow[19];

    G4Tubs*                     fSolidPMTJoinVac;
    G4LogicalVolume*            fLogicPMTJoinVac[19];

    G4Tubs*                     fSolidPMTBody;
    G4LogicalVolume*            fLogicPMTBody[19];

    G4Tubs*                     fSolidPMTBodyVac;
    G4LogicalVolume*            fLogicPMTBodyVac[19];

    G4Tubs*                     fSolidPMTTop;
    G4LogicalVolume*            fLogicPMTTop[19];

    G4Tubs*                     fSolidPMTTopVac;
    G4LogicalVolume*            fLogicPMTTopVac[19];

    G4Tubs*                     fSolidPMTLArDisk;
    G4LogicalVolume*            fLogicPMTLArDisk[19];

    G4Tubs*			fSolidPMTStem;
    G4LogicalVolume*		fLogicPMTStem[19];
    G4VPhysicalVolume*  	fPhysicPMTStem[19];
    G4VPhysicalVolume*          fPhysicPMTSurrounding[19];
    G4VPhysicalVolume*          fPhysicPMTHead[19];
    G4VPhysicalVolume*          fPhysicPMTHeadVac[19];
    G4VPhysicalVolume*          fPhysicPMTWindow[19];
    G4VPhysicalVolume*          fPhysicPMTJoinVac[19];
    G4VPhysicalVolume*          fPhysicPMTBody[19];
    G4VPhysicalVolume*          fPhysicPMTBodyVac[19];
    G4VPhysicalVolume*          fPhysicPMTTop[19];
    G4VPhysicalVolume*          fPhysicPMTTopVac[19];
    G4VPhysicalVolume*          fPhysicPMTLArDisk[19];


    // Surface
    G4OpticalSurface*           fOpPMTVacuumSurface;
    G4LogicalBorderSurface*     fPMTVacuumSurface[19];
    G4MaterialPropertiesTable*  fPMTVacuumSurfProp;

    G4OpticalSurface*           fOpPMTLArSurface;
    G4LogicalBorderSurface*     fPMTLArSurface[19];
    G4MaterialPropertiesTable*  fPMTLArSurfProp;

    G4OpticalSurface*           fOpTeflonLArSurface;
    G4LogicalBorderSurface*     fTeflonLArSurface[19];
    G4MaterialPropertiesTable*  fTeflonLArSurfProp;


};

#endif
