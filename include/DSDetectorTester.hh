#ifndef DSDetectorTester_H
#define DSDetectorTester_H

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Para.hh"
#include "G4Trap.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"



class	DSMaterial;

class   DSDetectorTester  {

  public:

    DSDetectorTester(G4VPhysicalVolume*);
   ~DSDetectorTester();

    G4VPhysicalVolume* GetDetectorComponent()  { return  fPhysicSteelTank; }
    

  private:

    void  DefineSurfaces();
    
    
    G4VPhysicalVolume*          fMotherVolume;
    
    G4Tubs*                     fSolidSteelTank;
    G4LogicalVolume*		fLogicSteelTank;
    G4VPhysicalVolume*		fPhysicSteelTank;

    G4Tubs*                     fSolidLAr;
    G4LogicalVolume*		fLogicLAr;
    G4VPhysicalVolume*		fPhysicLAr;

    G4Tubs*                     fSolidTeflon;
    G4LogicalVolume*		fLogicTeflon;
    G4VPhysicalVolume*		fPhysicTeflon;

    G4Tubs*                     fSolidTPB;
    G4LogicalVolume*		fLogicTPB;
    G4VPhysicalVolume*		fPhysicTPB;

};

#endif
/*
 * $Log: DSDetectorTester.hh,v $
 * Revision 1.1  2014/05/07 12:20:51  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.1  2013/05/31 13:02:18  dfranco
 * Added a detector tester, with simpplified geometry (configuration number = 4) to test optical properties of the materials
 *
 *
 */
