#ifndef DSDetectorWaterTank_H
#define DSDetectorWaterTank_H

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Para.hh"
#include "G4Trap.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"



class	DSMaterial;

class   DSDetectorWaterTank  {

  public:

    DSDetectorWaterTank(G4VPhysicalVolume*);
   ~DSDetectorWaterTank();

    G4VPhysicalVolume* GetDetectorComponent()  { return  fPhysicWaterVolume;}
    

  private:

    DSMaterial*                 dsmaterial;
    G4Material*			fWorldMaterial;

    G4VPhysicalVolume*          fMotherVolume;
    
    G4Tubs*                     fSolidSteelTank;
    G4LogicalVolume*		fLogicSteelTank;
    G4VPhysicalVolume*		fPhysicSteelTank;

    G4Tubs*                     fSolidWaterVolume;
    G4LogicalVolume*		fLogicWaterVolume;
    G4VPhysicalVolume*		fPhysicWaterVolume;

    G4Tubs*                     fSolidNVleg;
    G4LogicalVolume*		fLogicNVleg;
    G4VPhysicalVolume*		fPhysicNVleg;

    G4Box*                      fSolidNVlegUp;
    G4LogicalVolume*		fLogicNVlegUp;
    G4VPhysicalVolume*		fPhysicNVlegUp;

    G4Box*                      fSolidNVbarB;
    G4LogicalVolume*		fLogicNVbarB;
    G4VPhysicalVolume*		fPhysicNVbarB;

    G4Box*                      fSolidNVbarH;
    G4LogicalVolume*		fLogicNVbarH;
    G4VPhysicalVolume*		fPhysicNVbarH;

    G4Para*                     fSolidNVbarDR;
    G4LogicalVolume*		fLogicNVbarDR;
    G4VPhysicalVolume*		fPhysicNVbarDR;

    G4Para*                     fSolidNVbarDL;
    G4LogicalVolume*		fLogicNVbarDL;
    G4VPhysicalVolume*		fPhysicNVbarDL;

    G4Tubs*                     fSolidNVdisc;
    G4LogicalVolume*		fLogicNVdisc;
    G4VPhysicalVolume*		fPhysicNVdisc;

    G4Trap*                     fSolidNVwing;
    G4LogicalVolume*		fLogicNVwing;
    G4VPhysicalVolume*		fPhysicNVwing;

};

#endif
/*
 * $Log: DSDetectorWaterTank.hh,v $
 * Revision 1.1  2014/05/07 12:20:51  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.6  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
