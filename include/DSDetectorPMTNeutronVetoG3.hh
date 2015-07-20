#ifndef DSDetectorPMTNeutronVetoG3_H
#define DSDetectorPMTNeutronVetoG3_H

#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4Box.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Paraboloid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"


class   DSDetectorPMTNeutronVetoG3  {

  public:

    DSDetectorPMTNeutronVetoG3(G4VPhysicalVolume*);
   ~DSDetectorPMTNeutronVetoG3();

    G4VPhysicalVolume* GetDetectorComponent()  { return  fPhysicPhotocathode;}
    

  private:
  
    void                        DefineSurfaces();

    G4VPhysicalVolume*          fMotherVolume;

    G4Paraboloid*               fSolidPhotocathode;
    G4LogicalVolume*            fLogicPhotocathode;
    G4VPhysicalVolume* 		      fPhysicPhotocathode;
 
    G4Paraboloid*               fSolidPMThalfbulb;
    G4UnionSolid*               fSolidPMTbulb;
    G4LogicalVolume*            fLogicPMTbulb;
    G4VPhysicalVolume*		      fPhysicPMTbulb;

    G4Tubs*                     fSolidPMTcyl;
    G4SubtractionSolid*         fSolidPMTbase; 
    G4LogicalVolume*            fLogicPMTbase;
    G4VPhysicalVolume*		      fPhysicPMTbase;

    

    G4OpticalSurface*           fOpVPhotocathodeSurface;
    G4OpticalSurface*           fOpLumirrorSurface;
    G4OpticalSurface*           fOpPMTBackSurface;
  
    G4LogicalBorderSurface*     fPMTbaseSurface;
    G4LogicalBorderSurface*     fPMTbulbSurface;
    G4LogicalBorderSurface*     fPhotocathodeSurface;
  
};

#endif
/*
 * $Log: DSDetectorPMTNeutronVetoG3.hh,v $
 * Revision 1.1  2014/07/25 14:07:25  perassos
 * Implementation of the DSG3 (TPC + NV). Configuration #8
 *
 * Revision 1.1  2014/05/07 12:20:51  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.6  2013/05/25 07:58:22  swesterd
 * Got the veto PMT optical boundaries all working along with photocathode optical properties, added PMT quantum efficiency to DSTrackingAction, and added a function to DSTrackingAction that locates and quadratically interpolates points in data, for getting useful QEs
 *
 * Revision 1.5  2013/05/20 13:58:38  swesterd
 * updated the veto PMT geometry
 *
 * Revision 1.4  2013/05/14 04:20:21  swesterd
 * Fixed the veto PMT geometry
 *
 * Revision 1.3  2013/03/22 14:09:39  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
