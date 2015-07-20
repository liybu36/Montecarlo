#ifndef DSDetectorSource_H
#define DSDetectorSource_H

#include "G4LogicalVolume.hh"
#include "G4Polycone.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Tubs.hh"
#include "G4VPhysicalVolume.hh"

#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4MaterialPropertiesTable.hh"


class	DSMaterial;

class   DSDetectorSource  {

  public:

    DSDetectorSource(G4VPhysicalVolume*);
   ~DSDetectorSource();
    
    void SetSourcePosition(G4ThreeVector val)  { fPosition = val; }
    G4ThreeVector GetSourcePosition()          { return fPosition ; }
    


  private:
    void                        DefineSurfaces();

    G4VPhysicalVolume*          fMotherVolume;
    G4ThreeVector               fPosition ;
    
    
};

#endif
/*
 * $Log: DSDetectorSource.hh,v $
 * Revision 1.1  2014/11/06 17:39:51  dfranco
 * added source and fixed a bug: removed gaseous nitrogen from the neutron veto
 *
 *
 */
