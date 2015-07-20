#ifndef DSVisManager_h
#define DSVisManager_h 1

#ifdef G4VIS_USE

#include "G4VisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DSVisManager: public G4VisManager {

public:

  DSVisManager ();

private:

  void RegisterGraphicsSystems ();

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

#endif
/*
 * $Log: DSVisManager.hh,v $
 * Revision 1.1  2014/05/07 12:20:56  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.1  2013/03/25 16:56:10  dfranco
 * Added visualization class and macro
 *
 *
 */
