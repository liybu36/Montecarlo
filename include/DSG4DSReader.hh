//
//	DSG4DSReader.hh
//

#ifndef DSG4DSReader_h
#define DSG4DSReader_h 1
//#include "DSInputVertex.hh"
//#include "DSOutputStructure.hh"
#include "G4ThreeVector.hh"
#include "DSEventStructure.hh"
#include "globals.hh"
#include <stdio.h>
#include <iostream>

using namespace std;




class DSG4DSReader {
 private:
   DSG4DSReader();
 public:
  ~DSG4DSReader();
  static DSG4DSReader* Get();
  void    ReadHeader();
  G4bool  ReadEvent();
  void DumpHeader();
  void DumpEvent();
  void SetG4DSFile();

  void ClearAll();  


  

  
  // vector<DepositStructure>&           GetVDeposits() { return theDeposits             ; }  
 vector<DepositStructure>           GetVDeposits() { return theDeposits             ; }  
  
  DepositStructure&  GetDeposits()                    { return theDepositStructure		; }
 
  void SetDeposit(DepositStructure val)                            { theDepositStructure           = val;} 
  void    SetDeposits()                        { theDeposits.push_back(theDepositStructure)  ;} 
 

private:

  static DSG4DSReader *me;


  void                 SkipEvents();
  
  EventStructureDiskFormat fEvent;
  DaughterStructure fDaughter;
  DepositStructure fDeposit;
  UserStructure fUser;
  PhotonStructure fPhoton;
  PhotoElectronStructure fPhotonPM;
  PhotoElectronStructure fMuPhoton;
  PhotoElectronStructure fVetoPhoton;


 

  DepositStructure              theDepositStructure;
  vector<DepositStructure>       theDeposits;





};

#endif
/*
 * $Log: DSG4DSReader.hh,v $
 * Revision 1.1  2014/05/07 12:20:52  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.1  2014/03/11 09:56:25  meregaglia
 * Added generator starting from energy deposits
 *
 * Revision 1.1  2009-10-20 12:37:38  dfranco
 * Added a new generator in order to read a g4ds output file and begin a simulation
 * from the energy deposits. Useful for external background.
 *
 *
 */
