
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4String.hh"

#include "DSLogger.hh" 
#include "DSIO.hh"
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include "G4String.hh"

using namespace std;

DSIO* DSIO::me = 0;

// singleton
DSIO::DSIO(){

  fIsBinary  = 0;
  fBinaryFileName = "outtest";
  fStreamDSGeometryFileName  = "../data/detector/DSGeometry.dat";
  fStreamDSCryoProfFileName  = "../data/detector/DSCryostatProfiles.dat";
  fStreamDSG3CryoProfFileName  = "../data/detector/DSG3CryostatProfilesFromTechDrawings.dat";
  fStreamDSVPMTGeometryFileName = "../data/detector/VPMTGeometry.dat";
  fStreamDSOpticsFileName  = "../data/detector/DSOptics.dat";

}

DSIO* DSIO::Get() {
  if (!me) 
    me = new DSIO();
       
  return me;
}



G4String DSIO::CheckFileName(G4String mys) {
    std::ifstream file;
    std::ostringstream newfilename ; 
    newfilename  << mys << ".log"   ;
    file.open(newfilename.str().c_str(),std::ios::in);
    if(!file) {
      ChangeName(mys);
      newfilename.str(""); 
      newfilename  << mys;
    } else { 
      file.close();
       for(int i = 1; i< 1E4; i++) {
        newfilename.str("");
        newfilename  <<  mys << "_v" << i << ".log"  ;
        file.open(newfilename.str().c_str (),std::ios::in);
        if(!file) {  
	  newfilename.str("");
	  newfilename  <<  mys << "_v" << i;
	  ChangeName(newfilename.str()); 
	  DSLog(trace)  << "Output files already exist! Name changed: "<<  newfilename.str() << endl ;
	  break ; 
	} file.close();     
      }
    }
    fFileName = newfilename.str();
    return fFileName ;
}

void DSIO::ChangeName(G4String mys) {
  std::ostringstream newfilename ;
  
  newfilename << mys << ".log";
  SetLogFileName(newfilename.str());   
  
  newfilename.str("");
  newfilename << mys << ".fil";
  
  SetBinaryFileName(newfilename.str()); 

}



//--------------------------------------------------------------
// Binary File
//--------------------------------------------------------------

void DSIO::OpenBinaryFile(){
     fBinaryFile.open(fBinaryFileName.c_str(),ofstream::out | ofstream::binary);
}

void DSIO::CloseBinaryFile(){
     //fBinaryFile.flush();
     fBinaryFile.close();
     DSLog(routine) << "Binary File Closed" <<  endlog;      
}


//--------------------------------------------------------------
// Log Files
//--------------------------------------------------------------

void DSIO::OpenLogFiles() {
  if(!fStreamLogFile.is_open()) {
    DSLog(routine) << "Log file created " << endl ;

    fStreamLogFile.open (fLogFileName, ofstream::out);
    fStreamLogFile << "#############  G4DS Log File  #############" << endl ;

  }
} 

void DSIO::CloseLogFiles(){
  fStreamLogFile.close();
  DSLog(routine) << "Log Files Flushed And Closed" <<  endlog;      
}


//--------------------------------------------------------------
// G4DS Files
//--------------------------------------------------------------

void DSIO::OpenG4DSFile(){
     fG4DSFile.open(fG4DSFileName.c_str(),ifstream::in | ifstream::binary);
}

void DSIO::CloseG4DSFile(){
     fG4DSFile.close();
     DSLog(routine) << "G4DS File Closed" <<  endlog;      
}


// DSGeometry
ifstream&  DSIO::GetStreamDSGeometry() {
  if(!fStreamDSGeometry.is_open()) {
    fStreamDSGeometry.open(fStreamDSGeometryFileName,ifstream::in);
  }
  if(fStreamDSGeometry.eof()) {
    fStreamDSGeometry.close();
    fStreamDSGeometry.clear();
    fStreamDSGeometry.open(fStreamDSGeometryFileName,ifstream::in);
  }
  return fStreamDSGeometry ;
}


void DSIO::CloseStreamDSGeometry(){

  if( fStreamDSGeometry.is_open() )
    fStreamDSGeometry.close();
  DSLog(routine) << "Geometry Input File Closed" <<  endlog;      

}

// DSOptics
ifstream&  DSIO::GetStreamDSOptics() {
  if(!fStreamDSOptics.is_open()) {
    fStreamDSOptics.open(fStreamDSOpticsFileName,ifstream::in);
  }
  if(fStreamDSOptics.eof()) {
    fStreamDSOptics.close();
    fStreamDSOptics.clear();
    fStreamDSOptics.open(fStreamDSOpticsFileName,ifstream::in);
  }
  return fStreamDSOptics ;
}

void DSIO::CloseStreamDSOptics(){

  if( fStreamDSOptics.is_open() )
    fStreamDSOptics.close();
  DSLog(routine) << "Optics Input File Closed" <<  endlog;      

}



ifstream& DSIO::GetStreamDSCryostatProfile(){
  if( !fStreamDSCryoProfile.is_open() ) {
    fStreamDSCryoProfile.open(fStreamDSCryoProfFileName,ifstream::in);
  }
  if( fStreamDSCryoProfile.eof() ){
    fStreamDSCryoProfile.close();
    fStreamDSCryoProfile.clear(); 
    fStreamDSCryoProfile.open(fStreamDSCryoProfFileName,ifstream::in);
  }
  return fStreamDSCryoProfile ;

}


ifstream&  DSIO::GetStreamDSVPMTGeometry() {
  if(!fStreamDSVPMTGeometry.is_open()) {
    fStreamDSVPMTGeometry.open(fStreamDSVPMTGeometryFileName,ifstream::in);
  }
  if(fStreamDSVPMTGeometry.eof()) {
    fStreamDSVPMTGeometry.close();
    fStreamDSVPMTGeometry.clear();
    fStreamDSVPMTGeometry.open(fStreamDSVPMTGeometryFileName,ifstream::in);
  }
  return fStreamDSVPMTGeometry ;
}

void DSIO::CloseStreamDSVPMTGeometry(){

  if( fStreamDSVPMTGeometry.is_open() )
    fStreamDSVPMTGeometry.close();
  DSLog(routine) << "Veto PMT Geometry Input File Closed" <<  endlog;      

}


ifstream& DSIO::GetStreamDSG3CryostatProfile(){

  if( !fStreamDSG3CryoProfile.is_open() ) {
    fStreamDSG3CryoProfile.open(fStreamDSG3CryoProfFileName,ifstream::in);
  }
  if( fStreamDSG3CryoProfile.eof() ){
    fStreamDSG3CryoProfile.close();
    fStreamDSG3CryoProfile.clear(); 
    fStreamDSG3CryoProfile.open(fStreamDSG3CryoProfFileName,ifstream::in);
  }
  return fStreamDSG3CryoProfile ;

}


void DSIO::CloseStreamDSCryostatProfile(){

  if( fStreamDSCryoProfile.is_open() )
    fStreamDSCryoProfile.close();
  DSLog(routine) << "Cryostat Profile Input File Closed" <<  endlog;      

}


void DSIO::CloseStreamDSG3CryostatProfile(){

  if( fStreamDSG3CryoProfile.is_open() )
    fStreamDSG3CryoProfile.close();
  DSLog(routine) << "G3 Cryostat Profile Input File Closed" <<  endlog;      

}

/*
 * $Log: DSIO.cc,v $
 * Revision 1.3  2014/10/13 18:43:57  swesterd
 * fixed veto PMT positions
 *
 * Revision 1.2  2014/07/25 14:07:20  perassos
 * Implementation of the DSG3 (TPC + NV). Configuration #8
 *
 * Revision 1.1  2014/05/07 12:21:03  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.5  2014/03/11 16:49:56  dfranco
 * read optical parameters from external file (DSOptics.dat)
 *
 * Revision 1.4  2014/03/11 09:54:38  meregaglia
 * Added generator starting from energy deposits
 *
 * Revision 1.3  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
