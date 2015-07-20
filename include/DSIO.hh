#ifndef _BXIO_HH
#define _BXIO_HH 1

#include "G4String.hh"
#include <iostream>
#include <fstream>

using namespace std;

class DSIO  {
  private:
    DSIO();
  public:
    static DSIO* Get();
    
    virtual ~DSIO() {}

    G4String CheckFileName(G4String); 
    G4String GetFileName()                    { return fFileName ;} 
    void SetFileName(G4String a )             { fFileName = a ;}   
 
      // Binary File   
    ofstream&  GetBinaryFile()                { return fBinaryFile; }
    G4String GetBinaryFileName()              { return fBinaryFileName;}
    void OpenBinaryFile(); 
    void CloseBinaryFile(); 
    G4bool GetIsBinary()                      { return fIsBinary ;} 
    void SetIsBinary(G4bool a)                { fIsBinary = a ;} 
   
      // Log Files    
    G4String   GetLogFileName()               { return fLogFileName;}
    ofstream&  GetStreamLogFile()             { return fStreamLogFile;}
    void OpenLogFiles(); 
    void CloseLogFiles(); 
    
    //G4DS
    ifstream&  GetG4DSFile()                 { return fG4DSFile; }
    G4String GetG4DSFileName()               { return fG4DSFileName;}
    void OpenG4DSFile(); 
    void CloseG4DSFile(); 
    void   SetIsG4DS(G4bool a)               { fIsG4DS = a;}    
    G4bool IsG4DS()                          { return fIsG4DS ;} 

    void    SetDSGeometry(string val)         { fStreamDSGeometryFileName = val ;}
    string  GetDSGeometry()                   { return fStreamDSGeometryFileName ;}
    ifstream&  GetStreamDSGeometry();  
    ifstream&  GetStreamDSCryostatProfile();
    ifstream&  GetStreamDSG3CryostatProfile();

    void CloseStreamDSGeometry();
    void CloseStreamDSCryostatProfile();
    void CloseStreamDSG3CryostatProfile();

    void   SetDSVPMTGeometry(string val)       { fStreamDSVPMTGeometryFileName = val ;}
    string GetDSVPMTGeometry()                 { return fStreamDSVPMTGeometryFileName ;}
    ifstream&  GetStreamDSVPMTGeometry();
    void CloseStreamDSVPMTGeometry();
    
    ifstream&  GetStreamDSOptics();  
    void       SetDSOpticsFileName(string val) { fStreamDSOpticsFileName = val ;}
    void       CloseStreamDSOptics();
 
    void SetG4DSFile(G4String name)     { fG4DSFileName = name;}
 
  private:
  
    static DSIO *me;
   
    ifstream  fStreamDSGeometry;
    G4String  fStreamDSGeometryFileName;

    //optics tuning
    ifstream  fStreamDSOptics;
    G4String  fStreamDSOpticsFileName;
     
    ifstream  fStreamDSCryoProfile;
    G4String  fStreamDSCryoProfFileName; 

    ifstream  fStreamDSG3CryoProfile;
    G4String  fStreamDSG3CryoProfFileName;

    ifstream  fStreamDSVPMTGeometry;
    G4String  fStreamDSVPMTGeometryFileName;

    ofstream  fBinaryFile;    
    ofstream  fStreamLogFile;
    
    G4String  fFileName;
    G4String  fLogFileName;
    G4String  fBinaryFileName;
    G4String  fG4DSFileName;
    ifstream  fG4DSFile;    
    G4bool    fIsBinary ;
    G4bool    fIsG4DS;    
    
    void SetLogFileName(G4String name)        { fLogFileName = name;}
    void SetBinaryFileName(G4String name)     { fBinaryFileName = name;}

    
    void ChangeName(G4String); 
};

#endif
/*
 * $Log: DSIO.hh,v $
 * Revision 1.3  2014/10/13 18:43:46  swesterd
 * fixed veto PMT positions
 *
 * Revision 1.2  2014/07/25 14:07:25  perassos
 * Implementation of the DSG3 (TPC + NV). Configuration #8
 *
 * Revision 1.1  2014/05/07 12:20:53  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.5  2014/03/19 16:44:43  dfranco
 * update external configuration reader for optics tuning
 *
 * Revision 1.4  2014/03/11 16:50:00  dfranco
 * read optical parameters from external file (DSOptics.dat)
 *
 * Revision 1.3  2014/03/11 09:56:25  meregaglia
 * Added generator starting from energy deposits
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
