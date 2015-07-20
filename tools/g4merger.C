// --------------------------------------------------------------------------//
/** 
 * AUTHOR: D. Franco
 * CONTACT: dfranco@in2p3.fr
 * 
 * Generate a root file reading the binary file from the g4ds output
*/
// --------------------------------------------------------------------------//
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TObject.h"
#include "TRandom.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TH1.h"
#include "TRandom.h"
#include "TMath.h"
#include "TCanvas.h"
#include "string.h"
#include "TMinuit.h"
#include "TVector3.h"
#include "algorithm"
#include "../include/DSEventStructure.hh"


using namespace std;


long int _NEVENTS    = 10000000000 ;
int _COUNTER         = 0 ;
int _GLOBAL_COUNTER  = 0 ;
bool _NPE            = 0 ;
bool _S1ENE          = 0 ;
ofstream fout;

vector<float> vpdf, vpdf_time, vtheta, vlenght, vjitter_prob, vjitter_time;

HeaderStructure                theHeader;
EventStructureDiskFormat       theEvent;
vector<DepositStructure>       theDeposits;
vector<DaughterStructure>      theDaughters;
vector<UserStructure>          theUsers ;
vector<PhotonStructure>        thePhotons;
vector<PhotoElectronStructure> thePhotoElectrons ;
vector<PhotoElectronStructure> theVetoPhotoElectrons ;
vector<PhotoElectronStructure> theMuPhotoElectrons ;

struct cmp_photoelectron{ 
  bool operator() (const PhotoElectronStructure& a, const PhotoElectronStructure& b) {
    return a.Time<b.Time ;
  }
};
struct cmp_photon{ 
  bool operator() (const PhotonStructure& a, const PhotonStructure& b) {
    return a.Time<b.Time ;
  }
};
struct cmp_deposit{ 
  bool operator() (const DepositStructure& a, const DepositStructure& b) {
    return a.Time<b.Time ;
  }
};

bool _readHeader (ifstream *file) {
  int event_size;
  int event_size2;
  file->read ((char *)(&event_size), sizeof (int));  
  file->read ((char *)(&theHeader), sizeof( HeaderStructure)  );  
  file->read ((char *)(&event_size2), sizeof (int));  
  if (!(*file)) return false;
  return true;
} 

DepositStructure _readDeposit (ifstream *file) {
  DepositStructure theDeposit ;
  file->read ((char *)(&theDeposit), sizeof( DepositStructure) );  
  return theDeposit;
} 

DaughterStructure _readDaughter (ifstream *file) {
  DaughterStructure theDaughter ;
  file->read ((char *)(&theDaughter), sizeof( DaughterStructure) );  
  return theDaughter;
} 

UserStructure _readUser (ifstream *file) {
  UserStructure theUser ;
  file->read ((char *)(&theUser), sizeof( UserStructure) );  
  return theUser;
} 

PhotonStructure _readPhoton (ifstream *file) {
  PhotonStructure thePhoton ;
  file->read ((char *)(&thePhoton), sizeof( PhotonStructure) );       
  return thePhoton;
} 

PhotoElectronStructure _readPhotoElectron (ifstream *file) {
  PhotoElectronStructure thePhotoElectron ;
  file->read ((char *)(&thePhotoElectron), sizeof( PhotoElectronStructure) );       
  return thePhotoElectron;
}
 
PhotoElectronStructure _readVetoPhotoElectron (ifstream *file) {
  PhotoElectronStructure theVetoPhotoElectron ;
  file->read ((char *)(&theVetoPhotoElectron), sizeof( PhotoElectronStructure) );       
  return theVetoPhotoElectron;
} 

PhotoElectronStructure _readMuPhotoElectron (ifstream *file) {
  PhotoElectronStructure theMuPhotoElectron ;
  file->read ((char *)(&theMuPhotoElectron), sizeof( PhotoElectronStructure) );       
  return theMuPhotoElectron;
} 



bool _readEvent (ifstream *file) {
  theDeposits.clear();
  theDaughters.clear();
  theUsers.clear() ;
  thePhotons.clear();
  thePhotoElectrons.clear() ;
  theVetoPhotoElectrons.clear() ;
  theMuPhotoElectrons.clear() ;   

  if (!(*file)) return false;
  
  int event_size;
  int event_size2;
  file->read ((char *)(&event_size), sizeof (int)); 
  file->read ((char *)(&theEvent), sizeof( EventStructureDiskFormat)  ); 
  if(file->eof()) return false ;
  


  for(int i=0; i<theEvent.NDaughters; i++) theDaughters.push_back(_readDaughter(file));
  for(int i=0; i<theEvent.NDeposits; i++)  theDeposits.push_back(_readDeposit(file));
  for(int i=0; i<theEvent.NUsers; i++)     theUsers.push_back(_readUser(file));
  for(int i=0; i<theEvent.NPH; i++)        thePhotons.push_back(_readPhoton(file));
  for(int i=0; i<theEvent.NPE; i++)        thePhotoElectrons.push_back(_readPhotoElectron(file));
  for(int i=0; i<theEvent.VetoNPE; i++)    theVetoPhotoElectrons .push_back(_readMuPhotoElectron(file));
  for(int i=0; i<theEvent.MuNPE; i++)      theMuPhotoElectrons.push_back(_readMuPhotoElectron(file));
  file->read ((char *)(&event_size2), sizeof (int));  
  if(file->eof()) return false ;
  if(event_size != event_size2) return false ;
  //  cout << "Problem at event " << theEvent.EventID << endl;
  //}


  return true;

}

bool _skipEvent (ifstream *file) {
  int event_size;
  int event_size2;
  file->read ((char *)(&event_size), sizeof (int));
  char * buffer = new char [event_size]; 
  file->read (buffer,event_size);
  file->read ((char *)(&event_size2), sizeof (int));  
  if(file->eof()) return false ;
  return true;
}
//--------------------------------------------
//                    writer 
//--------------------------------------------
 
int writer() {

  int SIZE =   sizeof(EventStructureDiskFormat) 
             + theDaughters.size()*sizeof(DaughterStructure)
             + theDeposits.size() *sizeof(DepositStructure)
             + theUsers.size()    *sizeof(UserStructure)
             + thePhotons.size()	   *sizeof(PhotonStructure)
             + thePhotoElectrons.size()	   *sizeof(PhotoElectronStructure)
             + theVetoPhotoElectrons.size()   *sizeof(PhotoElectronStructure)
             + theMuPhotoElectrons.size()     *sizeof(PhotoElectronStructure);


  fout.write(reinterpret_cast<char*>(&SIZE)   ,sizeof( int  ));
  fout.write(reinterpret_cast<char*>(&theEvent)                   ,sizeof(EventStructureDiskFormat));
  for(int i = 0; i <theEvent.NDaughters; i++ ) 
     fout.write(reinterpret_cast<char*>(&theDaughters[i])        ,sizeof(DaughterStructure));
  for(int i = 0; i <theEvent.NDeposits; i++ ) 
     fout.write(reinterpret_cast<char*>(&theDeposits[i])         ,sizeof(DepositStructure));       
  for(int i = 0; i <theEvent.NUsers; i++ ) 
    fout.write(reinterpret_cast<char*>(&theUsers[i])             ,sizeof(UserStructure));
  for(int i = 0; i <theEvent.NPH; i++ ) 
    fout.write(reinterpret_cast<char*>(&thePhotons[i])           ,sizeof(PhotonStructure));
  for(int i = 0; i <theEvent.NPE; i++ ) 
    fout.write(reinterpret_cast<char*>(&thePhotoElectrons[i])    ,sizeof(PhotoElectronStructure));    
  for(int i = 0; i <theEvent.MuNPE; i++ ) 
    fout.write(reinterpret_cast<char*>(&theMuPhotoElectrons[i])  ,sizeof(PhotoElectronStructure));    
  for(int i = 0; i <theEvent.VetoNPE; i++ ) 
    fout.write(reinterpret_cast<char*>(&theVetoPhotoElectrons[i]),sizeof(PhotoElectronStructure));    
  fout.write(reinterpret_cast<char*>(&SIZE) ,sizeof( int    ));
  

}

//--------------------------------------------
//                    exe 
//--------------------------------------------
void exe (char filename[], int first) {
  ifstream *_bin_fstream = new ifstream (filename, std::ios::binary);
  _readHeader(_bin_fstream); 
  
  if(first == 1) {
    int SIZE =  sizeof(HeaderStructure);
    fout.write(reinterpret_cast<char*>(&SIZE)   ,sizeof( int  ));
    fout.write(reinterpret_cast<char*>(&theHeader),SIZE);
    fout.write(reinterpret_cast<char*>(&SIZE)   ,sizeof( int  ));
  }
  for(int i=0; i<_NEVENTS; i++) { 
    ++_GLOBAL_COUNTER ;  
    if(!_readEvent(_bin_fstream))          break ; 
    if(!(i % 10000)) cout << i << endl ;
    if(_bin_fstream->eof())                break ;
    if(_NPE   && int(thePhotoElectrons.size()) == 0) continue ;
    if(_S1ENE && theEvent.S1Energy == 0) continue ;
    theEvent.EventID = _COUNTER ;
    writer();
    _COUNTER++;
  } 

  cout << _COUNTER <<"/" <<_GLOBAL_COUNTER<<"=" <<float(_COUNTER)/_GLOBAL_COUNTER << endl ;
}


//--------------------------------------------
//                    main 
//--------------------------------------------
int main (int argc, char *argv[]) {

  if(argc == 1 || (argc > 1 && !string(argv[1]).find("help")) )  { 
  //if(file.find("help") < 10) {
    cout << "Usage: g4merger [OUTPUT] [OPTIONS] [FILE1] [FILE2] ... [FILEN]" <<endl ;
    cout <<endl ;
    cout << " Options: " << endl ;
    cout << " nevents=N:  max number (N) of events to process per file (default: 1e8)" <<endl ;    
    cout << " s1cut:   remove events with s1ene = 0" <<endl ;    
    cout << " npecut:  remove events with npe = 0" <<endl ;    
    cout << endl;
    cout << "Version v783.r4" << endl ;
    cout << "..bye ;) DF  (dfranco@in2p3.fr)" << endl ;
    return 0 ;
  }

  for(int i=2;i<argc;++i) {
    string argument = argv[i];
    if(argument.find(".fil") == string::npos) {
      if(!argument.find("nevents=")) {  
	argument.erase(0,8);  
	_NEVENTS = (int)  atoi(argument.c_str()) ;
	cout << "   max number of events to process per file: " << _NEVENTS << endl;
      } else if(!argument.find("npecut")) _NPE = 1;
      else if(!argument.find("s1cut")) _S1ENE = 1;
    }
  }

  fout.open(argv[1],ofstream::out | ofstream::binary);
  int first = 1;
  for(int i=2;i<argc;++i) {
    string filename = argv[i];
    if(filename.find(".fil") != string::npos) {
      cout << "filename: " << filename << endl ;
      exe(argv[i], first);
      first = 0 ;
    }
  }
  fout.close();
  
  cout << "number of written events: " << _COUNTER << endl ;
  
}

/*
 * $Log: g4merger.C,v $
 * Revision 1.4  2015/05/07 09:24:09  pagnes
 * configDarkSide.sh - More recent RadioactDecay and PhotEvap libs added in the non-specific cluster section
 *
 * Revision 1.3  2014/11/12 09:54:22  dfranco
 * update
 *
 * Revision 1.2  2014/11/04 15:21:51  dfranco
 * added a merger for the binary files
 *
 *
 */

