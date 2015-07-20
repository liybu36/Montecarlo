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

int MAXDAUGHTER = 70000;
int MAXDEPOSIT  = 100000;
int MAXNPH      = 100000;
int MAXNPE      = 2000000;
int MAXUSER     = 1000;


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
  
  if(theEvent.NDaughters > MAXDAUGHTER) { cout << "Fatal: NDaughters > MAXDAUGHTER : " << theEvent.NDaughters << " > " << MAXDAUGHTER << endl ; exit(0) ;}
  if(theEvent.NDeposits  > MAXDEPOSIT ) { cout << "Fatal: NDeposits = " << theEvent.NDeposits << " > MAXDEPOSIT" << endl ; exit(0) ;}
  if(theEvent.NUsers     > MAXUSER)     { cout << "Fatal: NUsers > MAXUSER" << endl ; exit(0) ;}
  if(theEvent.NPH        > MAXNPH)      { cout << "Fatal: NPH > MAXNPH" << endl ; exit(0) ;}
  if(theEvent.NPE        > MAXNPE)      { cout << "Fatal: NPE > MAXNPE" << endl ; exit(0) ;}
  if(theEvent.VetoNPE    > MAXNPE)      { cout << "Fatal: VetoNPE > MAXNPE" << endl ; exit(0) ;}
  if(theEvent.MuNPE      > MAXNPE)      { cout << "Fatal: MuNPE > MAXNPE" << endl ; exit(0) ;} 


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
//                    main 
//--------------------------------------------
int main (int argc, char *argv[]) {

  if(argc == 1 || (argc > 1 && !string(argv[1]).find("help")) )  { 
  //if(file.find("help") < 10) {
    cout << "Usage: g4rooter [FILE] [OPTIONS] [OUTPUT]" <<endl ;
    cout <<endl ;
    cout << " Options: " << endl ;
    cout << " nevents=N:  max number (N) of events to process (default: 10000)" <<endl ;    
    cout << " skipN=N:    skip the first N events (default: 0)" <<endl ;    
    cout << " kB=xxx:     set the Birks parameter (default: 0.012 cm/MeV)" <<endl ;    
    cout << " LY=xxx:     set the light yield (default: 500 p.e./MeV)" <<endl ;    
    cout << endl;
    cout << " Output: " << endl ;
    cout << " filename.root (default: FILE.root)" << endl ;
    cout << endl;
    cout << "Version v783.r4" << endl ;
    cout << "..bye ;) DF  (dfranco@in2p3.fr)" << endl ;
    return 0 ;
  }

  string file = argv[1] ;

  if(file.find(".fil") == string::npos) {
    cout << "file " <<  file << " not found.... Bye!" << endl;
    cout << "ps: the input file needs the .fil extension"  << endl;
    return 0;
  }
  string rootfile = file;
  rootfile.replace(file.find(".fil"),4,".root");

  int nevents     = 100000000;
  int skipNevents = 0;
  float kB        = 0.012;
  float LY        = 500.;

  int loop=2;
  while(argv[loop]) {
    string argument = argv[loop];
    if(!argument.find("nevents=")) {  
      argument.erase(0,8);  
      nevents = (int)  atoi(argument.c_str()) ;
      cout << "   max number of events to process: " << nevents << endl;
    }    
    if(!argument.find("skipN=")) {  
      argument.erase(0,6);  
      skipNevents = atoi(argument.c_str()) ;
      cout << "   events to skip: " << skipNevents << endl;
    }    
    if(!argument.find("kB=")) {  
      argument.erase(0,3);  
      kB = atof(argument.c_str()) ;
      cout << "   kB: " << kB << " cm/MeV" << endl;
    }    
    if(!argument.find("LY=")) {  
      argument.erase(0,3);  
      kB = atof(argument.c_str()) ;
      cout << "   LY: " << LY << " p.e./MeV" << endl;
    }    
    size_t found;
    found=argument.rfind(".root");
    if(found!=string::npos) { 
      rootfile = argument;
      cout << "   output root file: " << argument << endl ;
    }
    
    loop++;
  }
  
  TFile *ff = new TFile(rootfile.c_str(),"recreate");
  
  // event extra variables
  float radius = 0;
  
  // daughter variables
  int    Did[MAXDAUGHTER], Dpdg[MAXDAUGHTER], Dtrackid[MAXDAUGHTER], Dparenttrackid[MAXDAUGHTER],   // Dpid[MAXDAUGHTER],
    Dprocess[MAXDAUGHTER], Dmat[MAXDAUGHTER];
  float  Dene[MAXDAUGHTER], Dx[MAXDAUGHTER], Dy[MAXDAUGHTER], Dz[MAXDAUGHTER], 
         Dr[MAXDAUGHTER], Dpx[MAXDAUGHTER], Dpy[MAXDAUGHTER], Dpz[MAXDAUGHTER];
  double Dtime[MAXDAUGHTER];
  
  // deposit variables


/*  int    dep_pdg[MAXDEPOSIT], dep_mat[MAXDEPOSIT], dep_track[MAXDEPOSIT];
  float  dep_ene[MAXDEPOSIT], dep_x[MAXDEPOSIT], dep_y[MAXDEPOSIT], 
         dep_z[MAXDEPOSIT], dep_r[MAXDEPOSIT], dep_step[MAXDEPOSIT]; 
  double dep_time[MAXDEPOSIT];
*/

  int    *dep_pdg   = new int[MAXDEPOSIT]; 
  int    *dep_mat   = new int[MAXDEPOSIT];
  int    *dep_track = new int[MAXDEPOSIT];
  int    *dep_parenttrack = new int[MAXDEPOSIT];
  float  *dep_totalene  = new float[MAXDEPOSIT];
  float  *dep_kineticene  = new float[MAXDEPOSIT];
  float  *dep_ene   = new float[MAXDEPOSIT]; 
  float  *dep_x     = new float[MAXDEPOSIT];
  float  *dep_y     = new float[MAXDEPOSIT]; 
  float  *dep_z     = new float[MAXDEPOSIT]; 
  float  *dep_r     = new float[MAXDEPOSIT]; 
  float  *dep_step  = new float[MAXDEPOSIT]; 
  double *dep_time  = new double[MAXDEPOSIT];
   

  // user variables
  int    INT1[MAXUSER], INT2[MAXUSER];  
  float  FLOAT1[MAXUSER], FLOAT2[MAXUSER];
  double DOUBLE[MAXUSER];
  
  // pe variables
  int    *pe_pmt = new int[MAXNPE];     
  double *pe_time= new double[MAXNPE]; 
  int NPES1, NPES2 ; 

  // veto pe variables
  int    *veto_pe_pmt  = new int[MAXNPE];    
  double *veto_pe_time = new double[MAXNPE]; 
  
  // mu pe variables
  int    *mu_pe_pmt  = new int[MAXNPE];     
  double *mu_pe_time = new double[MAXNPE];   

  // ph variables
  int    ph_volume[MAXNPH], ph_pid[MAXNPH];
  float  ph_wl[MAXNPH], ph_x[MAXNPH], ph_y[MAXNPH], ph_z[MAXNPH];  
  double ph_time[MAXNPH]; 

  float s1ene,s2ene, veto_visene, mu_visene, ene, qene, qnpe, tpcene, vetoene, muene ;
  TTree *dstree = new TTree("dstree","The G4DS Root Tree");

  //  dstree->SetMaxVirtualSize(100000);
  
  dstree->SetMaxVirtualSize(10000000000);
  dstree->Branch("ev",             &theEvent.EventID,         "ev/I");
  dstree->Branch("pdg",            &theEvent.PDG,             "pdg/I");
  dstree->Branch("ene0",           &theEvent.Energy,          "ene0/F");  
  dstree->Branch("s1ene",          &theEvent.S1Energy,        "s1ene/F");     
  dstree->Branch("s2ene",          &theEvent.S2Energy,        "s2ene/F");     
  dstree->Branch("veto_visene",    &theEvent.VetoVisEnergy,   "veto_visene/F");    
  dstree->Branch("mu_visene",      &theEvent.MuVisEnergy,     "mu_visene/F");      

  dstree->Branch("tpcene",          &theEvent.TPCDepEnergy,   "tpcene/F");     
  dstree->Branch("vetoene",          &theEvent.VetoDepEnergy, "vetoene/F");     
  dstree->Branch("muene",          &theEvent.MuDepEnergy,     "muene/F");     
  dstree->Branch("ene",            &ene,		      "ene/F");        
  dstree->Branch("x",              &theEvent.Position[0],     "x/F"); 	       
  dstree->Branch("y",              &theEvent.Position[1],     "y/F"); 	       
  dstree->Branch("z",              &theEvent.Position[2],     "z/F"); 	       
  dstree->Branch("r",              &radius,		      "radius/F");     
  dstree->Branch("px",             &theEvent.Direction[0],    "px/F");	       
  dstree->Branch("py",             &theEvent.Direction[1],    "py/F");	       
  dstree->Branch("pz",             &theEvent.Direction[2],    "pz/F");	       
  dstree->Branch("bx",             &theEvent.CenterOfMass[0], "bx/F");	       
  dstree->Branch("by",             &theEvent.CenterOfMass[1], "by/F");	       
  dstree->Branch("bz",             &theEvent.CenterOfMass[2], "bz/F");	       

  dstree->Branch("npeS1",          &NPES1 ,	              "npeS1/I");        
  dstree->Branch("npeS2",          &NPES2 ,        	      "npeS2/I");        
  dstree->Branch("npe",            &theEvent.NPE ,	      "npe/I");        
  dstree->Branch("munpe" ,         &theEvent.MuNPE ,	      "munpe/I");      
  dstree->Branch("vnpe",           &theEvent.VetoNPE ,        "vnpe/I");       
  dstree->Branch("nph",            &theEvent.NPH,	      "nph/I");        
  dstree->Branch("ndaughters",     &theEvent.NDaughters,      "ndaughters/I"); 
  dstree->Branch("ndeposits" ,     &theEvent.NDeposits,       "ndeposits/I");  
  dstree->Branch("nusers",         &theEvent.NUsers,	      "nusers/I");     
  dstree->Branch("dau_id",         Did,                  "Did[ndaughters]/I");
  dstree->Branch("dau_pdg",        Dpdg,                 "Dpdg[ndaughters]/I");
  //  dstree->Branch("dau_pid",        Dpid,                 "Dpid[ndaughters]/I");
  dstree->Branch("dau_trackid",    Dtrackid,              "Dtrackid[ndaughters]/I");
  dstree->Branch("dau_parenttrackid",    Dparenttrackid,      "Dparenttrackid[ndaughters]/I");
  dstree->Branch("dau_process",    Dprocess,             "Dprocess[ndaughters]/I");
  dstree->Branch("dau_time",       Dtime,                "Dtime[ndaughters]/D");
  dstree->Branch("dau_ene",        Dene,                 "Dene[ndaughters]/F");   
  dstree->Branch("dau_x",          Dx,                   "Dx[ndaughters]/F");	 
  dstree->Branch("dau_y",          Dy,                   "Dy[ndaughters]/F");	 
  dstree->Branch("dau_z",          Dz,                   "Dz[ndaughters]/F") ;   
  dstree->Branch("dau_r",          Dr,                   "Dr[ndaughters]/F");	 
  dstree->Branch("dau_px",         Dpx,                  "Dpx[ndaughters]/F");   
  dstree->Branch("dau_py",         Dpy,                  "Dpy[ndaughters]/F") ;  
  dstree->Branch("dau_pz",         Dpz,                  "Dpz[ndaughters]/F");   
  dstree->Branch("dep_pdg",     dep_pdg,                  "dep_pdg[ndeposits]/I");    
  dstree->Branch("dep_mat",     dep_mat,                  "dep_mat[ndeposits]/I");    

  dstree->Branch("dep_track",   dep_track,                "dep_track[ndeposits]/I");    
  dstree->Branch("dep_parenttrack",   dep_parenttrack,                "dep_parenttrack[ndeposits]/I");    


  dstree->Branch("dep_time",    dep_time,                 "dep_time[ndeposits]/D");   
  dstree->Branch("dep_totalene",     dep_totalene,                  "dep_totalene[ndeposits]/F")  ;  
  dstree->Branch("dep_kineticene",     dep_kineticene,                  "dep_kineticene[ndeposits]/F")  ;  
  dstree->Branch("dep_ene",     dep_ene,                  "dep_ene[ndeposits]/F")  ;  
  dstree->Branch("dep_step",    dep_step,                 "dep_step[ndeposits]/F")  ;  
  dstree->Branch("dep_x",       dep_x,                    "dep_x[ndeposits]/F");      
  dstree->Branch("dep_y",       dep_y,                    "dep_y[ndeposits]/F");      
  dstree->Branch("dep_z",       dep_z,                    "dep_z[ndeposits]/F") ;     
  dstree->Branch("dep_r",       dep_r,                    "dep_r[ndeposits]/F") ;     
 
  dstree->Branch("userint1",    INT1,                     "int1[nusers]/I");  
  dstree->Branch("userint2",    INT2,                     "int2[nusers]/I");
  dstree->Branch("userfloat1",  FLOAT1,                   "float1[nusers]/F")  ; 
  dstree->Branch("userfloat2",  FLOAT2,                   "float2[nusers]/F");   
  dstree->Branch("userdouble0", DOUBLE,                   "double0[nusers]/D");  
 
  dstree->Branch("pe_time",     pe_time,                  "pe_time[npe]/D");
  dstree->Branch("pe_pmt",      pe_pmt,                   "pe_pmt[npe]/I");   
  
  dstree->Branch("vpe_time",  veto_pe_time,               "veto_pe_time[vnpe]/D");   
  dstree->Branch("vpe_pmt",   veto_pe_pmt,                "veto_pe_pmt[vnpe]/I"); 
    
  dstree->Branch("mupe_time",   mu_pe_time,               "mu_pe_time[munpe]/D");   
  dstree->Branch("mupe_pmt",    mu_pe_pmt,                "mu_pe_pmt[munpe]/I");   
  
  dstree->Branch("ph_volume",    ph_volume,               "ph_volume[nph]/I");
  dstree->Branch("ph_pid",       ph_pid,                  "ph_pid[nph]/I");
  dstree->Branch("ph_wl",        ph_wl,                   "ph_wl[nph]/F");
  dstree->Branch("ph_x",         ph_x,                    "ph_x[nph]/F");
  dstree->Branch("ph_y",         ph_y,                    "ph_y[nph]/F");
  dstree->Branch("ph_z",         ph_z,                    "ph_z[nph]/F");
  dstree->Branch("ph_time",      ph_time,                 "ph_time[nph]/D");
  
  
  // Open the binary file
  ifstream *_bin_fstream;
  _bin_fstream = new ifstream (file.c_str(), std::ios::binary);
  cout << endl ;
  cout << "Binary File: " << file << endl;
  cout << endl ;
  if (!(*_bin_fstream)) {
    std::cerr << "Cannot open file. Exiting...\n";
    exit(1);
  }
  
  int counter = 0;
  
  // Read Header
  _readHeader(_bin_fstream); 
  
  
  // Loop over the events
  for(int _i=0; _i<nevents; _i++) {  
  
    // Skip Events
    if(_i < skipNevents) {
      _skipEvent(_bin_fstream);
      continue;
    }
    // Read Event
    if(!_readEvent(_bin_fstream)) break ; 

    // Close the binary file if the end is reached
    if(_bin_fstream->eof()) break ;
    
    // Print counter

    if (!(counter % 100000) && counter>0) std::cout << counter <<" processed events " << " (event id = " <<  theEvent.EventID << ")"<< std::endl;
    counter++ ;
    
    // sort photon and photoelectron vectors by time
    std::sort(theDeposits.begin(),theDeposits.end(), cmp_deposit());
    std::sort(thePhotons.begin(),thePhotons.end(), cmp_photon());
    std::sort(thePhotoElectrons.begin(),thePhotoElectrons.end(), cmp_photoelectron());
    std::sort(theVetoPhotoElectrons.begin(),theVetoPhotoElectrons.end(), cmp_photoelectron());
    std::sort(theMuPhotoElectrons.begin(),theMuPhotoElectrons.end(), cmp_photoelectron());
    
    // initialize variables 
    ene         = 0;
    
   
    // Fill daughter variables
    for(int i=0;i<theEvent.NDaughters;++i) {
      Did[i]      = theDaughters[i].Id ;    
      Dpdg[i]     = theDaughters[i].PDG ;   
      //     Dpid[i]     = theDaughters[i].PID ;   
      Dtrackid[i]     = theDaughters[i].TrackID ;   
      Dparenttrackid[i]     = theDaughters[i].ParentTrackID ;   
      Dprocess[i] = theDaughters[i].Process ;
      Dtime[i]    = theDaughters[i].Time ;  
      Dene[i]     = theDaughters[i].Energy ;   
      Dx[i]       = theDaughters[i].Position[0] ;     
      Dy[i]       = theDaughters[i].Position[1] ;     
      Dz[i]       = theDaughters[i].Position[2] ;     
      Dpx[i]      = theDaughters[i].Direction[0] ;    
      Dpy[i]      = theDaughters[i].Direction[1] ;
      Dpz[i]      = theDaughters[i].Direction[2] ;
      Dr[i]       = sqrt(pow(Dx[i],2)+pow(Dy[i],2)+pow(Dz[i],2));
    }
    // Fill deposits variables
    for(int i=0;i<theEvent.NDeposits;++i) {
      dep_pdg[i]  = theDeposits[i].PID;
      dep_mat[i]  = theDeposits[i].Volume;
      dep_track[i]= theDeposits[i].Track;
      dep_parenttrack[i]= theDeposits[i].ParentTrack;
      dep_totalene[i]  = theDeposits[i].TotalEnergy;
      dep_kineticene[i]  = theDeposits[i].KineticEnergy;
      dep_ene[i]  = theDeposits[i].Energy;
      dep_step[i] = theDeposits[i].Step;
      dep_x[i]    = theDeposits[i].Position[0];
      dep_y[i]    = theDeposits[i].Position[1];
      dep_z[i]    = theDeposits[i].Position[2];
      dep_time[i] = theDeposits[i].Time;
      dep_r[i]    = sqrt(pow(dep_x[i],2)+pow(dep_y[i],2)+pow(dep_z[i],2));
      if(dep_mat[i] == 6) ene         += dep_ene[i];
    }    
    // Fill user variables
    for(int i=0;i<theEvent.NUsers;++i) {
      INT1[i]    = theUsers[i].UserInt1;
      INT2[i]    = theUsers[i].UserInt2;
      FLOAT1[i]  = theUsers[i].UserFloat1;
      FLOAT2[i]  = theUsers[i].UserFloat2; 
      DOUBLE[i]  = theUsers[i].UserDouble; 
    }
    NPES1 = 0 ; 
    NPES2 = 0 ; 
    // Fill pe variables
    for(int i=0;i<theEvent.NPE;++i) {
      pe_pmt[i]  = thePhotoElectrons[i].PMT;    
      pe_time[i] = thePhotoElectrons[i].Time;    
      if (pe_time[i] < 10000) NPES1++ ; 
      else NPES2++ ; 
    }
    // Fill veto pe variables
    for(int i=0;i<theEvent.VetoNPE;++i) {
      veto_pe_pmt[i]  = theVetoPhotoElectrons[i].PMT;    
      veto_pe_time[i] = theVetoPhotoElectrons[i].Time;    
    }    
     // Fill mu pe variables
    for(int i=0;i<theEvent.MuNPE;++i) {
      mu_pe_pmt[i]  = theMuPhotoElectrons[i].PMT;    
      mu_pe_time[i] = theMuPhotoElectrons[i].Time;    
    }
    // Fill photon variables
    for(int i=0;i<theEvent.NPH;++i) {
      ph_volume[i] = thePhotons[i].VolumeID; 
      ph_pid[i]    = thePhotons[i].PID; 
      ph_wl[i]     = thePhotons[i].Wavelength;  
      ph_x[i]      = thePhotons[i].Position[0];   
      ph_y[i]      = thePhotons[i].Position[1];   
      ph_z[i]      = thePhotons[i].Position[2];   
      ph_time[i]   = thePhotons[i].Time;    
    }
        
    // Fill the event
    dstree->Fill();
    
  }
  // Write the tree
  dstree->Write();
  
  // Close the root file
  ff->Close();

  // Close the binary file
  _bin_fstream->close();
  cout << "Rootfile " << rootfile.c_str() << " created! " << endl ;
  cout << "Bye...!" << endl ;
  return 0 ;
}


/*
 * $Log: g4rootered.C,v $
 * Revision 1.5  2014/11/13 16:47:11  dfranco
 * removed variables which were creating conflicts with the previous version of g4ds10
 *
 * Revision 1.4  2014/10/21 15:57:24  pagnes
 * g4rootered typos
 *
 * Revision 1.3  2014/10/13 18:43:59  swesterd
 * fixed veto PMT positions
 *
 * Revision 1.2  2014/07/23 14:54:25  pagnes
 * Updated array declaration and new (trivial) npeS1-npeS2 variables added
 *
 * Revision 1.1  2014/05/07 12:21:08  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.14  2013/10/20 16:30:20  swesterd
 * updated the veto scintillator optical properties based on Aldos measurements
 *
 * Revision 1.13  2013/10/01 06:26:25  swesterd
 * added waveform averager
 *
 * Revision 1.12  2013/08/27 04:07:02  swesterd
 * some fine tuning of the boron scintillator kB and scint yield, and some modifications to the DSG2 geometry
 *
 * Revision 1.11  2013/08/20 03:25:53  swesterd
 * added G2 TPC geoemtry (not complete) and added monoenergetic energy distribution to generator
 *
 * Revision 1.10  2013/08/06 13:58:22  dfranco
 * Added 3 variables to store the deposited energy in LAr, scintillator, and water. The last two are not yet implemented. g4rooter has been updated with 3 new variables: tpcene, vetoene, and muene
 *
 * Revision 1.9  2013/07/24 09:49:03  dfranco
 * Added S1 and S2 equivalent energies, the material index for each deposit, the command killS1S2 to kill photons and electrons generated by DSLight (after storing the equivalent energies)
 *
 * Revision 1.8  2013/06/13 10:04:58  dfranco
 * update g4rooter variable names
 *
 * Revision 1.7  2013/06/04 01:02:31  swesterd
 * other than the optical boundary of the trunks, the veto optics appear to be complete and up and running...modulo whatever I may have missed...
 *
 * Revision 1.6  2013/05/07 23:06:32  swesterd
 * added optical boundaries and Lumirror in the veto
 *
 * Revision 1.5  2013/04/04 09:22:54  dfranco
 * added variables to dstree
 *
 * Revision 1.4  2013/04/04 09:07:47  dfranco
 * added deposit step length, visibile (quenched) energy, original energy to dstree
 *
 * Revision 1.3  2013/03/22 13:23:21  dfranco
 * deposit times sorted
 *
 * Revision 1.2  2013/03/20 09:54:28  dfranco
 * New version of g4ds
 *
 */
