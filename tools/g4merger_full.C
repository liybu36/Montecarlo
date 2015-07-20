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
#include "TRandom3.h"
#include "../include/DSEventStructure.hh"


using namespace std;


TRandom3 * ran = new TRandom3();

long int _NEVENTS    = 10000000000 ;
int MAXDEPOSIT  = 20000;


bool isToBeAdded;
// cluser params
float dist_max_xy = 500. ;    //  cm 
float dist_max_z  = 0.2 ;  //  cm 
float deltaT_max  = dist_max_z*10000.  ;  // ns 

int LArIndex ;
int _COUNTER    = 0 ;
bool _NPE       = 0 ;
bool _S1ENE     = 0 ;
ofstream fout;

int ScintillatorIndex ; 

struct Cluster {

  double       genPartPDG ;   // PDG of the particle generating the energy deposit
  double       genPartZ ;     // Z of the particle generating the energy deposit
  double       Length ; 
  double       Radius ;
  double       Energy ;       
  double       kinEne ;       // kinetic energy of the particle generating the energy deposit
  double       dEdx ;
  double       nucl;
  double       elec; 
  double       T0 ;
  float        X0 ;
  float        Y0 ;
  float        Z0 ; 
  float        X1 ; 
  float        Y1 ; 
  float        Z1 ; 
  int          npe ;
  int	       nDeposits;
  int	       nElectrons ;
  int	       nPhotons ;
  int	       nExcitons ;

};

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
vector<Cluster>                Clusters ;

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
///quenching for TPC
double S1quench(double ene, double alpha) {
  
  double F = 200.; 
  double W = 19.5e-6;
  //double alpha = 0.21 ;
  double p[8];
  
  p[0] = 2.65965e-01;
  p[1] = 4.92038e-01;
  p[2] = 1.58933e-04; 
  p[3] = -3.46480e-02;
  p[4] = 7.67703e-01;
  p[5] = 6.31385e-01;
  p[6] = -5.64532e-02;
  p[7] = -7.86833e-06;

  double reco = p[0]*(1 - p[1]*exp((p[2]*F+p[6])*ene*1000))*exp(p[3]*pow(ene*1000,p[4]))+p[5]+p[7]*F;
  double epsilon =1;
  //cout << "reco "<<reco << " "<<ene<<endl ;
  if (ene ==0) return 0;

  double out = (epsilon*ene/W/(1+alpha) *(alpha + reco))*W/ene ;

  //cout<<" out "<<out<<endl;

  return out;//(epsilon*ene/W/(1+alpha) *(alpha + reco))*W/ene ;

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
  LArIndex = theHeader.LArIndex ; 
  
  
  // deposit variables
  int    dep_pdg[MAXDEPOSIT], dep_mat[MAXDEPOSIT], dep_id[MAXDEPOSIT];
  float  dep_ene[MAXDEPOSIT], dep_qene[MAXDEPOSIT],dep_x[MAXDEPOSIT], dep_y[MAXDEPOSIT], dep_Prompt[MAXDEPOSIT] , dep_Del[MAXDEPOSIT] , 
         dep_z[MAXDEPOSIT], dep_r[MAXDEPOSIT], dep_step[MAXDEPOSIT]; 
  double dep_time[MAXDEPOSIT];
  float f90like , cl_ene[MAXDEPOSIT], cl_x[MAXDEPOSIT], cl_y[MAXDEPOSIT], cl_z[MAXDEPOSIT], cl_t[MAXDEPOSIT]; 
  float cl_nucl[MAXDEPOSIT], cl_elec[MAXDEPOSIT];
 
  int ndepoTPC;
  int ndeptpc, nclus, cl_ndep[MAXDEPOSIT], ab_mat, cl_npe[MAXDEPOSIT];
  float s1ene,s2ene, veto_visene, mu_visene, ene, qene, qnpe, tpcene, vetoene,
    muene , depTPCTot , depVeto, npeVeto, npeTPC, npeTPC400, depTPC400 ,
    depVeto70, eneTPC, eneVeto, eneslideVeto, timeslideVeto, npeslideVeto, timeVeto; 

  if(first == 1) {
    int SIZE =  sizeof(HeaderStructure);
    fout.write(reinterpret_cast<char*>(&SIZE)   ,sizeof( int  ));
    fout.write(reinterpret_cast<char*>(&theHeader),SIZE);
    fout.write(reinterpret_cast<char*>(&SIZE)   ,sizeof( int  ));
  }
  for(int i=0; i<_NEVENTS; i++) {  
    if(!_readEvent(_bin_fstream))                    break ; 
    if(!(i % 10000)) cout << i << endl ;
    if(_bin_fstream->eof())                          break ;
    if(_NPE   && int(thePhotoElectrons.size()) == 0) continue ;
    if(_S1ENE && theEvent.S1Energy == 0)             continue ;
    
    
    float VetoLY  = 520.;    //photons per MeV
    float TPCLY   = 7100.;   //photons per MeV
    float qvalpha = 0.027;             // alpha in veto (40keVee for 1.47MeV alphas)
    float qvli    = 0.0 ;               // Li in veto  -  TODO         
    float qvp     = 0.3333;            // p in veto 
    float qarAr   = 0.25 ;             // Ar in LAr
    float quench ; 
  
    depTPCTot = 0;
    depVeto=0;
    npeVeto=0;
    npeTPC=0;
    npeTPC400=0;
    depTPC400=0;
    depVeto70=0;
    ndeptpc = 0 ;
    f90like = -1;
    eneTPC = 0; 
    eneVeto= 0;
    timeVeto = 1e12;
    eneslideVeto = 0;
    timeslideVeto = 1e12;
    npeslideVeto = 0;
   
    // For the Clustering
    Cluster LocalClus;
    Clusters.clear(); 
    
    float  nuceneT = 0 , lepeneT = 0;  
    float nucene, lepene, TotEne;
    float dep_dist_xy = 0,   dep_dist_z  = 0,  deltaT      = 0;     
    isToBeAdded = false ;
    
    int neutron_id = 0 ; 
    double tmp_maxpe = 0, tmp_time_maxpe = 0, t0 = 0;

    bool first = true;

    for(int i=0;i<theEvent.NDeposits;++i) {
      
      dep_id[i]   = i+1 ;
      dep_pdg[i]  = theDeposits[i].PID;
      dep_mat[i]  = theDeposits[i].Volume;
      dep_ene[i]  = theDeposits[i].Energy;
      dep_step[i] = theDeposits[i].Step;
      dep_x[i]    = theDeposits[i].Position[0];
      dep_y[i]    = theDeposits[i].Position[1];
      dep_z[i]    = theDeposits[i].Position[2];
      dep_time[i] = theDeposits[i].Time;
      dep_r[i]    = sqrt(pow(dep_x[i],2)+pow(dep_y[i],2));
      dep_qene[i] = 0;
      
      if(dep_mat[i] == 8 )ndepoTPC++;
      
      if(dep_time[i] > 400e3) continue ; 
      
      // CLUSTERING --> da fare per ogni cluster, energia quenchata + tempo!
      if (dep_mat[i] != LArIndex ) continue ; 
      
      if (dep_x[i]*dep_x[i]+dep_y[i]*dep_y[i] > 17.77*17.77 || dep_z[i] > 14.131 ||  dep_z[i]<-21.439)  continue ;
      
      

      double alpha = 0;

      // TPC quenching
      if(dep_pdg[i] - 1e9 > 100000)
	alpha = 1;
      else
	alpha = 0.21;

      double s1_dep = S1quench(dep_ene[i],alpha);
      double qdepo  = s1_dep * dep_ene[i];
  
      if(dep_pdg[i] -1e9 >  100000) 
	qdepo *= qarAr;     // Lindhard factor
	
	  


      dep_qene[i] = qdepo;

      isToBeAdded = true;
      
      for( int j = 0; j < int( Clusters.size() ); j++) {

       dep_dist_z =   abs( Clusters[j].Z0 - dep_z[i] ) ; 
       deltaT     =   abs( Clusters[j].T0 - dep_time[i] );  
        // clustering conditions
       if( dep_dist_z < dist_max_z &&  deltaT < deltaT_max ) {

         quench = 0 ;
         if (dep_pdg[i] - 1e9 >  100000  )       Clusters[j].nucl += dep_qene[i] ; 
         else if (dep_pdg[i] < 30  )  Clusters[j].elec += dep_qene[i] ; 

         //weighted average
         Clusters[j].Z0 = ( Clusters[j].Z0*Clusters[j].Energy + dep_z[i]*dep_qene[i] ) / (dep_qene[i]+Clusters[j].Energy)  ; 
         Clusters[j].X0 = ( Clusters[j].X0*Clusters[j].Energy + dep_x[i]*dep_qene[i] ) / (dep_qene[i]+Clusters[j].Energy)  ; 
         Clusters[j].Y0 = ( Clusters[j].Y0*Clusters[j].Energy + dep_y[i]*dep_qene[i] ) / (dep_qene[i]+Clusters[j].Energy)  ; 

         Clusters[j].Energy   += dep_qene[i] ;
         Clusters[j].nDeposits++;
 
         isToBeAdded = false;
         break;
       }  
      }   
      // Otherwise, new cluster
      if( isToBeAdded ){
      //    cout << "New Cluster Added" << endl; 
        quench = 0 ;
        if (dep_pdg[i]  >  100000  )      LocalClus.nucl = dep_qene[i] ; 
        else if (dep_pdg[i] - 1e9 < 30  ) LocalClus.elec = dep_qene[i] ; 
        LocalClus.genPartPDG =  dep_pdg[i] ;
        LocalClus.Energy     =  dep_qene[i];
        LocalClus.T0         =  dep_time[i];
        LocalClus.X0         =  dep_x[i] ;
        LocalClus.Y0         =  dep_y[i] ;	  
        LocalClus.Z0         =  dep_z[i] ;	  
        LocalClus.nDeposits  =  1;
        Clusters.push_back( LocalClus ) ;
      }
    }
    
    
    
    /// Cleaning of the clusters
    int count_cl_erased = 0 ;
    int ClToBeErased[200];
    for(int i=0; i< int(Clusters.size()); ++i) {
      if(Clusters[i].Energy*TPCLY < 20) Clusters[i].npe =  ran->Poisson(Clusters[i].Energy*TPCLY)  ;
      else                              Clusters[i].npe =  int(ran->Gaus(Clusters[i].Energy*TPCLY , sqrt(Clusters[i].Energy*TPCLY)) + 0.5) ;
      
      if(Clusters[i].npe < 5) {
        ClToBeErased[count_cl_erased] = i ;
        count_cl_erased++;
      }
    }
    
    for(int i=count_cl_erased-1;i>=0;--i) Clusters.erase(Clusters.begin() + ClToBeErased[i]);
    
    nclus = int(Clusters.size());
    
    if(nclus != 1) continue ;
    
    theEvent.EventID = _COUNTER ;
    writer();
    _COUNTER++;
  } 

  cout << _COUNTER << endl ;
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
 * $Log: g4merger_full.C,v $
 * Revision 1.1  2015/04/29 15:42:22  dfranco
 * added a clustering version of g4merger
 *
 * Revision 1.3  2014/11/12 09:54:22  dfranco
 * update
 *
 * Revision 1.2  2014/11/04 15:21:51  dfranco
 * added a merger for the binary files
 *
 *
 */

