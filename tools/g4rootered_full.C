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
#include "TRandom3.h"
#include "TMath.h"
#include "TCanvas.h"
#include "string.h"
#include "TMinuit.h"
#include "TVector3.h"
#include "algorithm"
#include "../include/DSEventStructure.hh"


using namespace std;

TRandom3 * ran = new TRandom3();

int MAXDAUGHTER = 7100;
int MAXDEPOSIT  = 20000;
int MAXNPH      = 70000;
int MAXNPE      = 200000;
int MAXUSER     = 1000;


vector<float> vpdf, vpdf_time, vtheta, vlenght, vjitter_prob, vjitter_time;

//clustering 
bool isToBeAdded;
// cluser params
float dist_max_xy = 500. ;    //  cm 
float dist_max_z  = 0.2 ;  //  cm 
float deltaT_max  = dist_max_z*10000.  ;  // ns 

int ScintillatorIndex ; 
int LArIndex; 

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


bool readNoise(TH1F *hh){
  double npe[12],cum[12];
  double norma  = 0 ;
  double norma2 = 0 ;
  double tmp = 0 ;
 
  ifstream fin("noise_2.dat");
  for(int i=0;i<11;++i) {
    fin >> npe[i] >> tmp ;

    if(i ==0 )    cum[i] = tmp;
    else          cum[i] = cum[i-1] + tmp ;  
    
    npe[i] += 5 ;
    norma  += tmp ;
  }
  fin.close();

  for(int i=0;i<11;++i) {
    cum[i] /= norma ;
    //cout << npe[i] << " " << cum[i] << endl;
  }

  for(int i=0;i<1e6;++i) {
    double myran = ran->Uniform();
    for(int k=0;k<11;++k) 
      
      if(cum[k] > myran ) {
        double mm ;
        double qq ;
        if(k == 0) {
          mm = (npe[k] - 0)/(cum[k] - 0);
	  qq = npe[k] - mm*cum[k];	  
	} else {	
          mm = (npe[k] - npe[k-1])/(cum[k] - cum[k-1]);
	  qq = npe[k] - mm*cum[k];
	}
        hh->Fill(myran*mm + qq);
        break ;
      }
  
  }


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
  for(int i=0; i<theEvent.VetoNPE; i++)    theVetoPhotoElectrons.push_back(_readMuPhotoElectron(file));
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



///quenching factors for VETO
// 18
double MeV = 1000.0;

double proton_energy[] = {0*MeV,0.029*MeV,0.094*MeV,0.2*MeV,0.34*MeV,
                          0.52*MeV, 0.72*MeV, 0.94*MeV,2*MeV,3*MeV,
                          4*MeV,6*MeV,10*MeV,20*MeV,30*MeV,
                          40*MeV,60*MeV,100*MeV};
// 18  
double proton_quenching[] = {0*MeV,0.003*MeV,0.005*MeV,0.009*MeV,0.020*MeV,
                             0.041*MeV,0.071*MeV,0.127*MeV,0.6*MeV,1*MeV,
			     1.6*MeV,3*MeV,6*MeV,13*MeV,20*MeV,30*MeV,45*MeV,70*MeV};
    
// 13
double alpha_energy[] = {0*MeV,.7*MeV,.85*MeV,1*MeV,1.2*MeV,
				       2*MeV,3*MeV,4*MeV,6*MeV,
				       8*MeV,10*MeV,20*MeV,30*MeV};
// 13
double alpha_quenching[] = {0*MeV,.017*MeV,.02*MeV,.025*MeV,.036*MeV,
				.06*MeV,.15*MeV,.22*MeV,.5*MeV,
				1*MeV,1.8*MeV,7*MeV,12*MeV};
    
//quenching for carbon nuclear recoils
//AW from Astropart. Phys. 16 (2002) 333-338 and NIM 33 (1965) 131-135
// 11
double carbon_energy[] = {0*MeV,0.046*MeV, 0.111*MeV, 0.229*MeV,
				       0.368*MeV, 0.500*MeV, 1.2*MeV, 2.0*MeV,
				       3.0*MeV, 4.0*MeV, 5.0*MeV};
// 11
double carbon_quenching[] = {0*MeV,0.0022*MeV, 0.0026*MeV, 0.0032*MeV,
				0.0044*MeV, 0.005*MeV, 0.007*MeV, 0.011*MeV,
				0.018*MeV, 0.025*MeV, 0.035*MeV};


double interpolator(double ene, int dim, double *x, double *y) {
  for(int i=0;i<dim;++i) {
    if(x[i] > ene) {
      double m = (y[i] - y[i-1])/(x[i] - x[i-1]);
      double q = y[i] - m*x[i];
      return m*ene + q ;
    }
  }
  
  return 0;
}

double qalpha(double ene) {
   return interpolator(ene,13,alpha_energy, alpha_quenching)/ene * 0.548 ;
}
double qproton(double ene) {
   return interpolator(ene,18,proton_energy, proton_quenching)/ene  ;
} 
double qcarbon(double ene) {
   return interpolator(ene,11,carbon_energy, carbon_quenching)/ene  ;
}
double qelectron (double ene)/*keV*/ {
  //RNS
  //kB = 0.012 cm/MeV from Borexino
  //Birks Quenching parameterized as in
  //"The ionization quench factor in liquid-scintillation counting standardizations
  //Malonda, Carles
  //Applied Radiation and Isotopes 51 (1999) 183-188

  double A1 = 0.32903;
  double A2 = 0.14404;
  double A3 = 0.08059;
  double A4 = -0.04536E-3;
  double A5 = 0.15623;
  double A6 = 0.06611;
  return ( (A1 + A2*TMath::Log(ene) + A3*TMath::Log(ene)*TMath::Log(ene) + A4*TMath::Log(ene)*TMath::Log(ene)*TMath::Log(ene))/
	   (1 + A5*TMath::Log(ene) + A6*TMath::Log(ene)*TMath::Log(ene) + A4*TMath::Log(ene)*TMath::Log(ene)*TMath::Log(ene))  );
}


double getqfactor(double ene,int pdg) {
  //return 1;
  
  if(pdg - 1e9 == 20040) return qalpha(ene); 
  //if(pdg - 1e9 == 20040) return 1./50.; ///alpha from B has 1.47 MeV, quenching between 30 and 60 keV
  else if(pdg == 11 || pdg == -11 || pdg == 22 ) return qelectron(ene);
  else if(pdg == 2212 || pdg == 2112) return qproton(ene);
  else if(pdg - 1e9 > 30000 ) return qcarbon(ene);
  //else if(pdg - 1e9 > 30000 ) return 0;
  
  
  // else 
  //cout<<"quenching not computed, pdg = "<<" "<<pdg<<endl;

  return 1 ;

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
  int    Did[MAXDAUGHTER], Dpdg[MAXDAUGHTER], Dpid[MAXDAUGHTER], Dprocess[MAXDAUGHTER];
  float  Dene[MAXDAUGHTER], Dx[MAXDAUGHTER], Dy[MAXDAUGHTER], Dz[MAXDAUGHTER], 
         Dr[MAXDAUGHTER], Dpx[MAXDAUGHTER], Dpy[MAXDAUGHTER], Dpz[MAXDAUGHTER];
  double Dtime[MAXDAUGHTER];
  
  // deposit variables
  int    dep_pdg[MAXDEPOSIT], dep_mat[MAXDEPOSIT], dep_id[MAXDEPOSIT];
  float  dep_ene[MAXDEPOSIT], dep_qene[MAXDEPOSIT],dep_x[MAXDEPOSIT], dep_y[MAXDEPOSIT], dep_Prompt[MAXDEPOSIT] , dep_Del[MAXDEPOSIT] , 
         dep_z[MAXDEPOSIT], dep_r[MAXDEPOSIT], dep_step[MAXDEPOSIT]; 
  double dep_time[MAXDEPOSIT];
 
  // user variables
  int    INT1[MAXUSER], INT2[MAXUSER];  
  float  FLOAT1[MAXUSER], FLOAT2[MAXUSER];
  double DOUBLE[MAXUSER];
  
  // pe variables
  int    pe_pmt[MAXNPE];     
  double pe_time[MAXNPE];  

  // veto pe variables
  int    veto_pe_pmt[MAXNPE];    
  double veto_pe_time[MAXNPE]; 
  
  // mu pe variables
  int    mu_pe_pmt[MAXNPE];     
  double mu_pe_time[MAXNPE];   

  // ph variables
  int    ph_volume[MAXNPH], ph_pid[MAXNPH];
  float  ph_wl[MAXNPH], ph_x[MAXNPH], ph_y[MAXNPH], ph_z[MAXNPH];  
  double ph_time[MAXNPH]; 

  int ndepoTPC;
  int ndeptpc, nclus, cl_ndep[MAXDEPOSIT], ab_mat, cl_npe[MAXDEPOSIT];
  float s1ene,s2ene, veto_visene, mu_visene, ene, qene, qnpe, tpcene, vetoene,
    muene , depTPCTot , depVeto, npeVeto, npeTPC, npeTPC400, depTPC400 ,
    depVeto70, eneTPC, eneVeto, eneslideVeto, timeslideVeto, npeslideVeto; 
  double ab_time, timeVeto;
  float f90like , cl_ene[MAXDEPOSIT], cl_x[MAXDEPOSIT], cl_y[MAXDEPOSIT], cl_z[MAXDEPOSIT], cl_t[MAXDEPOSIT]; 
  float cl_nucl[MAXDEPOSIT], cl_elec[MAXDEPOSIT];
  
  double tpromptVeto, tdelayedVeto,tlateVeto;
  //int    epromptVeto, edelayedVeto,elateVeto;
  
  int     prompt_npeVeto    ;
  int     prompt_npeNoise    ;
  double  prompt_timeVeto   ;
  double  prompt_zVeto   ;
  double  prompt_rVeto   ;

  int     late_npeVeto       ;
  double  late_timeVeto    ;
  int     del_npeVeto       ;
  double  del_timeVeto       ;
  
  
  
  TTree *dstree = new TTree("dstree","The G4DS Root Tree");
  dstree->SetMaxVirtualSize(100000);
  
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

  dstree->Branch("npe",            &theEvent.NPE ,	      "npe/I");        
  dstree->Branch("munpe" ,         &theEvent.MuNPE ,	      "munpe/I");      
  dstree->Branch("vnpe",           &theEvent.VetoNPE ,        "vnpe/I");       
  dstree->Branch("nph",            &theEvent.NPH,	      "nph/I");        
  dstree->Branch("ndaughters",     &theEvent.NDaughters,      "ndaughters/I"); 
  dstree->Branch("ndeposits" ,     &theEvent.NDeposits,       "ndeposits/I");  
  dstree->Branch("ndepositsTPC" ,  &ndepoTPC,                 "ndepositsTPC/I");
  dstree->Branch("nusers",         &theEvent.NUsers,	      "nusers/I");    
  
  dstree->Branch("nusers",         &theEvent.NUsers,	      "nusers/I");    
  
  
   

  dstree->Branch("dau_id",         Did,                  "Did[ndaughters]/I");
  dstree->Branch("dau_pdg",        Dpdg,                 "Dpdg[ndaughters]/I");
  dstree->Branch("dau_pid",        Dpid,                 "Dpid[ndaughters]/I");
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

  dstree->Branch("dep_id",      dep_id,                   "dep_id[ndeposits]/I");    
  dstree->Branch("dep_pdg",     dep_pdg,                  "dep_pdg[ndeposits]/I");    
  dstree->Branch("dep_mat",     dep_mat,                  "dep_mat[ndeposits]/I");    
  dstree->Branch("dep_time",    dep_time,                 "dep_time[ndeposits]/D");   
  dstree->Branch("dep_ene",     dep_ene,                  "dep_ene[ndeposits]/F")  ;  
  dstree->Branch("dep_qene",     dep_qene,                  "dep_qene[ndeposits]/F")  ;  
  dstree->Branch("dep_step",    dep_step,                 "dep_step[ndeposits]/F")  ;  
  dstree->Branch("dep_x",       dep_x,                    "dep_x[ndeposits]/F");      
  dstree->Branch("dep_y",       dep_y,                    "dep_y[ndeposits]/F");      
  dstree->Branch("dep_z",       dep_z,                    "dep_z[ndeposits]/F") ;     
  dstree->Branch("dep_r",       dep_r,                    "dep_r[ndeposits]/F") ;     
  dstree->Branch("depTPC",       &depTPCTot,	          "depTPCTot/F") ;
  dstree->Branch("depTPC400",     &depTPC400,	          "depTPC400/F") ;
  dstree->Branch("depVeto",      &depVeto,	          "depVeto/F") ;
  dstree->Branch("npeVeto",      &npeVeto,	          "npeVeto/F") ;
  dstree->Branch("npeTPC",      &npeTPC,	          "npeTPC/F") ;
  dstree->Branch("npeTPC400",      &npeTPC400,	          "npeTPC400/F") ;
  dstree->Branch("depVeto70",   &depVeto70,	          "depVeto70/F") ;
  dstree->Branch("ndeptpc",      &ndeptpc,	          "ndeptpc/I") ;
  dstree->Branch("f90like",      &f90like,	          "f90like/F") ;
  dstree->Branch("ab_time",    &ab_time,                   "ab_time/D");   
  dstree->Branch("ab_mat",    &ab_mat,                     "ab_mat/I");   
  dstree->Branch("nclus",         &nclus,	          "nclus/I") ;
  dstree->Branch("cl_ene",        cl_ene,	          "cl_ene[nclus]/F") ;
  dstree->Branch("cl_ndep",       cl_ndep,	          "cl_ndep[nclus]/I") ;
  dstree->Branch("cl_x",          cl_x,	                   "cl_x[nclus]/F") ;
  dstree->Branch("cl_y",          cl_y,	                   "cl_y[nclus]/F") ;
  dstree->Branch("cl_z",          cl_z,	                   "cl_z[nclus]/F") ;
  dstree->Branch("cl_t",          cl_t,	                  "cl_t[nclus]/F") ;
  dstree->Branch("cl_npe",        cl_npe,	                  "cl_npe[nclus]/I") ;
  dstree->Branch("cl_nucl",        cl_nucl,	                  "cl_nucl[nclus]/F") ;
  dstree->Branch("cl_elec",        cl_elec,	                  "cl_elec[nclus]/F") ;
  dstree->Branch("eneTPC",         &eneTPC,                   "eneTPC/F");    
  
  
  dstree->Branch("prompt_npeVeto",        &prompt_npeVeto,                  "prompt_npeVeto/I"); 
  dstree->Branch("prompt_npeNoise",        &prompt_npeNoise,                  "prompt_npeNoise/I"); 
  dstree->Branch("prompt_timeVeto",       &prompt_timeVeto,                 "prompt_timeVeto/D");      
  dstree->Branch("prompt_zVeto",       &prompt_zVeto,                 "prompt_zVeto/D");      
  dstree->Branch("prompt_rVeto",       &prompt_rVeto,                 "prompt_rVeto/D");      

  dstree->Branch("late_npeVeto",        &late_npeVeto,                  "late_npeVeto/I"); 
  dstree->Branch("late_timeVeto",       &late_timeVeto,                 "late_timeVeto/D");      
  dstree->Branch("del_npeVeto",        &del_npeVeto,                  "del_npeVeto/I"); 
  dstree->Branch("del_timeVeto",       &del_timeVeto,                 "del_timeVeto/D");      


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
  int counter_ene = 0;
 
  ab_time=0; 
  ab_mat = -1 ; 
  // Read Header
  _readHeader(_bin_fstream); 
  LArIndex = theHeader.LArIndex ; 
  ScintillatorIndex = theHeader.ScintillatorIndex ; 
  

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
    
    //if(s1ene < 0.00000001) continue ;
    
    
    // sort photon and photoelectron vectors by time
    std::sort(theDeposits.begin(),theDeposits.end(), cmp_deposit());
    std::sort(thePhotons.begin(),thePhotons.end(), cmp_photon());
    std::sort(thePhotoElectrons.begin(),thePhotoElectrons.end(), cmp_photoelectron());
    std::sort(theVetoPhotoElectrons.begin(),theVetoPhotoElectrons.end(), cmp_photoelectron());
    std::sort(theMuPhotoElectrons.begin(),theMuPhotoElectrons.end(), cmp_photoelectron());
    
    // initialize variables 
    ene         = 0;
    //TPC
    double tp0 = 2.65965e-01;
    double tp1 = 4.92038e-01;
    double tp2 = 1.58933e-04;   
    double tp3 = -3.46480e-02;
    double tp4 = 7.67703e-01;
    double tp5 = 6.31385e-01;
    double tp6 = -5.64532e-02;
    double tp7 = -7.86833e-06;
    double fDriftField = 200;
    
    // Fill daughter variables
    for(int i=0;i<theEvent.NDaughters;++i) {
      Did[i]      = theDaughters[i].Id ;    
      Dpdg[i]     = theDaughters[i].PDG ;   
      Dpid[i]     = theDaughters[i].TrackID ;   
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
    //float VetoLY  = 520.;    //photons per MeV
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

    ndepoTPC = 0;
    
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
    
    
    //if(nclus != 1) continue ;
    //if(nclus == 0) continue ;
  
    double T0;
    if(nclus!=0)
      T0 = Clusters[0].T0 ;
    else
      T0 = 0;
      

    bool   prompt_firstVeto  = true ;

    double tmp_maxpe_late    = 0; 
    double t0_late           = -100;

    double tmp_maxpe_delayed = 0 ;
    double t0_delayed        = -100;
    
    
    double elateVeto        = 0;
    double edelayedVeto      = 0;

    double enelateVeto       = 0 ;
    double enedelayedVeto    = 0 ;
    double delnpeVeto        = 0 ;
    double latenpeVeto        = 0 ;

  
    prompt_npeVeto   = 0;
    prompt_npeNoise   = 0;
    prompt_timeVeto  = -10000;
    prompt_zVeto  = -10000;
    prompt_rVeto  = -10000;
    late_npeVeto     = 0;
    late_timeVeto    = -10000;
    del_npeVeto      = 0;
    del_timeVeto     = -10000;



    for(int k=0;k<theEvent.NDeposits;++k) {

      if(dep_mat[k] != ScintillatorIndex) continue ;
      if(dep_time[k] > 400.e3) break ;
      
      double deltaT = dep_time[k] - T0 ;
      double quench = getqfactor(dep_ene[k],dep_pdg[k]);
      // prompt cut 
      if(deltaT > -10 && deltaT  < 200) {
        if(prompt_firstVeto) {
          prompt_timeVeto = dep_time[k] - T0;
          prompt_firstVeto = false ;
	  prompt_zVeto = dep_z[k];
	  prompt_rVeto = dep_r[k];
        }
        int vnpe = 0 ; 
	double vene = quench*VetoLY*dep_ene[k];
        if ( vene < 20 ) vnpe = ran->Poisson(vene)  ;
        else vnpe =  int(ran->Gaus( vene, sqrt(vene)) + 0.5) ;
	prompt_npeVeto += vnpe ;
	
      }

    
    
       
      
      if(deltaT > 8800) {
        tmp_maxpe_late	   = quench*dep_ene[k];
        t0_late  	   = dep_time[k];
      }
      if(deltaT < 8800 && deltaT > 200) {
        tmp_maxpe_delayed  = quench*dep_ene[k];
        t0_delayed         = dep_time[k];
      }

      for(int kk=k+1;kk<theEvent.NDeposits;++kk){
	
	if(dep_mat[kk] != ScintillatorIndex) continue ;
	
	// late
	if(dep_time[kk] - T0 > 8800)  {
	  if( dep_time[kk] - t0_late < 300.){
	    tmp_maxpe_late += getqfactor(dep_ene[kk],dep_pdg[kk])*dep_ene[kk];
	  } else break;
	}
	
	// delayed
	if(dep_time[kk] - T0 < 8800 && dep_time[kk] - T0 > 200)  {
	  if( dep_time[kk] - t0_delayed < 300.){
	    tmp_maxpe_delayed += getqfactor(dep_ene[kk],dep_pdg[kk])*dep_ene[kk];
	  } else break;
	}
	
      }
      
      
      if(tmp_maxpe_late > elateVeto){
	elateVeto  = tmp_maxpe_late;
	late_timeVeto  = t0_late - T0;
      }  

      if(tmp_maxpe_delayed > edelayedVeto){
	edelayedVeto  = tmp_maxpe_delayed;
	del_timeVeto  = t0_delayed - T0;
      }

     


    }
    
    if(edelayedVeto == 0)              del_npeVeto = 0 ;
    else if(edelayedVeto*VetoLY < 20 ) del_npeVeto = ran->Poisson(edelayedVeto*VetoLY);
    else                                 del_npeVeto = int(ran->Gaus(edelayedVeto*VetoLY,sqrt(edelayedVeto*VetoLY))+0.5 );

    if(elateVeto == 0)                 late_npeVeto = 0;
    else if(elateVeto*VetoLY < 20 )    late_npeVeto = ran->Poisson(elateVeto*VetoLY);
    else                                 late_npeVeto = int(ran->Gaus(elateVeto*VetoLY,sqrt(elateVeto*VetoLY))+0.5 );

    
    for (int i=0; i<Clusters.size() ; ++i) {
      cl_ene[i]   = Clusters[i].Energy; 
      cl_x[i]     = Clusters[i].X0; 
      cl_y[i]     = Clusters[i].Y0; 
      cl_z[i]     = Clusters[i].Z0; 
      cl_t[i]     = Clusters[i].T0; 
      cl_nucl[i]  = Clusters[i].nucl; 
      cl_elec[i]  = Clusters[i].elec; 
      cl_ndep[i]  = Clusters[i].nDeposits;  
      cl_npe[i]   = Clusters[i].npe;

    } 

    // Fill user variables
    for(int i=0;i<theEvent.NUsers;++i) {
      INT1[i]    = theUsers[i].UserInt1;
      INT2[i]    = theUsers[i].UserInt2;
      FLOAT1[i]  = theUsers[i].UserFloat1;
      FLOAT2[i]  = theUsers[i].UserFloat2; 
      DOUBLE[i]  = theUsers[i].UserDouble; 
    }
    // Fill pe variables
    double f90_light = 0.;
    double fall_light = 0.;
    for(int i=0;i<theEvent.NPE;++i) {
      pe_pmt[i]  = thePhotoElectrons[i].PMT;    
      pe_time[i] = thePhotoElectrons[i].Time;    
      //compute f90
      if(pe_time[i]-pe_time[0]<90)
        f90_light+=1;
      if(pe_time[i]-pe_time[0]<7000)
        fall_light+=1;
    }
    if(theEvent.NPE!=0)
      f90like = f90_light/fall_light;
    

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
 * $Log: g4rootered_full.C,v $
 * Revision 1.8  2015/01/22 09:56:46  dfranco
 * add f90
 *
 * Revision 1.7  2015/01/20 13:01:29  dfranco
 * add f90
 *
 * Revision 1.6  2015/01/20 10:57:12  dfranco
 * compute f90
 *
 * Revision 1.5  2014/11/07 15:22:32  dfranco
 * move TPC to right position for clustering
 *
 * Revision 1.4  2014/11/03 15:37:10  dfranco
 * update the clustering algorithm
 *
 * Revision 1.3  2014/10/23 11:40:20  dfranco
 * fix bug in clustering
 *
 * Revision 1.2  2014/07/09 13:06:12  pagnes
 * Generators in materials fixed
 *
 * Revision 1.1  2014/05/08 10:59:17  pagnes
 * Scintillator Index added in binary header
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
