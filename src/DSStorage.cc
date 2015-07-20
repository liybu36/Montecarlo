#include <iostream>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <sstream>
#include <fstream>
#include "TMath.h"
#include "TComplex.h"
#include "DSStorage.hh"
#include "DSLogger.hh"

using namespace std;
using namespace TMath;

DSStorage* DSStorage::me = 0;

// singleton
DSStorage::DSStorage(){
 
  fIsEnDepGenerator    = false;
  fOverlap             = 0;
  fEventCounter        = 1000;
  fExportGDML          = 0;
  fWritePhotons        = 0;
  fVerbosity           = -1;
  fWriteDeposits       = 0;
  fWriteThermalElectrons= 0;
  fWriteDaughters      = 0;
  fPMT                 = 1000000;
  fVetoPMT             = 1000000;
  fMuPMT               = 1000000;
  fRDMDecay            = false ;
  fRealPDGMeanLife     = 0;
  fPreAbsTime          = 0;
  fNDaughters          = 1;
  fRDMChain            = false ;
  fKillS1S2            = false ;
  fKillS2              = false ;
  fKillS1              = false ;
  fScaleS2             = 1.0;
  fDriftField          = 0;  // V/cm
  fExtractionField     = 0;  // V/cm
  fTMBfraction         = 0.5;
  fTimeCut             = 1*ms; 
  fScintillator        = 0;  // 0: BorateScintillator  ; 1 Gd Scintillator
  fFastSimulation      = 0; 
  fIsExternalLArScintillating = false ;
  fVetoYieldFactor     = 1.0;
  fTunedS1At200V       = true ;
  
  fThomasImelNullField = 0.099;
  fThomasImelEp0       = 0.156977;
  fThomasImelEp1       = -0.1;

  fDokeBirksNFp1       =  0.0003; 
  fDokeBirksNFp3       =  0.75;
  fDokeBirksEp1        =  0.07;
  fDokeBirksEp2        = -0.85;
  fDokeBirksEp3        =  0.55;
  
  fHolderRadius        =  60.*cm;   //Default value of source holder position
  fHolderZ             = 0.*cm;
  fHolderPhi           = 0.*degree;
  fHolderSource        = false;

  fIs5KGeometry        = false;
  fIs20KGeometry       = false;
  
  fWL[0] = 300.0; fRI[0] = 2.375; fEC[0] = 0.247 ;
  fWL[1] = 350.0; fRI[1] = 2.292; fEC[1] = 0.082 ;
  fWL[2] = 400.0; fRI[2] = 2.182; fEC[2] = 0.045 ;
  fWL[3] = 450.0; fRI[3] = 2.1  ; fEC[3] = 0.021 ;
  fWL[4] = 500.0; fRI[4] = 2.06 ; fEC[4] = 0.016 ;
  fWL[5] = 550.0; fRI[5] = 2.05 ; fEC[5] = 0.014 ;
  fWL[6] = 600.0; fRI[6] = 2.04 ; fEC[6] = 0.012 ;
  fWL[7] = 650.0; fRI[7] = 2.03 ; fEC[7] = 0.011 ;
  fWL[8] = 700.0; fRI[8] = 2.02 ; fEC[8] = 0.0105;
  fWL[9] = 750.0; fRI[9] = 2.01 ; fEC[9] = 0.0105;
  fWL[10]= 800.0; fRI[10]= 1.914; fEC[10]= 0.01  ;
  
  fSourcePosition  = G4ThreeVector(60.0*cm, 0, 0);
  
  ifstream fin("../data/detector/GridTransm.dat");

  double xx, tt;
  for (int i=0; i < 100; ++i ) {
    fin >> xx >> tt; 
    fAng.push_back( xx/360*2.*TMath::Pi()); 
    // PDM 5/21/2014: just store the value directly.
    //    fTrans.push_back(1-(1-tt)/2.);
    fTrans.push_back( tt );
  }
  fin.close();
}

DSStorage* DSStorage::Get() {
  if (!me) 
    me = new DSStorage();
       
  return me;
}

G4double DSStorage::GetGridParameters(G4double n1, G4double n3, G4double d, G4double lambda, G4double th0) {
  
// PDM 5/21/2014: No redefinition of angle of incidence needed.
//  th0 -= TMath::Pi()/2.;   
  double inc = 0;
  inc = th0;
 
  // PDM: having trouble with normal incidence: returns xa=xb=0.  
  if (inc == 0 ) return 1.0;

  // Interpolate from table
  std::vector<double>::iterator it = fAng.begin();
  std::vector<double>::iterator tr = fTrans.begin(); 
  while ( *it < inc ) {
     ++it; 
     ++tr;
  }
  float xa, xb, ya, yb;
  xb = *it; --it; xa = *it;
  yb = *tr; --tr; ya = *tr; 

  double m = (yb-ya)/(xb-xa);
  fGridParameter = m*inc+ya-xa*m;
  //pdm  std::cout<<"***DSStorage: xa,xb,ya,yb= "<<xa<<"  "<<xb<<"  "<<ya<<"  "<<yb<<std::endl;
  return fGridParameter;
}

G4double *DSStorage::GetITOParameters(G4double n1, G4double n3, G4double d, G4double lambda, G4double th0) {

// PDM 5/21/2014: No redefinition of angle of incidence needed.
//  th0 -= TMath::Pi()/2.;
// n1 materiale 1
// n3 materiale 2
// spessore in nm 
// lambda wl luce incidente
// th0 theta angolo incidenza in radianti rispetto alla normale



  //  Double_t i_n1, i_n3, i_th0, i_d, i_lambda;

  G4double n = 0; 
  G4double k = 0;

  TComplex n2;
  TComplex th1,th2,th3;
  TComplex ct1,st1,ct2,st2,ct3,st3;
  TComplex delta;


  // i_n1 = 1.233; i_n3 = 1.49;
  // i_n=2.182; i_k=0.045;
  // i_d=100.0; i_lambda = 420;
  // i_th0 = 0.0;

  //  std::cout << n << "  " << k << std::endl;

  int np = 11;
  if(lambda<fWL[0] || lambda>fWL[np-1]) {
    fITOParameters[0] = -1;
    fITOParameters[1] = -1;
    fITOParameters[2] = -1;
    return fITOParameters;
  }
  
  for(int i=0; i<np; i++) {
    if( (lambda>=fWL[i]) && (lambda<fWL[i+1]) ) {
      n = fRI[i+1]+(fRI[i]-fRI[i+1])*(lambda-fWL[i+1])/(fWL[i]-fWL[i+1]);
      k = fEC[i+1]+(fEC[i]-fEC[i+1])*(lambda-fWL[i+1])/(fWL[i]-fWL[i+1]);
      break;
    }
  }

  n2(n, k);

  th1(th0, 0.0);
  // PDM 7/25/14: Debugging to make this code give same answers as stand-alone Matlab code it was based on.
  //  Main problem was in definition of ct3, which is imaginary when there is total internal reflection.
  //  Secondary problem was in deriving ct2 and ct3 from st2 and st3, where there is apparently a difference in the sign
  //  convention between the root TComplex class and Matlab.  I got around this by an empirical fix (note the Conjugate 
  //  in ct3 but NOT in ct2), that I compared to the Matlab result for many cases.

  //  th2 = TComplex::ASin((n1*TMath::Sin(th0))/n2);
  //  th3(TMath::ASin((n1*TMath::Sin(th0))/n3), 0.0);

  ct1 = TComplex::Cos(th1);
  st1 = TComplex::Sin(th1);

  //  ct2 = TComplex::Cos(th2);
  //  st2 = TComplex::Sin(th2);
  st2 = st1*(n1/n2);
  ct2 = TComplex::Sqrt(1.-st2*st2);

  //  ct3 = TComplex::Cos(th3);
  //  st3 = TComplex::Sin(th3);
  st3 = st1*(n1/n3);
  ct3 = TComplex::Sqrt(1.-st3*st3);
  ct3 = TComplex::Conjugate(ct3);

  delta = (2.0*TMath::Pi()*d/lambda) * (n2 * ct2);

  TComplex t12p,t23p,r12p,r23p;
  TComplex t12s,t23s,r12s,r23s;
  
  t12s = (2.0 * n1 * ct1) / (n1*ct1 + n2*ct2);
  t23s = (2.0 * n2 * ct2) / (n2*ct2 + n3*ct3);

  t12p = (2.0 * n1 * ct1) / (n2*ct1 + n1*ct2);
  t23p = (2.0 * n2 * ct2) / (n3*ct2 + n2*ct3);

  r12s = (n1*ct1 - n2*ct2) / (n1*ct1 + n2*ct2);
  r23s = (n2*ct2 - n3*ct3) / (n2*ct2 + n3*ct3);
  
  r12p = (n2*ct1 - n1*ct2) / (n2*ct1 + n1*ct2);
  r23p = (n3*ct2 - n2*ct3) / (n3*ct2 + n2*ct3);

  TComplex M11p,M21p;
  TComplex Rptotal, Tptotal;

  M11p = TComplex::Exp(-TComplex::I()*delta) + r12p*r23p*TComplex::Exp(TComplex::I()*delta);
  M11p /= t12p*t23p;

  M21p = r12p*TComplex::Exp(-TComplex::I()*delta) + r23p*TComplex::Exp(TComplex::I()*delta);
  M21p /= t12p*t23p;

  Rptotal = M21p/M11p;
  Tptotal = 1.0/M11p;

  TComplex M11s,M21s;
  TComplex Rstotal, Tstotal;

  M11s = TComplex::Exp(-TComplex::I()*delta) + r12s*r23s*TComplex::Exp(TComplex::I()*delta);
  M11s /= t12s*t23s;

  M21s = r12s*TComplex::Exp(-TComplex::I()*delta) + r23s*TComplex::Exp(TComplex::I()*delta);
  M21s /= t12s*t23s;

  Rstotal = M21s/M11s;
  Tstotal = 1.0/M11s;

  Double_t Rs,Rp,Ts,Tp;

  Rs = Rstotal.Rho2();
  Rp = Rptotal.Rho2();
  Ts = Tstotal.Rho2() * (n3*ct3.Re())/(n1*ct1.Re());
  Tp = Tptotal.Rho2() * (n3*TComplex::Conjugate(ct3).Re())/(n1*ct1.Re());

  G4double RR = (Rs+Rp)/2.0;    // riflettivita'   -> riflessione speculare
  G4double TT = (Ts+Tp)/2.0;    // trasmissivita'  -> stessa direzione
  G4double AA = 1.0 - (TT+RR);  // aa assorbimento -> muore


  fITOParameters[0] = RR;
  fITOParameters[1] = TT;
  fITOParameters[2] = AA;
  return fITOParameters;

}



/*
 * $Log: DSStorage.cc,v $
 * Revision 1.12  2015/04/23 14:04:02  pagnes
 * DS20K geometry added (config 10)
 *
 * Revision 1.11  2015/01/14 16:58:35  dfranco
 * added source vial and correspondent commands by Laura and Bianca; manual updated
 *
 * Revision 1.10  2014/12/22 14:40:43  dfranco
 * added the option to activate the recombination probability at 200 V/cm (/ds/physics/tunedS1); this option is by default true; selecting a specific drift field automatically switch off the tunedS1 option
 *
 * Revision 1.9  2014/11/21 10:18:59  dfranco
 * added a command to scale the veto scintillation yield factor and fixed the visible energy variable in the veto
 *
 * Revision 1.8  2014/11/20 15:32:05  dfranco
 * added a command to remove scintillation process from liquid argon between TPC and cryostat
 *
 * Revision 1.7  2014/11/06 17:39:45  dfranco
 * added source and fixed a bug: removed gaseous nitrogen from the neutron veto
 *
 * Revision 1.6  2014/07/25 17:15:39  meyers
 * Fixed a bug in ITO code.  Now gives same results as the standalone code it was based on
 *
 * Revision 1.5  2014/07/23 14:52:41  pagnes
 * write thermal e- and kill S1 commands added
 *
 * Revision 1.4  2014/07/16 08:23:04  pagnes
 * QE scaling to 1.0 added (/ds/manager/fast_simulation xxx)
 *
 * Revision 1.3  2014/06/03 13:31:35  meyers
 * Migrate TPC grid and ITO optics updates to g4ds10
 *
 * Revision 1.17  2014/05/28 14:33:40  meyers
 * Remove 90-degree rotation in ITO angle of incidence
 *
 * Revision 1.16  2014/05/23 17:55:49  meyers
 * Update grid optical model
 *
 * Revision 1.15  2014/04/11 12:33:29  pagnes
 * command to set TMB/PC ratio inside veto added
 *
 * Revision 1.14  2014/03/19 16:37:27  dfranco
 * update external configuration reader for optics tuning
 *
 * Revision 1.13  2014/03/11 09:54:38  meregaglia
 * Added generator starting from energy deposits
 *
 * Revision 1.12  2014/01/29 13:13:41  perassos
 * Update of the electric field handling and of the Nuclear Recoils generator
 *
 * Revision 1.11  2014/01/07 14:10:36  perassos
 * Added the commands to set in the macfile the electric field and the Thomas-Imel parameters
 *
 * Revision 1.10  2013/11/19 10:33:20  perassos
 * Added methods to handle the electric field and the liquid/gas interface z coordinate
 *
 * Revision 1.9  2013/07/24 09:49:02  dfranco
 * Added S1 and S2 equivalent energies, the material index for each deposit, the command killS1S2 to kill photons and electrons generated by DSLight (after storing the equivalent energies)
 *
 * Revision 1.8  2013/06/15 07:16:33  dfranco
 * Added spacename TMath for using TComplex for the ITO optical surface
 *
 * Revision 1.7  2013/06/11 22:48:33  dfranco
 * Added ITO optical boundary. The Sernelius function is defined in DSStorage, and called by G4OpBoundaryProcess. New ITO bool variable added to DSDetectorDS50.cc as surface property (G4MaterialPropertyTable)
 *
 * Revision 1.6  2013/06/10 14:15:39  dfranco
 * Added two commands: /ds/physics/killS2 and /ds/physics/scaleS2 to kill or scale the S2 light
 *
 * Revision 1.5  2013/04/03 10:14:25  dfranco
 * Fixed bugs with RDM and RDMChain staking actions. The logic of Geant4 is changed. Different excited states of a nucleus correspond to new particles (trackID). Code adapted.
 *
 * Revision 1.4  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
