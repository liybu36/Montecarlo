#include <iostream>
#include <sstream>
#include <string>
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"


using namespace std;

void PlotDetectorRhoZ();
void PlotDetectorXY();




void DisplayEvents( char*fname ){


  TFile* f1 = new TFile( fname );
  

  TCanvas* c1 = new TCanvas("c1","canvas",1024,640);
  c1->Divide(2,1);

  c1->cd(1);
  PlotDetectorRhoZ();
  dstree->Draw("dep_z:sqrt(dep_x**2 + dep_y**2)", "", "same");
  

  c1->cd(2);
  PlotDetectorXY();
  dstree->Draw("dep_y:dep_x", "", "same");

}




void DisplayEvent( char* fname, int ev_num ){


  TFile* f1 = new TFile( fname );
  
  ostringstream cond;
  cond << "ev == " << ev_num;


  TCanvas* c1 = new TCanvas("c1","canvas",1024,640);
  c1->Divide(2,1);

  c1->cd(1);
  PlotDetectorRhoZ();
  dstree->Draw("dep_z:sqrt(dep_x**2 + dep_y**2)>>h1", cond.str().c_str(), "same");
  h1->Draw("same");


  c1->cd(2);
  PlotDetectorXY();
  dstree->Draw("dep_y:dep_x", cond.str().c_str(), "same");

}



void PlotDetectorRhoZ() {


  // External and Internal Cryostats
  float extcryo_r[10] = {  0,    7.6, 13.1, 27.7, 32.1,  32.1,  30.1,  24.2,  14.3,    0 };
  float extcryo_z[10] = { 28.4, 28.0, 27.0, 21.4, 14.1, -74.4, -79.7, -83.5, -87.0, -88.8 };

  float intcryo_r[10] = {  0,   10.9, 16.3, 23.2, 25.6,  25.6,  24.0,  14.7,   8.4,   0   };
  float intcryo_z[10] = { 16.2, 15.0, 13.5, 10.0,  1.9, -70.2, -75.0, -80.0, -81.5, -82.2 };

  for(int i = 0; i < 10; i++){
    extcryo_z[i] += 46.4;
    intcryo_z[i] += 46.4;
  }

  TGraph* extcryo = new TGraph( 10, extcryo_r, extcryo_z );
  extcryo->SetLineColor( kGray+1 );
  extcryo->SetLineWidth( 2 );
  extcryo->SetTitle(" z vs rho ");
  extcryo->GetXaxis()->SetTitle("rho [cm]");
  extcryo->GetYaxis()->SetTitle("z [cm]");
  extcryo->Draw("al");

  TGraph* intcryo = new TGraph( 10, intcryo_r, intcryo_z );
  intcryo->SetLineColor( kGray );
  intcryo->SetLineWidth( 2 );
  intcryo->Draw("lsame");



  // Reflector
  //float refl_r[2] = { 17.8,  17.8 };
  float refl_r[5] = { 17.8,  17.8, 20.3, 20.3, 17.8 };
  float refl_z[5] = { -8.6,  26.9, 26.9, -8.6, -8.6 };

  TGraph* refl = new TGraph( 5, refl_r, refl_z );
  refl->SetLineColor( kRed-7 );
  refl->SetLineWidth( 2 );
  refl->Draw("lsame");



  // Bell Top
  float bell_r[5] = {  0,   21.2, 21.2,  0,    0 };
  float bell_z[5] = { 26.9, 26.9, 27.5, 27.5, 26.9 };

  TGraph* bell = new TGraph( 2, bell_r, bell_z );
  bell->SetLineColor( kBlue-7 );
  bell->SetLineWidth( 2 );
  bell->Draw("lsame");



  // Cathode Window
  float cawi_z[5] = { -8.6, -8.6, -9.2, -9.2, -8.6  };

  TGraph* cawi = new TGraph( 2, bell_r, cawi_z );
  cawi->SetLineColor( kBlue-7 );
  cawi->SetLineWidth( 2 );
  cawi->Draw("lsame");



  // Gas Pocket
  float gp_r[4] = {  0,   17.7, 17.7,  0 };
  float gp_z[4] = { 26.8, 26.8, 25.8, 25.8 };

  TGraph* gp = new TGraph( 4, gp_r, gp_z );
  gp->SetLineColor( kBlue-10 );
  gp->SetLineWidth( 2 );
  gp->Draw("lsame");



  // Teflon Caps
  float teflcap_r[4]  = {   0,    23.5,  23.5,   0,  }; 
  float teflcapt_z[4] = {  28.8,  28.8,  35.9,  35.9 };
  float teflcapb_z[4] = { -10.6, -10.6, -17.6, -17.6 };

  TGraph* teflcapt = new TGraph( 4, teflcap_r, teflcapt_z );
  TGraph* teflcapb = new TGraph( 4, teflcap_r, teflcapb_z );
  teflcapt->SetLineColor( kRed-9 );
  teflcapt->SetLineWidth( 2 );
  teflcapt->Draw("lsame");
  teflcapb->SetLineColor( kRed-9 );
  teflcapb->SetLineWidth( 2 );
  teflcapb->Draw("lsame");



  // Teflon Support
  float teflsup_r[5] = {  21.6, 23.5,  23.5,  21.6, 21.6 }; 
  float teflsup_z[5] = {  28.8, 28.8, -10.6, -10.6, 28.8 };

  TGraph* teflsup = new TGraph( 5, teflsup_r, teflsup_z );
  teflsup->SetLineColor( kBlue-9 );
  teflsup->SetLineWidth( 2 );
  teflsup->Draw("lsame");

}




void PlotDetectorXY(){


  const int npts = 50;
  float extcryo_x[npts];
  float extcryo_y[npts];
  float intcryo_x[npts];
  float intcryo_y[npts];
  float etefsup_x[npts];
  float etefsup_y[npts];
  float itefsup_x[npts];
  float itefsup_y[npts];
  float ereflec_x[npts];
  float ereflec_y[npts];
  float ireflec_x[npts];
  float ireflec_y[npts];
  
  float ecr_r = 32.1;
  float icr_r = 25.6;
  float ets_r = 23.5;
  float its_r = 21.6;
  float erf_r = 20.3;
  float irf_r = 17.8;
 
  for( int i = 0; i < npts; i++){

    extcryo_x[i] = ecr_r * cos( (2*float(i)/(npts-1))*acos(-1) );
    extcryo_y[i] = ecr_r * sin( (2*float(i)/(npts-1))*acos(-1) );
    intcryo_x[i] = icr_r * cos( (2*float(i)/(npts-1))*acos(-1) );
    intcryo_y[i] = icr_r * sin( (2*float(i)/(npts-1))*acos(-1) );
    etefsup_x[i] = ets_r * cos( (2*float(i)/(npts-1))*acos(-1) );
    etefsup_y[i] = ets_r * sin( (2*float(i)/(npts-1))*acos(-1) );
    itefsup_x[i] = its_r * cos( (2*float(i)/(npts-1))*acos(-1) );
    itefsup_y[i] = its_r * sin( (2*float(i)/(npts-1))*acos(-1) );
    ereflec_x[i] = erf_r * cos( (2*float(i)/(npts-1))*acos(-1) );
    ereflec_y[i] = erf_r * sin( (2*float(i)/(npts-1))*acos(-1) );
    ireflec_x[i] = irf_r * cos( (2*float(i)/(npts-1))*acos(-1) );
    ireflec_y[i] = irf_r * sin( (2*float(i)/(npts-1))*acos(-1) );

  }

  TGraph* extcryo_xy = new TGraph( npts, extcryo_x, extcryo_y );
  extcryo_xy->SetLineColor( kGray+1 );
  extcryo_xy->SetLineWidth( 2 );
  extcryo_xy->SetTitle(" y vs x ");
  extcryo_xy->GetXaxis()->SetTitle("x [cm]");
  extcryo_xy->GetYaxis()->SetTitle("y [cm]");
  extcryo_xy->Draw("al");

  TGraph* intcryo_xy = new TGraph( npts, intcryo_x, intcryo_y );
  intcryo_xy->SetLineColor( kGray );
  intcryo_xy->SetLineWidth( 2 );
  intcryo_xy->Draw("lsame");

  TGraph* etefsup_xy = new TGraph( npts, etefsup_x, etefsup_y );
  etefsup_xy->SetLineColor( kBlue-9 );
  etefsup_xy->SetLineWidth( 2 );
  etefsup_xy->Draw("lsame");

  TGraph* itefsup_xy = new TGraph( npts, itefsup_x, itefsup_y );
  itefsup_xy->SetLineColor( kBlue-9 );
  itefsup_xy->SetLineWidth( 2 );
  itefsup_xy->Draw("lsame");

  TGraph* ereflec_xy = new TGraph( npts, ereflec_x, ereflec_y );
  ereflec_xy->SetLineColor( kRed-7 );
  ereflec_xy->SetLineWidth( 2 );
  ereflec_xy->Draw("lsame");

  TGraph* ireflec_xy = new TGraph( npts, ireflec_x, ireflec_y );
  ireflec_xy->SetLineColor( kRed-7 );
  ireflec_xy->SetLineWidth( 2 );
  ireflec_xy->Draw("lsame");


}
