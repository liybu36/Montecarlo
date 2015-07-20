#ifndef DSMaterial_h
#define DSMaterial_h

#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"

class DSMaterial  {

  private:
    DSMaterial();
    static DSMaterial *me;
  
  
  public:
    static DSMaterial* Get();
   ~DSMaterial();

    G4Material*  GetMaterial(G4String ); 
    void DefineMaterials();
    void DefineProperties();
    void DefineBoundaries();
    void SetMaterialIndexes(); 

    G4Material*  GetPPO()                                  { return fPPO		   ;}		 	  
    G4Material*  GetDMP()                                  { return fDMP		   ;}		 	  
    G4Material*  GetPC()                                   { return fPC 		   ;}		 	   
    G4Material*  GetWater()                                { return fWater		   ;}		 	  
    G4Material*  GetGlass()                                { return fGlass		   ;}		  
    G4Material*  GetSteel()                                { return fSteel		   ;}		  
    G4Material*  GetNylon()                                { return fNylon		   ;}		  
    G4Material*  GetNylonL()                               { return fNylonL		   ;}		 	  
    G4Material*  GetAir()                                  { return fAir		   ;}		 	  
    G4Material*  GetVacuum()                               { return fVacuum		   ;}		 	  
    G4Material*  GetBialkali()                             { return fBialkali		   ;}		 	  
    G4Material*  GetPaint()                                { return fPaint		   ;}		 	  
    G4Material*  GetBorexScintillator()                    { return fBorexScintillator     ;}		  
    G4Material*  GetDMPbuffer()                            { return fDMPbuffer  	   ;}		    
    G4Material*  GetQuartz()                               { return fQuartz		   ;}		 	  
    G4Material*  GetBe()                                   { return fBe 		   ;}	  	  
    G4Material*  GetPb()                                   { return fPb 		   ;}	  	  
    G4Material*  GetDerlin()                               { return fDerlin		   ;}		 	  
    G4Material*  GetRock()                                 { return fRock		   ;}		 	  
    G4Material*  GetConcrete()                             { return fConcrete		   ;}		 	  
    G4Material*  GetAluminum()                             { return fAluminum		   ;}		 	  
    G4Material*  GetAcrylic()                              { return fAcrylic		   ;}		 	  
    G4Material*  GetLiquidArgon()                          { return fLiquidArgon	   ;}		   
    G4Material*  GetNSLiquidArgon()                        { return fNSLiquidArgon	   ;}		   
    G4Material*  GetGaseousArgon()                         { return fGaseousArgon	   ;}		    
    G4Material*  GetLiquidXenon()                          { return fLiquidXenon	   ;}		    
    G4Material*  GetGaseousXenon()                         { return fGaseousXenon	   ;}		    
    G4Material*  GetStainlessSteel()                       { return fStainlessSteel	   ;}		    
    G4Material*  GetGridSteel()                            { return fGridSteel		   ;}		  
    G4Material*  GetMetalLead()                            { return fMetalLead  	   ;}		    	  
    G4Material*  GetTeflon()                               { return fTeflon		   ;}	  
    G4Material*  GetNorite()                               { return fNorite		   ;}	  
    G4Material*  GetHDPE()                                 { return fHDPE		   ;}		    	  
    G4Material*  GetStyrofoam()                            { return fStyrofoam  	   ;}		    	  
    G4Material*  GetMetalCopper()                          { return fMetalCopper	   ;}		    
    G4Material*  GetBoronScintillator()                    { return fBoronScintillator     ;}		    
    G4Material*  GetGdScintillator()                       { return fGdScintillator        ;}  
    G4Material*  GetTMB()                                  { return fTMB		   ;}		    	   
    G4Material*  GetShipSteel()                            { return fShipSteel  	   ;}		    	  
    G4Material*  GetTPB()                                  { return fTPB		   ;}		    	   
    G4Material*  GetThreeMFoil()                           { return fThreeMFoil 	   ;}		    
    G4Material*  GetFusedSilica()                          { return fFusedSilica	   ;}		    
    G4Material*  GetSodiumIodide()                         { return fSodiumIodide	   ;}		    
    G4Material*  GetITO()                                  { return fITO		   ;}		    	   
    G4Material*  GetLiquidNitrogen()                       { return fLiquidNitrogen	   ;}		    
    G4Material*  GetGaseousNitrogen()                      { return fGaseousNitrogen       ;}
    G4Material*  GetMetalTitanium()                        { return fMetalTitanium	   ;}		    
    G4Material*  GetKovar()                                { return fKovar		   ;}		    	  
    G4Material*  GetOldBialkali()                          { return fOldBialkali	   ;}		    
    G4Material*  GetKapton()                               { return fKapton		   ;}		    	  
    G4Material*  GetBlackHole()                            { return fBlackHole             ;}  
    G4Material*  GetPseudoArgon()                          { return fPseudoArgon           ;} 
    G4Material*  GetMetalSilicon()                         { return fMetalSilicon          ;} 
    G4Material*  GetGdWater()                              { return fGdWater               ;} 
     
    //Optical Boundary MPTs
    G4MaterialPropertiesTable* GetLumirrorMPT()                      { return fLumirrorMPT                      ;}
    G4MaterialPropertiesTable* GetElectropolishedStainlessSteelMPT() { return fElectropolishedStainlessSteelMPT ;}
    G4MaterialPropertiesTable* GetUntreatedStainlessSteelMPT()       { return fUntreatedStainlessSteelMPT       ;}
    G4MaterialPropertiesTable* GetAluminumFoilMPT()                  { return fAluminumFoilMPT                  ;}
    G4MaterialPropertiesTable* GetVPhotocathodeMPT()                 { return fVPhotocathodeMPT                 ;}
    G4MaterialPropertiesTable* GetPhotocathodeMPT()                  { return fPhotocathodeMPT                  ;}
    G4MaterialPropertiesTable* GetPMTBackMPT()                       { return fPMTBackMPT                       ;}

  private:
    
    G4double     GetLArRefIndex(G4double);
    G4double     GetGArRefIndex(G4double);
    G4double     GetLArEpsilon(G4double);
    G4double     GetGArEpsilon(G4double);
    G4double     GetLArRayLength(G4double);
    G4double     GetGArRayLength(G4double);
    
    G4Material*           fPPO;                       
    G4Material*           fDMP;                       
    G4Material*           fPC;                        
    G4Material*           fWater;                      
    G4Material*           fGlass;                      
    G4Material*           fSteel;                      
    G4Material*           fNylon;                      
    G4Material*           fNylonL;                      
    G4Material*           fAir;                       
    G4Material*           fVacuum;                      
    G4Material*           fBialkali;                      
    G4Material*           fPaint;                      
    G4Material*           fBorexScintillator;              
    G4Material*           fDMPbuffer;                      
    G4Material*           fQuartz;                      
    G4Material*           fBe;                        
    G4Material*           fPb;                        
    G4Material*           fDerlin;                      
    G4Material*           fRock;                      
    G4Material*           fConcrete;                      
    G4Material*           fAluminum;                      
    G4Material*           fAcrylic;                      
    G4Material*           fLiquidArgon; 
    G4Material*           fNSLiquidArgon;          
    G4Material*           fGaseousArgon;              
    G4Material*           fLiquidXenon;               
    G4Material*           fGaseousXenon;              
    G4Material*           fStainlessSteel;              
    G4Material*           fGridSteel;              
    G4Material*           fMetalLead;                      
    G4Material*           fTeflon;                      
    G4Material*           fNorite;                      
    G4Material*           fHDPE;                      
    G4Material*           fStyrofoam;                      
    G4Material*           fMetalCopper;               
    G4Material*           fBoronScintillator;              
    G4Material*           fTMB;                       
    G4Material*           fShipSteel;                      
    G4Material*           fTPB;                       
    G4Material*           fThreeMFoil;                
    G4Material*           fFusedSilica;               
    G4Material*           fSodiumIodide;              
    G4Material*           fITO;                       
    G4Material*           fLiquidNitrogen;              
    G4Material*           fGaseousNitrogen;
    G4Material*           fMetalTitanium;              
    G4Material*           fKovar;                      
    G4Material*           fOldBialkali;               
    G4Material*           fKapton; 
    G4Material*           fBlackHole;                  
    G4Material*           fGdOxide;                  
    G4Material*           fGdScintillator;
    G4Material*           fPseudoArgon;
    G4Material*           fMetalSilicon;
    G4Material*           fGdWater;
    
    // Material properties table
    G4MaterialPropertiesTable *fLumirrorMPT;
    G4MaterialPropertiesTable *fElectropolishedStainlessSteelMPT;
    G4MaterialPropertiesTable *fUntreatedStainlessSteelMPT;
    G4MaterialPropertiesTable *fAluminumFoilMPT;
    G4MaterialPropertiesTable *fVPhotocathodeMPT;
    G4MaterialPropertiesTable *fPhotocathodeMPT;
    G4MaterialPropertiesTable *fPMTBackMPT;
};

#endif
/*
 * $Log: DSMaterial.hh,v $
 * Revision 1.6  2015/04/17 14:51:29  dfranco
 * added water loaded with gadolinium
 *
 * Revision 1.5  2015/03/09 15:20:42  pagnes
 * DS 5tons geometry added (conf 9)
 *
 * Revision 1.4  2015/01/17 11:31:45  pagnes
 * PAr model added form optical tuning
 *
 * Revision 1.3  2014/11/20 15:32:12  dfranco
 * added a command to remove scintillation process from liquid argon between TPC and cryostat
 *
 * Revision 1.2  2014/05/07 14:27:26  dfranco
 * fixed some bugs and added GdScintillator
 *
 * Revision 1.1  2014/05/07 12:20:54  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.11  2014/04/10 15:15:11  pagnes
 * TMB added to BoronScintillator materialsrc/DSMaterial.cc
 *
 * Revision 1.10  2013/06/21 13:09:21  dfranco
 * Added DS10 optical surfaces
 *
 * Revision 1.9  2013/06/19 18:35:25  swesterd
 * added DSScintCelll and made tpc PMTs' QE and reflections work like veto PMTs
 *
 * Revision 1.8  2013/06/05 23:03:28  swesterd
 * moved optical boundary MPTs to DSMaterial and gave the trunks optical boundary properties consistent with untreated stainless steel
 *
 * Revision 1.7  2013/06/04 16:56:42  dfranco
 * Added fake BlackHole material to abosrb photons
 *
 * Revision 1.6  2013/04/19 16:24:14  dfranco
 * Added Rayleigh scattering to liquid and gaseous argon
 *
 * Revision 1.5  2013/04/19 13:36:49  dfranco
 * Added absorption length
 *
 * Revision 1.4  2013/03/22 17:48:02  dfranco
 * added refrective indexes for gasous and liquid argon materials
 *
 * Revision 1.3  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
