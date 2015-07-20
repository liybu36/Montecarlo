#ifndef _DSStorage_HH
#define _DSStorage_HH 1

#include "G4ThreeVector.hh"
#include "DSLogger.hh"
#include "G4String.hh"
#include "TComplex.h"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

class DSStorage  {
  private:
    DSStorage();
  public:
    static DSStorage* Get();
    
    virtual ~DSStorage() {}
    
     
  private:
  
    static DSStorage *me;
    
  
  public:
  
    void SetIsEnDepGenerator(G4bool val)  { fIsEnDepGenerator = val;  }
    G4int GetIsEnDepGenerator()          { return fIsEnDepGenerator; }           

    void SetCheckOverlap(G4int val)  { fOverlap = val;  }
    G4int GetCheckOverlap()          { return fOverlap; }           

    void SetEventCounter(G4int val)  { fEventCounter = val;  }
    G4int GetEventCounter()          { return fEventCounter; }           

    void SetExportGDML(G4bool val)   { fExportGDML = val;  }
    G4bool GetExportGDML()           { return fExportGDML; }           
   
    void SetWritePhotons(G4bool val)  { fWritePhotons = val;  }
    G4bool GetWritePhotons()          { return fWritePhotons; }           
   
    void SetVerbosity(G4int val)     { fVerbosity = val; }
    G4int GetVerbosity()             { return fVerbosity; }    
   
    void   SetBoronScintillatorIndex(G4int val)   {  fBoronScintillatorIndex = val; }
    G4int  GetBoronScintillatorIndex()            { return fBoronScintillatorIndex; }
    
    void   SetLiquidArgonIndex(G4int val)         {  fLiquidArgonIndex = val; }
    G4int  GetLiquidArgonIndex()                  { return fLiquidArgonIndex; }
    
    void   SetGaseousArgonIndex(G4int val)        {  fGaseousArgonIndex = val; } 
    G4int  GetGaseousArgonIndex()                 { return fGaseousArgonIndex; } 

    void   SetPMTMaterialIndex(G4int val)         {  fPMT = val; } 
    G4int  GetPMTMaterialIndex()                  { return fPMT; } 
    
    void   SetVetoPMTMaterialIndex(G4int val)     {  fMuPMT = val; } 
    G4int  GetVetoPMTMaterialIndex()              { return fMuPMT; } 
    
    void   SetMuPMTMaterialIndex(G4int val)       {  fVetoPMT = val; } 
    G4int  GetMuPMTMaterialIndex()                { return fVetoPMT; } 

    void   SetWriteDeposits(G4int val)            {  fWriteDeposits = val; } 
    G4int  GetWriteDeposits()                     { return fWriteDeposits; } 
 
    void   SetWriteThermalElectrons(G4int val)    {  fWriteThermalElectrons = val; }	     
    G4int  GetWriteThermalElectrons()             { return fWriteThermalElectrons; }	     

    void   SetWriteDaughters(G4int val)           {  fWriteDaughters = val; } 
    G4int  GetWriteDaughters()                    { return fWriteDaughters; } 

    void SetRDMDecay(G4bool val)                  {  fRDMDecay = val;  }
    G4bool GetRDMDecay()                          { return fRDMDecay; }           
    
    void SetRealPDGMeanLife(G4double val)         { fRealPDGMeanLife  = val; }
    G4double GetRealPDGMeanLife()                 { return fRealPDGMeanLife; }           

    void SetPreAbsTime(G4double val)              { fPreAbsTime  = val; }
    G4double GetPreAbsTime()                      { return fPreAbsTime; }           

    void   SetNDaughters(G4int val)               {  fNDaughters = val; } 
    G4int  GetNDaughters()                        { return fNDaughters; } 
 
    void SetRDMChain(G4bool val)                  { fRDMChain = val;  }
    G4bool GetRDMChain()                          { return fRDMChain; }    

    void SetKillS1S2(G4bool val)                  { fKillS1S2 = val;  }
    G4bool GetKillS1S2()                          { return fKillS1S2; }    

    void SetKillS2(G4bool val)                    { fKillS2 = val;  }
    G4bool GetKillS2()                            { return fKillS2; }    

    void SetKillS1(G4bool val)                    { fKillS1 = val;  }
    G4bool GetKillS1()                            { return fKillS1; }    

    void SetScaleS2(G4double val)                 { fScaleS2 = val;  }
    G4double GetScaleS2()                         { return fScaleS2; }    

    void   SetTMBfraction(G4double val)           { fTMBfraction = val;  } 
    G4double  GetTMBfraction()                    { return fTMBfraction; } 

    void   SetTimeCut(G4double val)               { fTimeCut = val;  } 
    G4double  GetTimeCut()                        { return fTimeCut; } 

    void   SetTrackCounter(G4int val)             { fTrackCounter = val;  } 
    G4int  GetTrackCounter()                      { return fTrackCounter; } 

    void SetEndOfEvent(G4bool val)                { fEndOfEvent = val;  }
    G4bool GetEndOfEvent()                        { return fEndOfEvent; }    

    void SetDSLightTrackSecondaries(G4bool val)   { fDSLightTrackSecondaries = val;  }
    G4bool GetDSLightTrackSecondaries()           { return fDSLightTrackSecondaries; }    

    void SetLArGArBoundaryPosZ(G4double val)      { fLArGArBoundaryPositionZ = val;  }         
    G4double GetLArGArBoundaryPosZ()              { return fLArGArBoundaryPositionZ; }

    void SetDriftField( G4double val )            { fDriftField = val;  }
    G4double GetDriftField()                      { return fDriftField; }

    void SetExtractionField( G4double val )       { fExtractionField = val;  }
    G4double GetExtractionField()                 { return fExtractionField; }

    void SetThomasImelNullField( G4double val )   { fThomasImelNullField = val;  }
    G4double GetThomasImelNullField()             { return fThomasImelNullField; }

    void SetThomasImelEp0( G4double val )         { fThomasImelEp0 = val;  }
    G4double GetThomasImelEp0()                   { return fThomasImelEp0; }

    void SetThomasImelEp1( G4double val )         { fThomasImelEp1 = val;  }
    G4double GetThomasImelEp1()                   { return fThomasImelEp1; }

    void SetDokeBirksNFp1( G4double val )         { fDokeBirksNFp1 = val;  }
    G4double GetDokeBirksNFp1()                   { return fDokeBirksNFp1; }

    void SetDokeBirksNFp3( G4double val )         { fDokeBirksNFp3 = val;  }
    G4double GetDokeBirksNFp3()                   { return fDokeBirksNFp3; }

    void SetDokeBirksEp1( G4double val )         { fDokeBirksEp1 = val;  }
    G4double GetDokeBirksEp1()                   { return fDokeBirksEp1; }

    void SetDokeBirksEp2( G4double val )         { fDokeBirksEp2 = val;  }
    G4double GetDokeBirksEp2()                   { return fDokeBirksEp2; }

    void SetDokeBirksEp3( G4double val )         { fDokeBirksEp3 = val;  }
    G4double GetDokeBirksEp3()                   { return fDokeBirksEp3; }

    G4double *GetITOParameters(G4double, G4double, G4double, G4double, G4double);
    G4double GetGridParameters(G4double, G4double, G4double, G4double, G4double);

    void SetNumberOfHits(G4int val)              { fNumberOfHits = val;  }
    G4int GetNumberOfHits()                      { return fNumberOfHits; }

    void SetScintillator(G4int val)              { fScintillator = val;  }
    G4int GetScintillator()                      { return fScintillator; }

    void SetFastSimulation(G4int val)            { fFastSimulation = val;  }
    G4int GetFastSimulation()                    { return fFastSimulation; }

    void SetSourcePosition(G4ThreeVector val)    { fSourcePosition = val; }
    G4ThreeVector GetSourcePosition()            { return fSourcePosition ; }
 
    void SetIsExternalLArScintillating(G4bool val)   { fIsExternalLArScintillating = val;  }
    G4bool GetIsExternalLArScintillating()           { return fIsExternalLArScintillating; }    
    
    void SetVetoYieldFactor( G4double val )         {  fVetoYieldFactor = val;  }
    G4double GetVetoYieldFactor()                   { return fVetoYieldFactor ; }
   
    void SetTunedS1At200V( G4bool val )               {  fTunedS1At200V = val;  }
    G4bool GetTunedS1At200V()                         { return  fTunedS1At200V; }
    
    void SetHolderRadius( G4double val )              {  fHolderRadius = val;  }
    G4double GetHolderRadius()                        { return  fHolderRadius; }
    
    void SetHolderZ( G4double val )                   {  fHolderZ = val;  }
    G4double GetHolderZ()                             { return  fHolderZ; }
    
    void SetHolderPhi( G4double val )                 {  fHolderPhi = val;  }
    G4double GetHolderPhi()                           { return  fHolderPhi; }
    
    void SetHolderSource( G4bool val )                 {  fHolderSource = val;  }
    G4bool GetHolderSource()                           { return  fHolderSource; }
 
    void Set5KGeometry( G4bool val )                 {  fIs5KGeometry = val;  }
    G4bool Get5KGeometry()                           { return  fIs5KGeometry; }

    void Set20KGeometry( G4bool val )                 {  fIs20KGeometry = val;  }
    G4bool Get20KGeometry()                           { return  fIs20KGeometry; }

  private:
    G4bool            fIsEnDepGenerator; 
    G4int             fOverlap; 	    	
    G4int             fEventCounter;	    	
    G4bool            fExportGDML;	    	
    G4bool            fWritePhotons;	    	
    G4int             fVerbosity;	    	
    G4int             fBoronScintillatorIndex;
    G4int             fLiquidArgonIndex;
    G4int             fGaseousArgonIndex;
    G4int             fWriteThermalElectrons;
    G4int             fWriteDeposits;
    G4int             fWriteDaughters;
    G4int             fPMT    ;
    G4int             fVetoPMT;
    G4int             fMuPMT  ;
    G4bool            fRDMDecay;
    G4double          fRealPDGMeanLife;
    G4double          fPreAbsTime;
    G4int             fNDaughters  ;
    G4bool            fRDMChain;
    G4bool            fKillS1S2 ;
    G4bool            fKillS2 ;
    G4bool            fKillS1 ;
    G4double          fScaleS2 ;
    G4double          fITOParameters[3];
    G4double          fGridParameter;
    vector<double>    fAng;
    vector<double>    fTrans;
    G4double          fTMBfraction;  
    G4double          fTimeCut;  

    G4double          fWL[11];
    G4double          fRI[11];
    G4double          fEC[11];   
    G4int             fTrackCounter ;  
    G4bool            fEndOfEvent ; 
    G4bool            fDSLightTrackSecondaries;
    G4double          fLArGArBoundaryPositionZ;  // z coordinate of the LAr-GAr interface - detector configuration dependent
    G4double          fDriftField;               // Electric field (V/cm)
    G4double          fExtractionField;          // Electric field (V/cm)

    G4double          fThomasImelNullField;      // Thomas Imel parameter at null field
    G4double          fThomasImelEp0;            // Thomas Imel parameter at non-null field = p0 * E**p1
    G4double          fThomasImelEp1;

    G4double          fDokeBirksNFp1;
    G4double          fDokeBirksNFp3;
    G4double          fDokeBirksEp1;
    G4double          fDokeBirksEp2;
    G4double          fDokeBirksEp3;

    G4int             fNumberOfHits;
    G4int             fScintillator;
    G4int             fFastSimulation;
    
    G4ThreeVector     fSourcePosition ;
    G4bool            fIsExternalLArScintillating;
    G4double          fVetoYieldFactor ;
    
    G4bool            fTunedS1At200V;
    
    G4double          fHolderRadius; 
    G4double          fHolderZ;
    G4double          fHolderPhi;
    G4bool            fHolderSource;

    G4bool            fIs5KGeometry;
    G4bool            fIs20KGeometry;
};

#endif
/*
 * $Log: DSStorage.hh,v $
 * Revision 1.11  2015/04/23 14:04:07  pagnes
 * DS20K geometry added (config 10)
 *
 * Revision 1.10  2015/01/14 16:58:42  dfranco
 * added source vial and correspondent commands by Laura and Bianca; manual updated
 *
 * Revision 1.9  2014/12/22 14:40:49  dfranco
 * added the option to activate the recombination probability at 200 V/cm (/ds/physics/tunedS1); this option is by default true; selecting a specific drift field automatically switch off the tunedS1 option
 *
 * Revision 1.8  2014/11/21 10:19:06  dfranco
 * added a command to scale the veto scintillation yield factor and fixed the visible energy variable in the veto
 *
 * Revision 1.7  2014/11/20 15:32:12  dfranco
 * added a command to remove scintillation process from liquid argon between TPC and cryostat
 *
 * Revision 1.6  2014/11/06 17:39:51  dfranco
 * added source and fixed a bug: removed gaseous nitrogen from the neutron veto
 *
 * Revision 1.5  2014/10/13 18:43:46  swesterd
 * fixed veto PMT positions
 *
 * Revision 1.4  2014/07/23 14:52:38  pagnes
 * write thermal e- and kill S1 commands added
 *
 * Revision 1.3  2014/07/16 08:23:12  pagnes
 * QE scaling to 1.0 added (/ds/manager/fast_simulation xxx)
 *
 * Revision 1.2  2014/05/07 14:27:26  dfranco
 * fixed some bugs and added GdScintillator
 *
 * Revision 1.1  2014/05/07 12:20:55  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.16  2014/04/11 12:33:35  pagnes
 * command to set TMB/PC ratio inside veto added
 *
 * Revision 1.15  2014/04/11 10:20:44  perassos
 * Added generation in materials
 *
 * Revision 1.14  2014/03/19 16:37:36  dfranco
 * update external configuration reader for optics tuning
 *
 * Revision 1.13  2014/03/11 09:56:25  meregaglia
 * Added generator starting from energy deposits
 *
 * Revision 1.12  2014/01/29 13:13:42  perassos
 * Update of the electric field handling and of the Nuclear Recoils generator
 *
 * Revision 1.11  2014/01/07 14:10:32  perassos
 * Added the commands to set in the macfile the electric field and the Thomas-Imel parameters
 *
 * Revision 1.10  2013/11/19 10:33:25  perassos
 * Added methods to handle the electric field and the liquid/gas interface z coordinate
 *
 * Revision 1.9  2013/08/01 14:28:52  dfranco
 * added energy loss data from SRIM/TRIM
 *
 * Revision 1.8  2013/07/24 09:48:58  dfranco
 * Added S1 and S2 equivalent energies, the material index for each deposit, the command killS1S2 to kill photons and electrons generated by DSLight (after storing the equivalent energies)
 *
 * Revision 1.7  2013/06/11 22:48:36  dfranco
 * Added ITO optical boundary. The Sernelius function is defined in DSStorage, and called by G4OpBoundaryProcess. New ITO bool variable added to DSDetectorDS50.cc as surface property (G4MaterialPropertyTable)
 *
 * Revision 1.6  2013/06/10 14:15:41  dfranco
 * Added two commands: /ds/physics/killS2 and /ds/physics/scaleS2 to kill or scale the S2 light
 *
 * Revision 1.5  2013/04/03 10:14:24  dfranco
 * Fixed bugs with RDM and RDMChain staking actions. The logic of Geant4 is changed. Different excited states of a nucleus correspond to new particles (trackID). Code adapted.
 *
 * Revision 1.4  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
