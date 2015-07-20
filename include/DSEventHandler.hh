#ifndef _DSEVENTHANDLER_HH
#define _DSEVENTHANDLER_HH 1

#include "DSEventStructure.hh"
#include "G4ThreeVector.hh"

#include "DSLogger.hh"
#include "G4String.hh"
#include <iostream>
#include <fstream>

using namespace std;

class DSEventHandler  {
  private:
    DSEventHandler();
  public:
    static DSEventHandler* Get();
    
    virtual ~DSEventHandler() {}
    
     
  private:
  
    static DSEventHandler *me;
    
  
  public:
  
    HeaderStructure& GetHeader()        { return theHeaderStructure               ; }
    
    void   SetEvents(int val)           { theHeaderStructure.Events          = val; }
    void   SetRun(int val)              { theHeaderStructure.Run             = val; }
    void   SetParticle(int val)         { theHeaderStructure.PDG             = val; }
    void   SetRate(float val)           { theHeaderStructure.Rate            = val; }  
    void   SetLArIndex(int val)         { theHeaderStructure.LArIndex         = val; }
    void   SetScintillatorIndex(int val){ theHeaderStructure.ScintillatorIndex= val; }
    void   SetDetectorFlag(int val)     { theHeaderStructure.DetectorFlag    = val; }
    

    int     GetEvents()	                { return theHeaderStructure.Events	  ; }	  
    int     GetRun() 	                { return theHeaderStructure.Run	          ; }
    int     GetParticle()	        { return theHeaderStructure.PDG	          ; }
    int     GetScintillatorIndex()      { return theHeaderStructure.ScintillatorIndex     ; }
    int     GetLArIndex()	        { return theHeaderStructure.LArIndex     ; }
    float   GetRate()  	                { return theHeaderStructure.Rate	  ; }	  
    int     GetDetectorFlag()	        { return theHeaderStructure.DetectorFlag  ; }	  
    
    
    
    EventStructure&             GetEvent()               { return theEventStructure             ; }   
    EventStructureDiskFormat&   GetEventDiskFormat()     { return theEventStructureDiskFormat   ; }   
    
    void    SetEventID(int val)                  { theEventStructure.EventID	    = val; }
    void    SetPDG(int val)	                 { theEventStructure.PDG	    = val; }
    void    SetTime(double val)                  { theEventStructure.Time	    = val; }
    void    SetEnergy(float val)                 { theEventStructure.Energy	    = val; }
    void    SetS1Energy(float val)               { theEventStructure.S1Energy       = val; }
    void    SetS2Energy(float val)               { theEventStructure.S2Energy       = val; }
    void    SetVetoVisEnergy(float val)          { theEventStructure.VetoVisEnergy  = val; }
    void    SetMuVisEnergy(float val)            { theEventStructure.MuVisEnergy    = val; }
    void    SetTPCDepEnergy(float val)           { theEventStructure.TPCDepEnergy   = val; }
    void    SetVetoDepEnergy(float val)          { theEventStructure.VetoDepEnergy  = val; }
    void    SetMuDepEnergy(float val)            { theEventStructure.MuDepEnergy    = val; }
    void    SetPosition(G4ThreeVector   );       
    void    SetCenterOfMass(G4ThreeVector );     
    void    SetDirection(G4ThreeVector  ) ;     
    void    SetNDaughters(int val)               { theEventStructure.NDaughters     = val; }
    void    SetNDeposits(int val)                { theEventStructure.NDeposits      = val; }
    void    SetNPE(int val)	                 { theEventStructure.NPE            = val; }
    void    SetMuNPE(int val)	                 { theEventStructure.MuNPE          = val; }
    void    SetVetoNPE(int val)	                 { theEventStructure.VetoNPE        = val; }
    void    SetNPH(int val)                      { theEventStructure.NPH            = val; }
    void    SetNUsers(int val)                   { theEventStructure.NUsers	    = val; }
    void    SetDaughters()	                 { theEventStructure.theDaughters.push_back(theDaughterStructure)                  ; }
    void    SetPhotoElectrons()	                 { theEventStructure.thePhotoElectrons.push_back(thePhotoElectronStructure)        ; }
    void    SetVetoPhotoElectrons()	         { theEventStructure.theVetoPhotoElectrons.push_back(theVetoPhotoElectronStructure); }
    void    SetMuPhotoElectrons()	         { theEventStructure.theMuPhotoElectrons.push_back(theMuPhotoElectronStructure)    ; }
    void    SetPhotons()	                 { theEventStructure.thePhotons.push_back(thePhotonStructure)                      ; }
    void    SetDeposits()	                 { theEventStructure.theDeposits.push_back(theDepositStructure)                    ; }
    void    SetUsers()	                         { theEventStructure.theUsers.push_back(theUserStructure)                          ; }
    
    int            GetEventID()	     	         { return theEventStructure.EventID	                    ; }   
    int            GetPDG()		     	 { return theEventStructure.PDG	                            ; }   
    double         GetTime()   		         { return theEventStructure.Time	                    ; }   
    float          GetEnergy() 	     	         { return theEventStructure.Energy	                    ; }   
    float          GetS1Energy() 	     	 { return theEventStructure.S1Energy	                    ; }   
    float          GetS2Energy() 	     	 { return theEventStructure.S2Energy	                    ; }   
    float          GetVetoVisEnergy() 	     	 { return theEventStructure.VetoVisEnergy	            ; }   
    float          GetMuVisEnergy() 	     	 { return theEventStructure.MuVisEnergy	                    ; }   
    float          GetTPCDepEnergy() 	     	 { return theEventStructure.TPCDepEnergy		    ; }   
    float          GetVetoDepEnergy() 	     	 { return theEventStructure.VetoDepEnergy		    ; }   
    float          GetMuDepEnergy() 	     	 { return theEventStructure.MuDepEnergy		            ; }   
    G4ThreeVector  GetPosition()		 { return CopyArrayToVector(theEventStructure.Position)     ; }		
    G4ThreeVector  GetCenterOfMass()		 { return CopyArrayToVector(theEventStructure.CenterOfMass) ; }		   
    G4ThreeVector  GetDirection()		 { return CopyArrayToVector(theEventStructure.Direction)    ; }		   
    int            GetNDaughters()	         { return theEventStructure.NDaughters                      ; }   
    int            GetNDeposits()	         { return theEventStructure.NDeposits                       ; }   
    int            GetNPE()		         { return theEventStructure.NPE	                            ; }   
    int            GetVetoNPE()		         { return theEventStructure.VetoNPE	                    ; }   
    int            GetMuNPE()		         { return theEventStructure.MuNPE	                    ; }   
    int            GetNPH()	                 { return theEventStructure.NPH	                            ; }     
    int            GetNUsers()	        	 { return theEventStructure.NUsers	                    ; }     
    
    vector<DaughterStructure>&          GetVDaughters()          { return theEventStructure.theDaughters            ; }
    vector<PhotoElectronStructure>&     GetVPhotoElectrons()     { return theEventStructure.thePhotoElectrons       ; }
    vector<PhotoElectronStructure>&     GetVMuPhotoElectrons()   { return theEventStructure.theMuPhotoElectrons     ; }
    vector<PhotoElectronStructure>&     GetVVetoPhotoElectrons() { return theEventStructure.theVetoPhotoElectrons   ; }
    vector<PhotonStructure>&            GetVPhotons()            { return theEventStructure.thePhotons              ; }    
    vector<DepositStructure>&           GetVDeposits()           { return theEventStructure.theDeposits             ; }    
    vector<UserStructure>&              GetVUsers()              { return theEventStructure.theUsers	            ; }
         
    
    DaughterStructure& GetDaughters()                   { return theDaughterStructure		; }
    
    void   SetDId(int val)                              { theDaughterStructure.Id	 = val  ; }
  //  void   SetDPID(int val)                             { theDaughterStructure.PID	 = val  ; }
    void   SetDTrackID(int val)                             { theDaughterStructure.TrackID	 = val  ; }
    void   SetDParentTrackID(int val)                             { theDaughterStructure.ParentTrackID	 = val  ; }

    void   SetDPDG(int val)                             { theDaughterStructure.PDG	 = val  ; }
    void   SetDProcess(int val)                         { theDaughterStructure.Process	 = val  ; }
    void   SetDTime(double val)                         { theDaughterStructure.Time	 = val  ; }
    void   SetDEnergy(float val)                        { theDaughterStructure.Energy	 = val  ; }
    void   SetDPosition(G4ThreeVector );       
    void   SetDDirection(G4ThreeVector );      

    int    GetDId()                                     { return theDaughterStructure.Id	; } 
  //    int    GetDPID()                                     { return theDaughterStructure.PID	; } 
    int    GetDTrackID()                                     { return theDaughterStructure.TrackID	; } 
    int    GetDParentTrackID()                                     { return theDaughterStructure.ParentTrackID	; } 

    int    GetDPDG()                                    { return theDaughterStructure.PDG	; } 
    int    GetDProcess()                                { return theDaughterStructure.Process	; } 
    double GetDTime()                                   { return theDaughterStructure.Time	; } 
    float  GetDEnergy()                                 { return theDaughterStructure.Energy	; } 
    G4ThreeVector   GetDDirection()                     { return CopyArrayToVector(theDaughterStructure.Position)  ; } 
    G4ThreeVector   GetDPosition()                      { return CopyArrayToVector(theDaughterStructure.Direction) ; } 
     
  
    PhotoElectronStructure&   GetPhotoElectrons()       { return thePhotoElectronStructure      ; }
    
    void   SetPhotoElectronPMT(int val)                 { thePhotoElectronStructure.PMT    = val; }
    void   SetPhotoElectronTime(double val)             { thePhotoElectronStructure.Time   = val; }

    int    GetPhotoElectronPMT()                        { return thePhotoElectronStructure.PMT    ; }
    double GetPhotoElectronTime()                       { return thePhotoElectronStructure.Time	  ; }
    


    PhotoElectronStructure&   GetMuPhotoElectrons()     { return theMuPhotoElectronStructure      ; }
    
    void   SetMuPhotoElectronPMT(int val)               { theMuPhotoElectronStructure.PMT    = val; }			       
    void   SetMuPhotoElectronTime(double val)           { theMuPhotoElectronStructure.Time   = val; }			       

    int    GetMuPhotoElectronPMT()                      { return theMuPhotoElectronStructure.PMT    ; } 		       
    double GetMuPhotoElectronTime()                     { return theMuPhotoElectronStructure.Time	; }		       
    

    PhotoElectronStructure&   GetVetoPhotoElectrons()   { return theVetoPhotoElectronStructure      ; } 		       
    
    void   SetVetoPhotoElectronPMT(int val)             { theVetoPhotoElectronStructure.PMT    = val; } 			    
    void   SetVetoPhotoElectronTime(double val)         { theVetoPhotoElectronStructure.Time   = val; } 			    

    int    GetVetoPhotoElectronPMT()                    { return theVetoPhotoElectronStructure.PMT    ; }			    
    double GetVetoPhotoElectronTime()                   { return theVetoPhotoElectronStructure.Time   ; }			    


    PhotonStructure&   GetPhotons()                     { return thePhotonStructure	  ;    }		       
    
    void   SetPhotonVolumeID(int val)                   { thePhotonStructure.VolumeID	= val; }		       
    void   SetPhotonTime(double val)                    { thePhotonStructure.Time	= val; }		       
    void   SetPhotonPID(int val)                        { thePhotonStructure.PID        = val;  }		       
    void   SetPhotonWavelength(float val)               { thePhotonStructure.Wavelength = val;  }		       
    void   SetPhotonPosition(G4ThreeVector);

    int    GetPhotonVolumeID()                          { return thePhotonStructure.VolumeID	; }		       
    double GetPhotonTime()                              { return thePhotonStructure.Time	; }		       
    int    GetPhotonPID()                               { return thePhotonStructure.PID 	; }		       
    float  GetPhotonWavelength()                        { return thePhotonStructure.Wavelength  ;  }		       
    G4ThreeVector  GetPhotonPosition()                  { return CopyArrayToVector(thePhotonStructure.Position) ; }    
         
    
    DepositStructure&  GetDeposits()                    { return theDepositStructure		; }

    void   SetDepPID(int val)                           { theDepositStructure.PID          = val; }    

    void   SetDepTrack(int val)                         { theDepositStructure.Track        = val; }    
    void   SetDepParentTrack(int val)                         { theDepositStructure.ParentTrack        = val; }    

    void   SetDepVolume(int val)                        { theDepositStructure.Volume       = val; }    
    void   SetDepTotalEnergy(float val)                      { theDepositStructure.TotalEnergy	   = val; }
    void   SetDepKineticEnergy(float val)                      { theDepositStructure.KineticEnergy	   = val; }
    void   SetDepEnergy(float val)                      { theDepositStructure.Energy	   = val; }
    void   SetDepStep(float val)                        { theDepositStructure.Step	   = val; }
    void   SetDepTime(double val)                       { theDepositStructure.Time	   = val; }
    void   SetDepPosition(G4ThreeVector) ;
    
    int    GetDepPID()                                   { return theDepositStructure.PID       ; }    

    int    GetDepTrack()                                 { return theDepositStructure.Track     ; }    
    int    GetDepParentTrack()                                 { return theDepositStructure.ParentTrack     ; }    
    int    GetDepVolume()                                { return theDepositStructure.Volume    ; }    
    float  GetDepTotalEnergy()                                { return theDepositStructure.TotalEnergy    ; }
    float  GetDepKineticEnergy()                                { return theDepositStructure.KineticEnergy    ; }
    float  GetDepEnergy()                                { return theDepositStructure.Energy    ; }
    float  GetDepStep()                                  { return theDepositStructure.Step      ; }
    double GetDepTime()                                  { return theDepositStructure.Time      ; }
    G4ThreeVector  GetDepPosition()                      { return CopyArrayToVector(theDepositStructure.Position)  ; }
    


    UserStructure&  GetUsers()                          { return theUserStructure	      ; }

    void  SetUserInt1(int val) 	   	                { theUserStructure.UserInt1 = val     ; }   
    void  SetUserInt2(int val) 	   	                { theUserStructure.UserInt2 = val     ; }   
    void  SetUserFloat1(float val) 	                { theUserStructure.UserFloat1 = val   ; }   
    void  SetUserFloat2(float val) 	                { theUserStructure.UserFloat2 = val   ; }   
    void  SetUserDouble(double val)	                { theUserStructure.UserDouble = val   ; }   
 
    int     GetUserInt1()                               { return  theUserStructure.UserInt1   ; }     
    int     GetUserInt2()                               { return  theUserStructure.UserInt2   ; }     
    float   GetUserFloat1()                             { return  theUserStructure.UserFloat1 ; }     
    float   GetUserFloat2()                             { return  theUserStructure.UserFloat2 ; }     
    double  GetUserDouble()                             { return  theUserStructure.UserDouble ; }     

    
    
    void SetEventStructureDiskFormat(EventStructureDiskFormat val)   { theEventStructureDiskFormat   = val;}
    void SetDaughter(DaughterStructure val)                          { theDaughterStructure          = val;}   
    void SetPhotoElectron(PhotoElectronStructure val)                { thePhotoElectronStructure     = val;} 
    void SetMuPhotoElectron(PhotoElectronStructure val)              { theMuPhotoElectronStructure   = val;} 
    void SetVetoPhotoElectron(PhotoElectronStructure val)            { theVetoPhotoElectronStructure = val;} 
    void SetPhoton(PhotonStructure val)                              { thePhotonStructure            = val;} 
    void SetDeposit(DepositStructure val)                            { theDepositStructure           = val;}    
    void SetUser(UserStructure val)                                  { theUserStructure              = val;} 
   
    

    void ClearHeader();
    void ClearEvent();
    void ClearPhotoElectron();
    void ClearMuPhotoElectron();
    void ClearVetoPhotoElectron();
    void ClearPhoton();
    void ClearDaughter();
    void ClearDeposit();
    void ClearUser();
    void ClearAll();

    void DumpHeader();
    void DumpEvent();
    void DumpPhotoElectron();
    void DumpMuPhotoElectron();
    void DumpVetoPhotoElectron();
    void DumpPhoton();
    void DumpDaughter();
    void DumpDeposit();
    void DumpUser();
    void DumpAll();

    G4ThreeVector                 CopyArrayToVector(G4float*);
   
  private:
    EventStructureDiskFormat      theEventStructureDiskFormat;
    EventStructure                theEventStructure;
    PhotoElectronStructure        thePhotoElectronStructure;
    PhotoElectronStructure        theMuPhotoElectronStructure;
    PhotoElectronStructure        theVetoPhotoElectronStructure;
    PhotonStructure               thePhotonStructure;
    DaughterStructure             theDaughterStructure;
    DepositStructure              theDepositStructure;
    HeaderStructure               theHeaderStructure;
    UserStructure                 theUserStructure;

    G4int                         fVerbosity;
};

#endif
/*
 * $Log: DSEventHandler.hh,v $
 * Revision 1.4  2014/11/13 16:47:09  dfranco
 * removed variables which were creating conflicts with the previous version of g4ds10
 *
 * Revision 1.3  2014/10/13 18:43:46  swesterd
 * fixed veto PMT positions
 *
 * Revision 1.2  2014/05/08 10:59:21  pagnes
 * Scintillator Index added in binary header
 *
 * Revision 1.1  2014/05/07 12:20:51  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.8  2013/10/22 14:03:07  swesterd
 * getvphotoelectrons was returning veto photoelectrons instead of tpc photoelectrons
 *
 * Revision 1.7  2013/08/20 03:25:49  swesterd
 * added G2 TPC geoemtry (not complete) and added monoenergetic energy distribution to generator
 *
 * Revision 1.6  2013/08/06 13:58:18  dfranco
 * Added 3 variables to store the deposited energy in LAr, scintillator, and water. The last two are not yet implemented. g4rooter has been updated with 3 new variables: tpcene, vetoene, and muene
 *
 * Revision 1.5  2013/07/24 09:48:57  dfranco
 * Added S1 and S2 equivalent energies, the material index for each deposit, the command killS1S2 to kill photons and electrons generated by DSLight (after storing the equivalent energies)
 *
 * Revision 1.4  2013/04/04 09:04:17  dfranco
 * added step length info to the deposit structure
 *
 * Revision 1.3  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
