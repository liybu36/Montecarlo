#include "DSEventHandler.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "DSLogger.hh"
#include "DSParameters.hh"
#include "DSStorage.hh"
#include "DSTrackingAction.hh"
#include "G4TrackingManager.hh"         
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4TrackStatus.hh"
#include <vector>

using namespace std;



DSTrackingAction::DSTrackingAction() { ; }
/*
void DSTrackingAction::PreUserTrackingAction(const G4Track* aTrack) {
  
  if(aTrack->GetDefinition()->GetPDGEncoding() != 50)
  cout << "Track " << aTrack->GetDefinition()->GetParticleName() << " "
       << aTrack->GetVertexKineticEnergy ()/MeV << " "
 	<< aTrack->GetTrackID() << " " 
 	<< aTrack->GetParentID() << " " 
       <<endl ;
    
  
  if(   aTrack->GetDefinition()->GetPDGEncoding() == 5 
     || aTrack->GetDefinition()->GetPDGEncoding() == 6 
     || aTrack->GetDefinition()->GetPDGEncoding() == 27 
     || aTrack->GetDefinition()->GetPDGEncoding() == 28 
     || aTrack->GetDefinition()->GetPDGEncoding() == 43 
     || aTrack->GetDefinition()->GetPDGEncoding() == 44   
  ) aTrack->SetTrackStatus(fStopAndKill );
  
  
  
  //cout << fKillTrackAndSecondaries << endl ;
  //if(aTrack->GetGlobalTime () > 1.*ms) aTrack->SetTrackStatus(  3);
  
      
  if(     DSStorage::Get()->GetWriteDaughters()  
       && DSStorage::Get()->GetRDMDecay() == 0 
       && aTrack->GetParentID()  <= DSStorage::Get()->GetNDaughters() // default = 1
       && aTrack->GetTrackID()   > 1
       && DSStorage::Get()->GetIsEnDepGenerator() == false 
       && aTrack->GetDefinition()->GetPDGEncoding() != 50 ) {

    int index = aTrack->GetCreatorProcess ()->GetProcessType ()*1000 + aTrack->GetCreatorProcess ()->GetProcessSubType ();
    DSEventHandler::Get()->SetDId(int(DSEventHandler::Get()->GetVDaughters().size()));
    //   DSEventHandler::Get()->SetDPID(aTrack->GetTrackID());
    DSEventHandler::Get()->SetDTrackID(aTrack->GetTrackID());
    DSEventHandler::Get()->SetDParentTrackID(aTrack->GetParentID());

    DSEventHandler::Get()->SetDPDG(aTrack->GetDefinition()->GetPDGEncoding());
    DSEventHandler::Get()->SetDProcess(index);
    DSEventHandler::Get()->SetDTime(aTrack->GetGlobalTime()/ns);
    DSEventHandler::Get()->SetDEnergy(aTrack->GetVertexKineticEnergy ()/MeV);
    DSEventHandler::Get()->SetDPosition(aTrack->GetVertexPosition()/m);
    DSEventHandler::Get()->SetDDirection(aTrack->GetVertexMomentumDirection() );
    DSEventHandler::Get()->SetDaughters();


    if(aTrack->GetDefinition()->GetPDGEncoding() != 50)
          DSLog(development) << " Save daughters: "
		<< aTrack->GetDefinition()->GetParticleName() << " " 
		<< aTrack->GetTrackID() << " " 
 		<< aTrack->GetParentID() << " " 
	       << aTrack->GetVertexKineticEnergy ()/MeV << " " 
	       << aTrack->GetGlobalTime()/ns
               << endlog ;
  }
  
}
*/
/////////////////////////////////////////////////////////////////////////



void DSTrackingAction::PostUserTrackingAction(const G4Track* aTrack) {

   if(   aTrack->GetTrackStatus() == fStopAndKill   
      && aTrack->GetDefinition()->GetPDGEncoding() == 50 
      && DSLogger::GetSeverity() == DSLogger::development )
     DSLog(development) << " EoT " << aTrack->GetDefinition()->GetParticleName() << " "
	  << aTrack->GetStep()->GetPostStepPoint()->GetPosition()/cm<< " "
	<< aTrack->GetVolume()->GetName()
	<<endlog ;
   
  // save tpc photoelectrons  
  if(  DSStorage::Get()->GetPMTMaterialIndex() == (G4int) aTrack->GetVolume()->GetLogicalVolume()->GetMaterial()->GetIndex() 
    && aTrack->GetDefinition()->GetPDGEncoding() == 50
    && aTrack->GetTrackStatus() == fStopAndKill ) {
   
    if (DSStorage::Get()->Get20KGeometry()) {
      DSEventHandler::Get()->SetPhotonVolumeID(aTrack->GetVolume()->GetLogicalVolume()->GetMaterial()->GetIndex ()); 
      //DSEventHandler::Get()->SetPhotonVolumeID(aTrack->GetVolume()->GetCopyNo ()  );
      DSEventHandler::Get()->SetPhotonTime(aTrack->GetGlobalTime()/ns); 
      DSEventHandler::Get()->SetPhotonPID(aTrack->GetParentID());  
      DSEventHandler::Get()->SetPhotonWavelength(h_Planck*c_light*1E12/(aTrack->GetKineticEnergy()/eV));  
      DSEventHandler::Get()->SetPhotonPosition(aTrack->GetStep()->GetPostStepPoint()->GetPosition());    
      DSEventHandler::Get()->SetPhotons();
    } else {
      G4String pmtname =  aTrack->GetVolume()->GetName();
      if(!pmtname.find("TPMT_")) {
	// Get PMT Q.E.
	//vector <G4double> myPmtQEWL = DSParameters::Get()->GetPmtQeWl();
	//vector <G4double> myPmtQE   = DSParameters::Get()->GetPmtQe();

	// Computer the photon wavelength
	//G4double myPhotonWL = nm*h_Planck*c_light*1E12/(aTrack->GetKineticEnergy()/eV);
	//G4double myQE = 0.;
	//myQE = interpolate(&myPmtQE[0],&myPmtQEWL[0],myPhotonWL,(int)myPmtQE.size());
	//DSLog(development) << "WL = " << myPhotonWL << "tQE = " << myQE << endlog;

	// apply the QE cut
	if(G4UniformRand() < DSParameters::Get()->GetTPCQE(h_Planck*c_light*1E12/(aTrack->GetKineticEnergy()/eV))) {      
     //   if(G4UniformRand() < myQE) {
          pmtname.erase(0,5);
          DSEventHandler::Get()->SetPhotoElectronPMT(atoi(pmtname.c_str()));
          DSEventHandler::Get()->SetPhotoElectronTime(aTrack->GetGlobalTime()/ns);
          DSEventHandler::Get()->SetPhotoElectrons();

          if(DSStorage::Get()->GetVerbosity() > 3) {
            DSLog(development) << "TPC PhotoElectron: " 
        	 << "PMT " << DSEventHandler::Get()->GetPhotoElectronPMT()
        	 << "T   " << DSEventHandler::Get()->GetPhotoElectronTime() << " ns " 
        	 << aTrack->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName() << " "  
        	 << aTrack->GetVolume()->GetName() << " "  
        	 << aTrack->GetVolume()->GetCopyNo () << " "  
        	 << aTrack->GetVolume()->GetLogicalVolume()->GetName() << " "  
        	 << aTrack->GetStep()->GetPostStepPoint()->GetPosition()/m << " " 
        	 << aTrack->GetStep()->GetPostStepPoint()->GetPosition().mag()/m << " " 
        	 << endlog ;
	  }     
        }
      }
    }   
  }
  
  
  // save veto photoelectrons  
  if( DSStorage::Get()->GetVetoPMTMaterialIndex() == (G4int) aTrack->GetVolume()->GetLogicalVolume()->GetMaterial()->GetIndex() 
    && aTrack->GetDefinition()->GetPDGEncoding() == 50
      && aTrack->GetTrackStatus() == fStopAndKill ) {
    DSLog(development) << "in veto PMT" << endlog;
    G4String pmtname =  aTrack->GetVolume()->GetName();
    if(!pmtname.find("VPMT_")) {
      
      // Get PMT Q.E.
      vector <G4double> myVPmtQEWL = DSParameters::Get()->GetVPmtQeWl();
      vector <G4double> myVPmtQE   = DSParameters::Get()->GetVPmtQe();

      // Compute the photon wavelength
      G4double myPhotonWL = nm*h_Planck*c_light*1E12/(aTrack->GetKineticEnergy()/eV);
      G4double myQE = 0.;
      myQE = interpolate(&myVPmtQE[0],&myVPmtQEWL[0],myPhotonWL,(int)myVPmtQE.size());
      DSLog(development) << "WL = " << myPhotonWL << "\tQE = " << myQE << endlog;
      // apply the QE cut
      if(G4UniformRand() <  myQE) {
	pmtname.erase(0,5);
	DSEventHandler::Get()->SetVetoPhotoElectronPMT(atoi(pmtname.c_str()));
	DSEventHandler::Get()->SetVetoPhotoElectronTime(aTrack->GetGlobalTime()/ns);
	DSEventHandler::Get()->SetVetoPhotoElectrons();
	
	if(DSStorage::Get()->GetVerbosity() > 3) {
	  DSLog(development) << "Veto PhotoElectron: " 
			   << "PMT " << DSEventHandler::Get()->GetVetoPhotoElectronPMT()
			   << "T   " << DSEventHandler::Get()->GetVetoPhotoElectronTime() << " ns " 
			   << aTrack->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName() << " "  
			   << aTrack->GetVolume()->GetName() << " "  
			   << aTrack->GetVolume()->GetCopyNo () << " "  
			   << aTrack->GetVolume()->GetLogicalVolume()->GetName() << " "  
			   << aTrack->GetStep()->GetPostStepPoint()->GetPosition()/m << " " 
			   << aTrack->GetStep()->GetPostStepPoint()->GetPosition().mag()/m << " " 
			   << endlog ;
	}
      }
    }   
  }

  // save water tank photoelectrons  
  if(  DSStorage::Get()->GetMuPMTMaterialIndex() == (G4int) aTrack->GetVolume()->GetLogicalVolume()->GetMaterial()->GetIndex() 
    && aTrack->GetDefinition()->GetPDGEncoding() == 50
    && aTrack->GetTrackStatus() == fStopAndKill ) {
    G4String pmtname =  aTrack->GetVolume()->GetName();
    if(!pmtname.find("WPMT_")) {
      pmtname.erase(0,5);
      DSEventHandler::Get()->SetMuPhotoElectronPMT(atoi(pmtname.c_str()));
      DSEventHandler::Get()->SetMuPhotoElectronTime(aTrack->GetGlobalTime()/ns);
      DSEventHandler::Get()->SetMuPhotoElectrons();
    
    
      if(DSStorage::Get()->GetVerbosity() > 3) {
        DSLog(development) << "WT PhotoElectron: " 
             << "PMT " << DSEventHandler::Get()->GetMuPhotoElectronPMT()
             << "T   " << DSEventHandler::Get()->GetMuPhotoElectronTime() << " ns " 
             << aTrack->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName() << " "  
             << aTrack->GetVolume()->GetName() << " "  
             << aTrack->GetVolume()->GetCopyNo () << " "  
             << aTrack->GetVolume()->GetLogicalVolume()->GetName() << " "  
             << aTrack->GetStep()->GetPostStepPoint()->GetPosition()/m << " " 
             << aTrack->GetStep()->GetPostStepPoint()->GetPosition().mag()/m << " " 
             << endlog ;
      }
    }
    
  }
  
  
  // save photons
  if(   DSStorage::Get()->GetWritePhotons()
    && aTrack->GetTrackStatus() == fStopAndKill
    && aTrack->GetDefinition()->GetPDGEncoding() == 50
    && !DSStorage::Get()->Get20KGeometry()     ) {
    if(DSStorage::Get()->GetVerbosity() > 4) {
      DSLog(development) << "Photon: " << aTrack->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName() << " "  
           << aTrack->GetVolume()->GetLogicalVolume()->GetName() << " "  
           << aTrack->GetStep()->GetPostStepPoint()->GetPosition()/m << " " 
           << aTrack->GetStep()->GetPostStepPoint()->GetPosition().mag()/m << " " 
           << endlog ;
    }
    DSEventHandler::Get()->SetPhotonVolumeID(aTrack->GetVolume()->GetLogicalVolume()->GetMaterial()->GetIndex ()); 
    //DSEventHandler::Get()->SetPhotonVolumeID(aTrack->GetVolume()->GetCopyNo ()  );
    DSEventHandler::Get()->SetPhotonTime(aTrack->GetGlobalTime()/ns); 
    DSEventHandler::Get()->SetPhotonPID(aTrack->GetParentID());  
    DSEventHandler::Get()->SetPhotonWavelength(h_Planck*c_light*1E12/(aTrack->GetKineticEnergy()/eV));  
    DSEventHandler::Get()->SetPhotonPosition(aTrack->GetStep()->GetPostStepPoint()->GetPosition());    
    
    
    DSEventHandler::Get()->SetPhotons();
    
  }
}


//binary search
int DSTrackingAction::binSearch(double* list, double x, int minIndex, int maxIndex) {
  if(x < list[minIndex+1] && x >= list[maxIndex-1])
    return minIndex;
  if(x < list[(minIndex+maxIndex)/2])
    return binSearch(list, x, minIndex, (minIndex+maxIndex)/2);
  return binSearch(list, x, (minIndex+maxIndex)/2, maxIndex);
}


int DSTrackingAction::search(double* list, double x, int listSize) {
  if(x <= list[0])
    return 0;
  if(x >= list[listSize-1])
    return listSize-1;
  return binSearch(list, x, 0, listSize-1);
}

//quadratically interpolate a point in an array
double DSTrackingAction::interpolate(double* ys, double *xs, double x, int size) {
  // Find the greatest element <= x
  int myIndex = search(xs,x,size);
  for(int i = 0; i < size; i++)
  // Check if it is below or above the range. If so, return endpoint
  if(myIndex == 0)
    return ys[0];
  if(myIndex == size-1)
    return ys[size-1];
  // Check if x falls on a known point and no interpolation is needed
  if(xs[myIndex] == x)
    return ys[myIndex];

  double denom1 = (xs[myIndex-1]-xs[myIndex])*(xs[myIndex-1]-xs[myIndex+1]);
  double denom2 = (xs[myIndex-1]-xs[myIndex])*(xs[myIndex]-xs[myIndex+1]);
  double denom3 = (xs[myIndex-1]-xs[myIndex+1])*(xs[myIndex]-xs[myIndex+1]);

  double b = ys[myIndex-1]/denom1 - ys[myIndex]/denom2 + ys[myIndex+1]/denom3;
  double c = -(xs[myIndex]+xs[myIndex+1])*ys[myIndex-1]/denom1 + (xs[myIndex-1]+xs[myIndex+1])*ys[myIndex]/denom2 - (xs[myIndex-1]+xs[myIndex])*ys[myIndex+1]/denom3;
  double d = xs[myIndex]*xs[myIndex+1]*ys[myIndex-1]/denom1 - xs[myIndex-1]*xs[myIndex+1]*ys[myIndex]/denom2 + xs[myIndex-1]*xs[myIndex]*ys[myIndex+1]/denom3;

  return b*x*x+c*x+d;
}


/*
 * $Log: DSTrackingAction.cc,v $
 * Revision 1.6  2015/04/28 10:15:52  pagnes
 * SiPM photoelectrons stored in DS20k (nph structure)
 *
 * Revision 1.5  2014/12/14 13:15:29  dfranco
 * cleaning of the code
 *
 * Revision 1.4  2014/11/20 13:04:58  dfranco
 * removed daughter information from PreUserStackingAction
 *
 * Revision 1.3  2014/11/13 16:47:05  dfranco
 * removed variables which were creating conflicts with the previous version of g4ds10
 *
 * Revision 1.2  2014/10/13 18:43:57  swesterd
 * fixed veto PMT positions
 *
 * Revision 1.1  2014/05/07 12:21:05  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.21  2014/03/11 09:54:38  meregaglia
 * Added generator starting from energy deposits
 *
 * Revision 1.20  2013/08/05 03:13:52  swesterd
 * some fine tuning of bscint and veto parameters
 *
 * Revision 1.19  2013/06/22 07:21:21  dfranco
 * Fixed a bug in the photoelectron absorption in DS50. Back to the previous QE method
 *
 * Revision 1.18  2013/06/19 18:35:28  swesterd
 * added DSScintCelll and made tpc PMTs' QE and reflections work like veto PMTs
 *
 * Revision 1.17  2013/06/06 23:00:19  swesterd
 * moved veto PMT numbers from copy number to name and made a separate logical volume for each veto PMT. Assigned each physical volume in the veto a unique copy number 1wxyz, where w=1 for PMT bulbs, w=2 for PMT bases, w=3 for photocathodes, and w=0 for everything else
 *
 * Revision 1.16  2013/06/05 16:28:30  dfranco
 * Improved optics in the test detector
 *
 * Revision 1.15  2013/06/04 14:11:36  dfranco
 * Added a function returning the TPC QE in DSParameters, applied in the tracking action
 *
 * Revision 1.14  2013/05/30 12:34:53  dfranco
 * Fixing the optical properties of the TPC. Not yet concluded
 *
 * Revision 1.13  2013/05/29 10:55:58  dfranco
 * Change DSLog level from debugging to development
 *
 * Revision 1.12  2013/05/27 23:59:02  swesterd
 * added a (currently commented out) Lumirror sheath to the cryostat and introduced DSOpBoundaryProcess to try to figure out why the boundaries are being screwy, with some edits so that it can handle constant and vector properties with freaking out
 *
 * Revision 1.11  2013/05/25 07:58:23  swesterd
 * Got the veto PMT optical boundaries all working along with photocathode optical properties, added PMT quantum efficiency to DSTrackingAction, and added a function to DSTrackingAction that locates and quadratically interpolates points in data, for getting useful QEs
 *
 * Revision 1.10  2013/05/20 13:58:41  swesterd
 * updated the veto PMT geometry
 *
 * Revision 1.9  2013/05/20 09:02:58  dfranco
 * PDG encoding for OpticalPhotons re-set to 50
 * *
 *
 */

