#include "DSSteppingAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "DSLogger.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "DSStorage.hh"
#include "DSEventHandler.hh"

using namespace std;


DSSteppingAction::DSSteppingAction() {;}


void DSSteppingAction::UserSteppingAction(const G4Step* theStep) { 


  
//  if(DSLogger::GetSeverity() == DSLogger::development
//     && theStep->GetPostStepPoint()) {
//    DSLog(development) 
//if(theStep->GetPostStepPoint())
/* cout     << " " <<theStep->GetTrack()->GetDefinition()->GetParticleName() << " "
      << " " <<theStep->GetTrack()->GetDefinition()->GetPDGEncoding() << " "
      <<  " E: " <<theStep->GetTrack()->GetKineticEnergy()/keV<< " keV; "
      << " Edep: " << theStep->GetTotalEnergyDeposit ()/keV << " "  
      //<<  " " <<theStep->GetTrack()->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName()<< " "
      <<  " " <<theStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial()->GetName()<< " "
      <<  " " <<theStep->GetPostStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial()->GetName()<< " "
      //<<  " " <<theStep->GetTrack()->GetVolume()->GetName()<< " "
      <<  " " << theStep->GetPostStepPoint()->GetPosition()/cm << " " 
      <<  " step " <<   theStep->GetStepLength()/um << " "   
      <<  " ID: " <<theStep->GetTrack()->GetTrackID() << " " 
      <<  " Parent ID: "    <<theStep->GetTrack()->GetParentID() << " " 
      <<  " gtime: " <<theStep->GetTrack()->GetGlobalTime() << " "
     // <<  " ltime: " <<theStep->GetTrack()->GetLocalTime()/ns  << " "
     // <<  " steps: " <<theStep->GetTrack()->GetCurrentStepNumber() << " "
      <<theStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName ()
//    << endlog ;
<<endl; 
//  }
*/
  if(theStep->GetTrack()->GetGlobalTime () > DSStorage::Get()->GetTimeCut() ) theStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
  
  if(   theStep->GetTrack()->GetDefinition()->GetPDGEncoding() == 5 
     || theStep->GetTrack()->GetDefinition()->GetPDGEncoding() == 6 
     || theStep->GetTrack()->GetDefinition()->GetPDGEncoding() == 27 
     || theStep->GetTrack()->GetDefinition()->GetPDGEncoding() == 28 
     || theStep->GetTrack()->GetDefinition()->GetPDGEncoding() == 43 
     || theStep->GetTrack()->GetDefinition()->GetPDGEncoding() == 44   
  ) theStep->GetTrack()->SetTrackStatus(fStopAndKill );
  
 // Write total deposited energy  
 
 // Write deposits 
 
 
  if(   theStep->GetTrack()->GetDefinition()->GetPDGEncoding() != 50
     && theStep->GetTotalEnergyDeposit () >0) {
     
    if(theStep->GetTrack()->GetVolume()->GetCopyNo () == 111)
      DSEventHandler::Get()->SetTPCDepEnergy(DSEventHandler::Get()->GetTPCDepEnergy() + theStep->GetTotalEnergyDeposit ()/MeV);


    if(theStep->GetTrack()->GetVolume()->GetCopyNo () == 911)
      DSEventHandler::Get()->SetVetoDepEnergy(DSEventHandler::Get()->GetVetoDepEnergy() + theStep->GetTotalEnergyDeposit ()/MeV);

    if(theStep->GetTrack()->GetVolume()->GetCopyNo () == 811)
      DSEventHandler::Get()->SetMuDepEnergy(DSEventHandler::Get()->GetMuDepEnergy() + theStep->GetTotalEnergyDeposit ()/MeV);

    if(DSParameters::Get()->GetSaturation() && theStep->GetTrack()->GetVolume()->GetCopyNo () == 911)
      DSEventHandler::Get()->SetVetoVisEnergy(DSParameters::Get()->GetSaturation()->VisibleEnergyDeposition(theStep)/MeV);

    //if(DSParameters::Get()->GetSaturation() && theStep->GetTrack()->GetVolume()->GetCopyNo () == 10001)
    //  DSEventHandler::Get()->SetVetoVisEnergy(DSParameters::Get()->GetSaturation()->VisibleEnergyDeposition(theStep)/MeV);


    if(DSParameters::Get()->GetSaturation() && theStep->GetTrack()->GetVolume()->GetCopyNo () == 811)
      DSEventHandler::Get()->SetMuVisEnergy(DSParameters::Get()->GetSaturation()->VisibleEnergyDeposition(theStep)/MeV);

    if(DSStorage::Get()->GetWriteDeposits()) {
      DSEventHandler::Get()->SetDepPID(theStep->GetTrack()->GetDefinition()->GetPDGEncoding());

      DSEventHandler::Get()->SetDepTrack(theStep->GetTrack()->GetTrackID());
      DSEventHandler::Get()->SetDepParentTrack(theStep->GetTrack()->GetParentID());


      DSEventHandler::Get()->SetDepVolume(theStep->GetTrack()->GetVolume()->GetLogicalVolume()->GetMaterial()->GetIndex());
      DSEventHandler::Get()->SetDepTotalEnergy(theStep->GetPreStepPoint()->GetKineticEnergy ()/MeV);
      DSEventHandler::Get()->SetDepKineticEnergy(theStep->GetTrack()->GetKineticEnergy ()/MeV);
      DSEventHandler::Get()->SetDepEnergy(theStep->GetTotalEnergyDeposit ()/MeV);
      DSEventHandler::Get()->SetDepStep(theStep->GetStepLength()/um);
      DSEventHandler::Get()->SetDepTime(theStep->GetTrack()->GetGlobalTime ()/ns);
      DSEventHandler::Get()->SetDepPosition(theStep->GetPostStepPoint()->GetPosition()/cm);
      DSEventHandler::Get()->SetDeposits();

    }
  }
}



/*
 * $Log: DSSteppingAction.cc,v $
 * Revision 1.5  2015/04/28 11:46:43  pagnes
 * tpcene (sum of the deposits in ActiveLAr volume) un-commented
 *
 * Revision 1.4  2014/11/21 10:19:00  dfranco
 * added a command to scale the veto scintillation yield factor and fixed the visible energy variable in the veto
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
 * Revision 1.13  2013/08/27 07:13:14  swesterd
 * add visible energy for the neutron veto
 *
 * Revision 1.12  2013/08/06 13:58:20  dfranco
 * Added 3 variables to store the deposited energy in LAr, scintillator, and water. The last two are not yet implemented. g4rooter has been updated with 3 new variables: tpcene, vetoene, and muene
 *
 * Revision 1.11  2013/07/24 09:49:02  dfranco
 * Added S1 and S2 equivalent energies, the material index for each deposit, the command killS1S2 to kill photons and electrons generated by DSLight (after storing the equivalent energies)
 *
 * Revision 1.10  2013/07/07 09:52:39  dfranco
 * removed a cout
 *
 * Revision 1.9  2013/05/31 13:43:51  dfranco
 * Change the DSLog(development) to add pre and post step volumes
 *
 * Revision 1.8  2013/05/29 10:55:58  dfranco
 * Change DSLog level from debugging to development
 *
 * Revision 1.7  2013/05/27 11:28:20  dfranco
 * Segmentation fault bug fixed, for the TPC case only. Useless files removed. Fixed the RINDEX, which was not defined in the whole range, for the scintillator.
 *
 *
 */
