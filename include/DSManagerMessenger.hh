#ifndef DSManagerMessenger_h
#define DSManagerMessenger_h 1

#include "G4UImessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "globals.hh"
#include "G4UIcmdWithAString.hh"

using namespace std;

class DSManager;
class G4UIcmdWithAString;

class DSManagerMessenger: public G4UImessenger {

  public:
    DSManagerMessenger(DSManager*);
   ~DSManagerMessenger();

  void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    
    G4UIdirectory          *fDirectory ;
    G4UIcmdWithAString     *fDSLogCmd;
    G4UIcmdWithAString     *fDSOpticsCmd;
    G4UIcmdWithABool       *fOverLapCmd;
    G4UIcmdWithABool       *fGDMLCmd;
    G4UIcmdWithABool       *fWritePhotonsCmd;
    G4UIcmdWithAnInteger   *fEventCounterCmd;
    G4UIcmdWithAnInteger   *fVerbosityCmd ;
    G4UIcmdWithABool       *fWriteDaughtersCmd;
    G4UIcmdWithABool       *fWriteDepositsCmd;
    G4UIcmdWithABool       *fWriteThermalElectronsCmd;
    G4UIcmdWithAnInteger   *fDaughterDepthCmd ;
    G4UIcmdWithADouble     *fTMBfractionCmd;
    G4UIcmdWithADoubleAndUnit     *fTimeCutCmd;
    G4UIcmdWithAnInteger   *fFastSimulationCmd ;
};

#endif
/*
 * $Log: DSManagerMessenger.hh,v $
 * Revision 1.3  2014/07/23 14:52:38  pagnes
 * write thermal e- and kill S1 commands added
 *
 * Revision 1.2  2014/07/16 08:23:13  pagnes
 * QE scaling to 1.0 added (/ds/manager/fast_simulation xxx)
 *
 * Revision 1.1  2014/05/07 12:20:53  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.5  2014/04/11 12:33:35  pagnes
 * command to set TMB/PC ratio inside veto added
 *
 * Revision 1.4  2014/03/19 16:37:36  dfranco
 * update external configuration reader for optics tuning
 *
 * Revision 1.3  2013/03/22 16:24:40  dfranco
 * added a command to set the "daughter depth" level to store in the output file
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
