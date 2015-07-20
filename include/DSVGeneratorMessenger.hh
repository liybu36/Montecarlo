
#ifndef DSVGeneratorMessenger_h
#define DSVGeneratorMessenger_h 1
#include "DSIO.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
#include "G4ParticleGun.hh"

class DSPrimaryGeneratorAction;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class DSVGenerator;

using namespace std;

class DSVGeneratorMessenger: public G4UImessenger {

  public:
    DSVGeneratorMessenger(DSVGenerator*  );
   ~DSVGeneratorMessenger();

  void SetNewValue(G4UIcommand*, G4String);
    
    
  private:
    DSVGenerator*                        generator;   
    G4UIdirectory*                       fDirectory;
    G4UIcmdWith3VectorAndUnit*           fPositionCmd; 
    G4UIcmdWith3Vector*                  fDirectionCmd;
    G4UIcmdWithADoubleAndUnit*           fEnergyCmd;
    G4UIcmdWithAString*                  fParticleCmd;
    G4UIcmdWithADoubleAndUnit*           fSphereBulkCmd;
    G4UIcmdWithADoubleAndUnit*           fSphereSurfCmd;
    G4UIcmdWith3VectorAndUnit*           fSphereCentreCmd; 
    G4UIcmdWithADoubleAndUnit*           fSphereBulkRadMinCmd;
    G4UIcmdWithAString*                  fConfineCmd ;
    G4UIcmdWithAString*                  fEnergyDisTypeCmd;
    G4UIcmdWithADoubleAndUnit*           fSetEminCmd;
    G4UIcmdWithADoubleAndUnit*           fSetEmaxCmd;   
    G4UIcmdWithADouble*                  fSetAlphaCmd;  
    G4UIcmdWithADoubleAndUnit*           fSetTempCmd;   
    G4UIcmdWithADoubleAndUnit*           fSetEzeroCmd;  
    G4UIcmdWithADouble*                  fSetGradientCmd;
    G4UIcmdWithADouble*                  fSetInterCeptCmd;
    G4UIcmdWithAString*                  fDisPosTypeCmd ;
    G4UIcmdWithAString*                  fDisPosShapeCmd ;
    G4UIcmdWithADoubleAndUnit*           fSetHalfXCmd;   
    G4UIcmdWithADoubleAndUnit*           fSetHalfYCmd;   
    G4UIcmdWithADoubleAndUnit*           fSetHalfZCmd;   
    G4UIcmdWithADoubleAndUnit*           fRadiusCmd;   
    G4UIcmdWithADoubleAndUnit*           fRadius0Cmd;   
    G4UIcmdWithAString*                  fEnergyFileCmd;
    G4UIcmdWithABool*                    fUniformTPCCmd;    
    G4UIcmdWithABool*                    fUniformGasPocketCmd;    
    G4UIcmdWithABool*                    fTPCCenterCmd;    
    G4UIcmdWithAnInteger*                fGenInCryostatsCmd;
    G4UIcmdWithAnInteger*                fGenInLiquidArgonCmd;
    G4UIcmdWithAnInteger*                fGenInTeflonCmd;
    G4UIcmdWithAnInteger*                fGenInFusedSilicaCmd;
    G4UIcmdWithAnInteger*                fGenInPMTPhotocathodeCmd;
    G4UIcmdWithABool*                    fGenInG2Cryostats;
    G4UIcmdWithAnInteger*                fGenInPMTStemCmd;
    G4UIcmdWithAnInteger*                fGenInHolderSourceCmd;
    G4UIcmdWithAnInteger*                fNumberOfParticlesCmd;                               
  
  G4UIcmdWithAString*               fSetAngDistTypeCmd;
  G4UIcmdWithADoubleAndUnit*        fSetMinThetaCmd;
  G4UIcmdWithADoubleAndUnit*        fSetMaxThetaCmd;
  G4UIcmdWithADoubleAndUnit*        fSetMinPhiCmd;
  G4UIcmdWithADoubleAndUnit*        fSetMaxPhiCmd;
  G4UIcmdWithADoubleAndUnit*        fSetBeamSigmaInAngRCmd;
  G4UIcmdWithADoubleAndUnit*        fSetBeamSigmaInAngXCmd;
  G4UIcmdWithADoubleAndUnit*        fSetBeamSigmaInAngYCmd;
  G4UIcmdWith3Vector*               fUserDefAngThetaCmd;
  G4UIcmdWith3Vector*               fUserDefAngPhiCmd;
  G4UIcmdWith3VectorAndUnit*        fSetFocusPointCmd;
  G4UIcmdWith3Vector*               fSetAngRot1Cmd;
  G4UIcmdWith3Vector*               fSetAngRot2Cmd;

};

#endif
/*
 * $Log: DSVGeneratorMessenger.hh,v $
 * Revision 1.5  2015/01/14 16:58:42  dfranco
 * added source vial and correspondent commands by Laura and Bianca; manual updated
 *
 * Revision 1.4  2014/10/22 14:04:20  dfranco
 * add method to generate events in Liqud Argon
 *
 * Revision 1.3  2014/07/18 13:54:52  dfranco
 * Added a new position generator in the Gas Pocket; calibrated the numnber of photons generated per electron in S2; added a new command to generate multiple particles in the same position
 *
 * Revision 1.2  2014/05/08 11:00:42  pagnes
 * Added generator in PMTs stem
 *
 * Revision 1.1  2014/05/07 12:20:56  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.10  2014/04/22 09:56:50  perassos
 * Added the option to generate bgd from the G2 cryostats
 *
 * Revision 1.9  2014/04/11 10:20:44  perassos
 * Added generation in materials
 *
 * Revision 1.8  2014/02/13 12:03:21  dfranco
 * added new commands for spatial distributions and fixed the manual
 *
 * Revision 1.7  2014/01/29 13:13:42  perassos
 * Update of the electric field handling and of the Nuclear Recoils generator
 *
 * Revision 1.6  2014/01/07 14:11:25  perassos
 * Added the command to simulate Ar recoils
 *
 * Revision 1.5  2013/10/15 16:39:37  dfranco
 * added the possibility to define the energy spectrum reading an ascii file
 *
 * Revision 1.4  2013/05/02 16:52:36  dfranco
 * Added Radius and Radius0 methods to the position distributions
 *
 * Revision 1.3  2013/05/02 16:36:14  dfranco
 * Added several position distributions
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
