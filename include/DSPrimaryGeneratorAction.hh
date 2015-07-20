#ifndef DSPrimaryGeneratorAction_h
#define DSPrimaryGeneratorAction_h
#include "G4VUserPrimaryGeneratorAction.hh"
#include "DSVGenerator.hh"
#include <iostream>
using namespace std;
class G4Event;
class DSPrimaryGeneratorActionMessenger;



class DSPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
  public:   

      DSPrimaryGeneratorAction();
      ~DSPrimaryGeneratorAction();

  public:
      virtual void GeneratePrimaries(G4Event*);
      void SetDSGenerator(DSVGenerator *gen) { generator = gen;}



  private:
      DSPrimaryGeneratorActionMessenger* fMessenger;
      DSVGenerator*                      generator ;

};
#endif 
/*
 * $Log: DSPrimaryGeneratorAction.hh,v $
 * Revision 1.1  2014/05/07 12:20:54  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
