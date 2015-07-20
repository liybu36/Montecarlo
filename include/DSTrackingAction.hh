#ifndef DSTRACKINGACTION_H
#define DSTRACKINGACTION_H


#include "G4UserTrackingAction.hh"


class DSTrackingAction : public G4UserTrackingAction {

  public:

	DSTrackingAction();
	virtual ~DSTrackingAction(){};
	//virtual void PreUserTrackingAction(const G4Track*);
	virtual void PostUserTrackingAction(const G4Track*);

  private:
  
       //binary search
       int binSearch(double* , double , int , int );
       int search(double* , double , int );
       //quadratically interpolate a point in an array
       double interpolate(double* , double*, double, int );

};

#endif
/*
 * $Log: DSTrackingAction.hh,v $
 * Revision 1.3  2014/12/14 13:15:33  dfranco
 * cleaning of the code
 *
 * Revision 1.2  2014/11/20 13:05:05  dfranco
 * removed daughter information from PreUserStackingAction
 *
 * Revision 1.1  2014/05/07 12:20:55  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.5  2013/06/19 18:35:25  swesterd
 * added DSScintCelll and made tpc PMTs' QE and reflections work like veto PMTs
 *
 * Revision 1.4  2013/05/25 07:58:22  swesterd
 * Got the veto PMT optical boundaries all working along with photocathode optical properties, added PMT quantum efficiency to DSTrackingAction, and added a function to DSTrackingAction that locates and quadratically interpolates points in data, for getting useful QEs
 *
 * Revision 1.3  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
