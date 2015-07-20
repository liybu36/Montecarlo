#ifndef _BXMANAGER_HH
#define _BXMANAGER_HH

#include "G4RunManager.hh"
//---------------------------------------------------------------------------//

class  DSManagerMessenger;


//---------------------------------------------------------------------------//

class DSManager : public G4RunManager
{
private:
  //default constructor
  DSManager();

public:
  static DSManager* Get();


  //copy constructor
  DSManager(const DSManager &);


  //public interface



  //protected members
protected:


  //private  members
private:
     //destructor
  virtual ~DSManager();

   static DSManager *Manager;
  // Pointers to main objects.
  //G4RunManager                    *fG4RunManager;
  DSManagerMessenger              *fDSMessenger;
};
#endif
/*
 * $Log: DSManager.hh,v $
 * Revision 1.1  2014/05/07 12:20:53  dfranco
 * Migration to Geant4.10.p01
 *
 * Revision 1.2  2013/03/22 14:09:40  dfranco
 * Added the cvs logger code at the end of each file
 *
 */
