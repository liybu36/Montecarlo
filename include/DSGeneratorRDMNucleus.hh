#ifndef DSGeneratorRDMNucleus_h
#define DSGeneratorRDMNucleus_h 1
#include "globals.hh"
#include "iostream"
////////////////////////////////////////////////////////////////////////////////
//
class DSGeneratorRDMNucleus
{
  // class description
  // The G4DSGeneratorRDMNucleus class is used to contain information identifying an
  // isotope (a,z,e)
  //
  // class description - end
public: // with description
  DSGeneratorRDMNucleus ();
  //    Default constructor
  //
  DSGeneratorRDMNucleus (G4int a, G4int z, G4double e);
  //    Constructor defining new isotope with A,Z.E
  //
  ~DSGeneratorRDMNucleus();
  //  Destructor
  
private:
  G4int a;
  G4int z;
  G4double e;

  //
  //
  // INLINE DECLARATIONS/DEFINITIONS:
  //
public: // with description
  inline  G4int GetA () const {return a;}
  //    Returns the value of a
  inline  G4int GetZ () const {return z;}
  //    Returns the value of z
  inline  G4double GetE () const {return e;}
  //    Returns the value of e

  //
  //
  // DECLARATIONS OF FRIENDS OF THE CLASS.
  //
  friend std::ostream &operator << (std::ostream &s, const DSGeneratorRDMNucleus &q);

};
////////////////////////////////////////////////////////////////////////////////
#endif



/*
 * Revision 1.1  2013/03/22 13:50:09  dfranco
 * Added radioactive decay and radioactive decay chain generators. They require the correspondent stacking actions. Two mac files are included as examples
 *
 */


