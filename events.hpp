// events.hpp
//
// This file develops a set of classes that represent earthquake or
// explosion seismic events as "sources" of elastodynamic phonons.
//
// The classes developed here derive from the base class developed in
// sources.hpp.  Including this file in your code module will
// automatically include sources.hpp as well.
//
//
#ifndef EVENTS_H_
#define EVENTS_H_
//
#include "sources.hpp"
#include "tensors.hpp"

//////
// CLASSES: -- Forward Declarations --
//
//   FROM OTHER HEADERS:  (Referenced here by pointer only - no
//                         need to include full header.)

class MediumCell;  /* Defined in media.hpp */


//////
// TYPEDEFS:
//

//////
// CLASSES: Definitions
//

//////
// CLASS:   ::: ShearDislocation  :::
//          Shear Dislocation Seismic Event Phonon Source
//
// FROM:    PhononSource  (Public Inheritance)
//
// ENCAPSULATES:
//
//   Moment-Tensor source of P- and S-waves from explosion or
//   shear-dislocation events.  Stores P and S probability
//   distributions, and can spit out a randomly generated take-off
//   angle if asked.
//
class ShearDislocation : public PhononSource {
protected:

  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Private Member Variables  (ShearDislocation Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  //
  //     (All the members of PhononSource PLUS:)
  //

  R3::XYZ  mLoc;      // Location of the event
  Real     mAmp_P;    // Initial amplitude multiplier for P-waves
  Real     mAmp_S;    // Initial amplitude multiplier for S-waves
   
  MediumCell * mpCell;  // MediumCell in which the Source is
                        // located. Set by the Link() method.


public:

  // ::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (ShearDislocation Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::

  ShearDislocation(Tensor::Tensor MT, R3::XYZ Loc);


  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Set Methods  (ShearDislocation Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::

  void Link(MediumCell * pCell) {
    mpCell = pCell;
  }


  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Do-Something Methods  (ShearDislocation Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::

  Phonon GenerateEventPhonon();

};


///
#endif //#ifndef EVENTS_H_
//
