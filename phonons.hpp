// phonons.hpp
//
// This file develops the Phonon class, which represents a phonon in
// flight. The types of information managed by this class include the
// phonon's position, direction-of-travel, and polarization type and
// orientation.  Additionally, this class handles transformations,
// meaning changes in direction and orientation as might be required
// by intereaction with a reflecting/transmitting surface or with a
// scatterer.  The transformation is specified via a "relative
// phonon", which is the same as a regular phonon except the direction
// and orientation componenents are considered deltas to apply to some
// other phonon.  See commentary surrounding the transform() method.
//
// This file does not participate in any namespaces. (All constructs
// are defined in the global namespace.)
//
#ifndef PHONONS_H_
#define PHONONS_H_
//
#include <cmath>
#include "geom.hpp"
#include "raytype.hpp"

//////
// CLASSES: -- Forward Declarations --
//
//   FROM OTHER HEADERS:  (Referenced here by pointer only - no
//                         need to include full header.)

class TravelRec;   /* Defined in media.hpp */
class MediumCell;  /* Defined in media.hpp */
class CellFace;    /* Defined in media.hpp */
class Scatterer;   /* Defined in scatterers.hpp */


////////////////////////////////////////////////////////////////////////
//**********************************************************************
//* %%%%                                                           %%%%
//* %%%%    CLASS Definitions                                      %%%%
//* %%%%                                                           %%%%
//*         Including:
//*
//*         o  class Phonon
//*
//********************************************************************** 
//

//////
// CLASS:   ::::  Phonon  ::::
// ENCAPS:  Elastic phonons.
//
// DEPENDENCIES:
//  The Phonon class knows about depends on:
//      o  Geometric primatives (defined in geom.hpp)
//      o  MediumCell objects   (defined in medium.hpp)   [by ref only]
//      o  Scatterer objects    (defined in scatters.hpp) [by ref only]
//
class Phonon {

  friend class DataReporter;  // For reporting to the outside world


private:

  // ::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Private Member Variables  (Phonon Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::

  unsigned long      // Serial identifier of phonon. Used in reporting
               mSID; // phonon paths. Not guaranteed unique if number
                     // of phonons exceeds limits of unsigned long
                     // (counter will roll over, per c++ standards for
                     // unsigned integer types).

  Real mTimeAlive;   // Phonon propagation sim-time (seconds). Used to
                     // report arrival time when phonon hits a
                     // collection face. Phonon will "die" if this
                     // time exceeds a time-to-live limit set in class
                     // static 'cm_ttl' member.

  Real mPathLength;  // Cumulative path length travelled by this
                     // phonon

  Real mAmplitude;   // Phonon amplitude

  R3::XYZ      mLoc; // Phonon's current Location
                     // 

  S2::ThetaPhi mDir; // Direction of Travel 
                     // 
                     // Also establishes r^hat direction in the
                     // coordinate system in which polarization is
                     // referenced: The mutually orthogonal
                     // {"outwards", "southwards", "eastwards"} define
                     // P, SV, and SH polarization directions. See
                     // comments for mPol for more.
                     // 

  Real         mPol; // Polarization angle in radians
                     // 
                     // '0' means theta^hat direction ("southwards"),
                     // 'pi/2' means phi^hat direction ("eastwards"),
                     // in the coordinate system defined w.r.t. mDir.
                     // 

  raytype     mType; // Polarization type: P, SH, SV are passed to
                     // constructor, but internally we only retain P
                     // or S.

  MediumCell * mpCell; // Points to the MediumCell in which the phonon
                       // is currently located.  Set by the InsertInto()
                       // method, which is called for propagatable
                       // phonons upon generation (typically from an
                       // EventSource object) or upon transition from
                       // one cell to another.  Relative phonons do
                       // not need to have this value set.

  Scatterer * mpScat;  // Points to the Scatterer that is active
                       // inside the current MediumCell.  Set
                       // automatically by the InsertInto() method.

private:

  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Class-Static Member Variables  (Phonon Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::

  static unsigned long cm_phonon_counter;  // Keep track inside class
                                           // of how many phonons have
                                           // been constructed.

  static Real cm_ttl;         // Phonon time-to-live. Phonon will
                              // "die" if propagation sim-time exceeds
                              // this amount.

  static Real cm_min_theta;   // These determine the colattitude bounds
  static Real cm_max_theta;   // for the direction of phonons.  Setting
                              // them to something other than 0 and Pi
                              // can create a small window-of-avoidance
                              // around the singular straight-up and
                              // straight-down directions. They are
                              // given initial values in phonons.cpp,
                              // but can be overridden with the public
                              // static method SetVerySmallTheta().

public:

  // ::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (Phonon Class) :::
  // ::::::::::::::::::::::::::::::::::::

  // Note: the distinction between SH and SV is important only on
  // construction of the phonon, but internally the polarization is
  // kept as a continuos value representing angle from vertical (see
  // below for clarification).  Thus, after construction, only this
  // value matters for S phonons. (And the value can vary smoothly as
  // propagation may rotate it, making the V-H distinction
  // non-stationary) As such, internally we record the type only as P
  // or S.

  Phonon(const R3::XYZ & loc, const S2::ThetaPhi & dir, raytype pol_t) :
    mTimeAlive(0),
    mPathLength(0),
    mAmplitude(1.0),
    mLoc (loc), 
    mDir (dir), 
    mPol (pol_t == RAY_SH ? Geometry::Pi * 0.5 : 0.0),
    mType (pol_t == RAY_P ? RAY_P : RAY_S), 
    mpCell (0)
  { 
    mSID = cm_phonon_counter++;
    nudge_if_singular();
  }

  Phonon(const S2::ThetaPhi & dir, raytype pol_t) :
    mTimeAlive(0),
    mPathLength(0),
    mAmplitude(1.0),
    mLoc (0,0,0), 
    mDir (dir), 
    mPol (pol_t == RAY_SH ? Geometry::Pi * 0.5 : 0.0),
    mType (pol_t == RAY_P ? RAY_P : RAY_S),
    mpCell (0)
  {
    mSID = cm_phonon_counter++;
    nudge_if_singular();
  }


public:

  // ::::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Set Methods  (Phonon Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::

  void SetLocation(R3::XYZ loc) {mLoc = loc;}
  void SetPolarization(Real pol) {mPol = pol;}

  static void SetVerySmallTheta(Real smalltheta) {
    cm_min_theta = smalltheta;
    cm_max_theta = Geometry::Pi - smalltheta;
  }

  static void SetTimeToLive(Real ttl) {
    cm_ttl = ttl;
  }


  // :::::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Read Methods  (Phonon Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::

  R3::XYZ      GetLocation() const  {return mLoc;}
  S2::ThetaPhi GetDirection() const {return mDir;}
  Real  GetTimeAlive() const {return mTimeAlive;}
  Real  GetAmplitude() const {return mAmplitude;}
  Real  GetPolarization() const {return mPol;}
  raytype GetRaytype() const {return mType;} // Returns RAY_P or RAY_S
  raytype GetFullRaytype() const {  
    // Extra check to distinguish SV from SH polarization. Uses the
    // somewhat-cryptic "ternary operator", which works like this:
    // (<cond> ? <result if true> : <result if false>)
    using Geometry::Pi90;
    using Geometry::Pi45;
    using std::abs;
    return (mType==RAY_P ? RAY_P 
                         : (abs(abs(mPol)-Pi90) < Pi45) ? RAY_SH 
                                                        : RAY_SV);
  }

  Real Velocity() const;
  //            Returns velocity at current location of phonon
  //

  R3::XYZ DirectionOfMotion() const;
  //            Returns the direction of particle motion, which is not
  //            necessarily the same as the direction of propagation
  //            (though it is for P-phonons of course).
  //


  // :::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Object-Manipulation Methods  (Phonon Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::::

  void Transform(const Phonon & rel); // Take 'rel' as a "relative
                                      // phonon" coding a deflection
                                      // and rotation. Rotatate
                                      // ourselves into rel.

  void PseudoReflect(const CellFace * pFace);
                      // A surface reflection method for
                      // prototyping purposes.  Not realy
                      // physical.  But does handle P<->S
                      // conversions.

  void Refract(const CellFace * pFace);
                      // Handles the transmission of a phonon from one
                      // cell into the adjacent cell, via a CellFace.
                      // This is actually a dispatcher function, which
                      // calls the appropriate specific function
                      // depending on whether ray-bending, or full
                      // reflection/transmission handling, are needed.

  void Refraction_Bend(const CellFace * pFace);
                      // Called by Refract(). Handles the case where
                      // ray bends across a velocity jump in a
                      // stair-step model.  Handles refraction only.
                      // No P/S conversions, and not full R/T
                      // treatment. Always transmits, except for
                      // post-critical incidence, which always
                      // reflects.

  void Refraction_FullRT(const CellFace * pFace);
                      // Called by Refract(). Handles the case where
                      // full Reflection/Transmission handling with
                      // raytype conversions needs to be handled.
                      // This would occur on faces joining cells that
                      // have a grid-indicated velocity discontinuity.
                      // E.g., a Moho interface, as indicated by
                      // attribute duplicity on the grid nodes.
 
  void Refraction_Continuous(const CellFace * pFace);
                      // Called by Refract(). Handles case where
                      // velocity is continuous across interface, and
                      // thus raypath needs no deflection, such as
                      // happens in cells that support velocity
                      // gradients.

  void InsertInto(MediumCell * pCell);  // Links the phonon to pCell,
                                        // updating any and all member
                                        // data needed for the phonon
                                        // to "reside" in the given
                                        // cell.

  void Propagate();   // When this method is called, a phonon will
                      // "propagate" through a "model" until such time
                      // as the phonon dies or leaves the model, at
                      // which point the method returns. Progress is
                      // reported through the DataOut object. This
                      // method is usually called from a loop that
                      // first generates a phonon from a Source event.


private:

  // :::::::::::::::::::::::::::::::::::::::
  // ::: Private Methods  (Phonon Class) :::
  // :::::::::::::::::::::::::::::::::::::::

  void Move(const TravelRec &);   // Advances a Phonon according to a
                                  // Travel Record, updating all
                                  // internal members appropriately.

  void nudge_if_singular() {
    // Adjust the direction-of-travel colatitude (theta) of the phonon
    // if it is outside the min_theta and max_theta bounds.  These
    // static members are typically set to tiny amounts > 0 and < Pi,
    // which causes the straight up and down directions (for which
    // polarization direction is ambiguous) to be avoided.
    //
    if (mDir.Theta() < cm_min_theta) mDir.SetTheta(cm_min_theta);
    if (mDir.Theta() > cm_max_theta) mDir.SetTheta(cm_max_theta);
  }

}; // class Phonon;
////

///
#endif //#ifndef PHONONS_H_
//
