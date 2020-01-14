// media_cellface.hpp
//
// The CellFace class hierarchy offers a set of classes that define the
// bounding surfaces of MediumCell objects.  CellFace is a pure virtual base
// class defining an interface for interacting with the bounding surfaces of
// various geometries.
//
//               CellFace
//        ./--------^--------\.
//     PlaneFace         SphereFace
//
//
// Additionaly a few supporting classes are defined here:
//
//     o   class GCAD_RetVal
//
//
#ifndef MEDIA_CELLFACE_H_
#define MEDIA_CELLFACE_H_
//
#include "geom.hpp"

//////
// CLASSES: -- Forward Declarations --
//
//   FROM OTHER HEADERS:  (Referenced here by pointer only -
//                         no need to include full header.)

class MediumCell; /* Defined in media.hpp */
class RTCoef;     /* Defined in rtcoef.hpp */


//////
// CLASSES: Definitions
//
// INCLUDING:
//
//   o  class GCAD_RetVal   (Helper class to pack up a complex return value)
//   o  class CellFace
//

//////
// CLASS:  GCAD_RetVal
//
//  This class encapsulates the return value of the GetCircArcDistToFace()
//  member function.
//
class GCAD_RetVal{

protected:

  Real mEntry;
  Real mExit;
  Real mHalf;
  bool mContinuous;
  //std::vector<Real> mForwardRegion;

public:

  GCAD_RetVal(Real enter, Real exit, Real half, bool continuous);

  GCAD_RetVal();

  Real Entry() const {return mEntry;}
  Real Exit() const {return mExit;}
  Real Half() const {return mHalf;}
  bool Continuous() const {return mContinuous;}

  // std::vector<Real> GetForwardRegion() const {return mForwardRegion;}

  bool Inside (const Real & theta) const;
  bool IsProper (const GCAD_RetVal & a, const GCAD_RetVal & b, const GCAD_RetVal & c) const;

};


//////
// CLASS:  CellFace
//
// ENCAPS:
//
//   This class encapsulates a single "face" of a polyhedral
//   MediumCell, and is used to store attributes of that face, such as
//   its interconnectivity to the faces of other MediumCells, spatial
//   orientation of the face, reflection/transmission properties, etc.
//
class CellFace {
public:

  // ::::::::::::::::::::::::::::::::::::::
  // ::: Enumerations  (CellFace Class) :::
  // ::::::::::::::::::::::::::::::::::::::

  enum face_id_e {            // Provides a numbering scheme to ID the
                              // faces on a given MediumCell

    F_TOP    = 0,             // Used by two-faced cells, e.g. the
    F_BOTTOM = 1,             // Cylinder cells

    FACE_A   = 0,             // Used by four-faced cells, e.g. the
    FACE_B   = 1,             // tetrahedral cells
    FACE_C   = 2,
    FACE_D   = 3
  };

protected:

  // :::::::::::::::::::::::::::::::::::::
  // ::: Member Data  (CellFace Class) :::
  // :::::::::::::::::::::::::::::::::::::

  // ::::::
  // :: Face Flags:  (These code what "happens"
  // :                at a given interface)
  //

  bool  mCollect;     // True if data-collection should occur when a
                      //     phonon interacts with this face.
  bool  mReflect;     // True if this face represents a hard-reflection
                      //     (e.g, free surface) interface.
  bool  mAdjoin;      // True if MediumCell adjoins another cell through
                      //     this face.
  bool  mGridDiscon;  // True if the grid specified an attribute
                      //   discontinuity across this interface, as
                      //   indicated by grid duplicity on the defining
                      //   nodes of the face. (Or in other words, if
                      //   the nodes have two attribute sets, one
                      //   applying "just above" the node and one
                      //   "just under".)

  // ::::::
  // :: Linkage:     (These members encode the connectivity
  // :                of faces and cells)
  //

  CellFace * mpOther;     // Points to corresponding CellFace on
                          // the adjoining cell.

  MediumCell * mpCell;    // Cell object to which this face belongs.

  // ::::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (CellFace Class) :::
  // ::::::::::::::::::::::::::::::::::::::
  //
  //   (Callable only by derived class constructors.)
  //
  CellFace(MediumCell * pOwner) :
    mCollect ( false  ),
    mReflect ( false  ),
    mAdjoin  ( false  ),
    mGridDiscon ( false ),
    mpOther  ( nullptr ),
    mpCell   ( pOwner ) {}


public:

  // ::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Set Methods  (CellFace Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::

  void SetCollect(bool state) {mCollect = state;}
  void SetReflect(bool state) {mReflect = state;}
  void SetDiscontinuous(bool state) {mGridDiscon = state;}

  void LinkTo(CellFace & other, bool disc);
                    // Handles book-keeping of linking two complimentary
                    // CellFaces together. The caller is responsible for
                    // determining whether the CellFace linkage spans a
                    // grid-indicated discontinuity.

  void LinkTo(CellFace & other);
                    // This version assumes that the discontinuity state
                    // has already been set. If either face is marked
                    // discontiuous, then this method will ensure that
                    // both will be marked as discontinuous.


  // ::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Get Methods  (CellFace Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::

  bool IsCollectionFace() const {return mCollect;}
  bool IsReflectionFace() const {return mReflect;}
  bool GridDiscontinuity() const {return mGridDiscon;}
  bool HasNeighbor() const {return mAdjoin;}


  // :::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Linkage Reference Methods  (CellFace Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::::

  MediumCell & Cell() {
    return *mpCell;
  }

  MediumCell & OtherCell() {
    return *(mpOther->mpCell);
  }

  // :::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Interrogative Methods  (CellFace Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::

  virtual R3::XYZ Normal(R3::XYZ loc) const = 0;

  virtual Real DistToExitFace(const R3::XYZ & loc, const R3::XYZ & dir) const = 0;
                    // Computes the directed (along a ray) distance
                    // from the point to the face.  Gives finite value
                    // if ray crosses from inside-to-out (exiting).
                    // Returns special values (+/-infinity) if the ray
                    // either crosses from out-to-in (entering) or
                    // else fails to intersect (runs parallel).

  Real VelocityJump(const R3::XYZ & loc) const;
                    // Returns a number characterizing the fractional
                    // velocity jump across the interface.

  RTCoef GetRTBasis(const R3::XYZ & loc, const R3::XYZ & dir) const;
                    // Returns an RTCoef object with the basis information
                    // (which depends on the incidence angle to the CellFace).
                    // NOTE: This depends on virtual method Normal(), but does
                    // not itself need to be virtual, so it's not.

};


//////
// CLASS:  PlaneFace
// FROM:   CellFace
//
//  Planar CellFace with a fixed normal direction.  Bounds RCUCylinder
//  and Tetra model cells.
//
class PlaneFace : public CellFace {
protected:

  // ::::::
  // :: Geometry:    (Determines the positioning and
  // :                orientation of the face)
  //

  R3::XYZ mNormal;        // Surface normal pointing outward from
                          // MediumCell (unit magnitude)

  R3::XYZ mPoint;         // A location known to be on the cellface.


public:

  // :::::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (PlaneFace Class) :::
  // :::::::::::::::::::::::::::::::::::::::

  PlaneFace (const R3::XYZ & N1, const R3::XYZ & N2, const R3::XYZ & N3,
             MediumCell * powner);

  PlaneFace (const R3::XYZ & N1, const R3::XYZ & N2, const R3::XYZ & N3,
             const R3::XYZ & N4, MediumCell * powner);
                    // First three nodes define corners of a face.
                    // The fourth node determines which point is
                    // outside.

  // ::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Interrogative Methods  (PlaneFace Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::

  virtual R3::XYZ Normal(R3::XYZ loc) const override {return mNormal;}

  virtual Real DistToExitFace(const R3::XYZ & loc, const R3::XYZ & dir) const override;

  Real DistToPoint(const R3::XYZ & loc) const;
                    // Computes direct (shortest) distance from the
                    // face to the given point.

  GCAD_RetVal GetCircArcDistToFace(const Real & R, const R3::XYZ & C, const R3::Matrix & S) const;
                    // Similar to DistToExitFace except it computes
                    // distance along a circular arc path.

};


//////
// CLASS:  SphereFace
// FROM:   CellFace
//
//  Spherical CellFace centered on the origin. Surface normal can either point
//  away from origin (making it an outer surface to a spherical shell) or
//  toward origin (making it an inner surface).  Bounds SphereShell family of
//  MediumCell classes.  The operations of this class depend intrinsicially on
//  the sphere being centered on the origin.
//
class SphereFace : public CellFace {
protected:

  // ::::::
  // :: Geometry:    (Determines the positioning and
  // :                orientation of the face)
  //

  Real   mRadius;   // Radius of sphere. Positive => outward pointing surface
                    // normal, negative => inward facing normal.
  Real     mRad2;   // Precomputed Radius-Squared


public:

  // :::::::::::::::::::::::::::::::::::::::
  // ::: Constructors  SphereFace Class) :::
  // :::::::::::::::::::::::::::::::::::::::

  SphereFace (Real radius, face_id_e toporbottom, MediumCell * powner);
        // Construct SphereFace as an outward-normal (toporbottom==F_TOP) or
        // inward-normal (toporbottom==F_BOTTOM) spherical surface of radius
        // 'radius'.  Argument 'radius' must be non-negative (but can be zero).

  // :::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Interrogative Methods  (SphereFace Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::

  virtual R3::XYZ Normal(R3::XYZ loc) const override;
        // Returns surface normal at point 'loc'.  Does not check if 'loc' is
        // actually on surface (generally doesn't matter).  Does not complain
        // if loc is (0,0,0) but instead returns an arbitrary unit vector in
        // that case.

  virtual Real DistToExitFace(const R3::XYZ & loc, const R3::XYZ & dir) const;

  //Real DistToPoint(const R3::XYZ & loc) const;
                    // Computes direct (shortest) distance from the
                    // face to the given point.

  //GCAD_RetVal GetCircArcDistToFace(const Real & R, const R3::XYZ & C, const R3::Matrix & S) const;
                    // Similar to DistToExitFace except it computes
                    // distance along a circular arc path.

};


///
#endif //#ifndef MEDIA_H_
//