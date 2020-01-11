// media_cellface.hpp
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

  // ::::::
  // :: Geometry:    (Determines the positioning and
  // :                orientation of the face)
  //

  R3::XYZ mNormal;        // Surface normal pointing outward from
                          // MediumCell (unit magnitude)

  R3::XYZ mPoint;         // A location known to be on the cellface.


public:

  // ::::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (CellFace Class) :::
  // ::::::::::::::::::::::::::::::::::::::

  CellFace (const R3::XYZ & N1, const R3::XYZ & N2, const R3::XYZ & N3,
            MediumCell * powner);

  CellFace (const R3::XYZ & N1, const R3::XYZ & N2, const R3::XYZ & N3,
            const R3::XYZ & N4, MediumCell * powner);
                    // First three nodes define corners of a face.
                    // The fourth node determines which point is
                    // outside.


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

  R3::XYZ Normal() const {return mNormal;}

  bool IsCollectionFace() const {return mCollect;}
  bool IsReflectionFace() const {return mReflect;}
  bool GridDiscontinuity() const {return mGridDiscon;}
  bool HasNeighbor() const {return mAdjoin;}


  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Component Reference Methods  (CellFace Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::

  MediumCell & Cell() {
    return *mpCell;
  }

  MediumCell & OtherCell() {
    return *(mpOther->mpCell);
  }


  // :::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Interrogative Methods  (CellFace Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::

  Real VelocityJump(const R3::XYZ & loc) const;
                    // Returns a number characterizing the fractional
                    // velocity jump across the interface.

  Real DistToPoint(const R3::XYZ & loc) const;
                    // Computes direct (shortest) distance from the
                    // face to the given point.

  Real DistToExitFace(const R3::XYZ & loc, const R3::XYZ & dir) const;
                    // Computes the directed (along a ray) distance
                    // from the point to the face.  Gives finite value
                    // if ray crosses from inside-to-out (exiting).
                    // Returns special values (+/-infinity) if the ray
                    // either crosses from out-to-in (entering) or
                    // else fails to intersect (runs parallel).

  GCAD_RetVal GetCircArcDistToFace(const Real & R, const R3::XYZ & C, const R3::Matrix & S) const;
                    // Similar to DistToExitFace except it computes
                    // distance along a circular arc path.

  // ::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Do-Something Methods  (CellFace Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::

  RTCoef GetRTBasis(const R3::XYZ & loc, const R3::XYZ & dir) const;
                    // Returns an RTCoef object with the basis
                    // information (which depends on the incidence
                    // angle to the CellFace) set

};




///
#endif //#ifndef MEDIA_H_
//
