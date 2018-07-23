// media.hpp
//
// This file develops the concepts and constructs of "elastic media"
// for seismic energy propagation.
//
// The general modeling strategy is to build an earth model out of a
// collection of medium "cells" that interconnect with each other in
// such a way as to fill the model space.  The cells can come in a
// variety of types.  Different types give different capabilities in
// terms of the complexity of the model that can be built out of those
// types.
//
// Two major aspects of modeling are determined by these model cells:
//
//  o Geometry:
//
//    Cells can be developed in a variety of shapes.  Our end goal is
//    to develop tetrahedral cells which can be combined into very
//    sophisticated 3D earth models.  However, initially we will
//    develop much simpler cells, such as spheres that represent an
//    entire universe (ie, no interconnectivity to other cells.)
//
//  o Velocity:
//
//    The manner in which seismic velocity is defined and how it
//    varies within a cell will also to be determined by the classes
//    developed here.  For example, one variation of the tetrahedral
//    cell will specify linear velocity gradients, which are
//    unambiguously determined by specifying velocities on the four
//    corners of the tetrahedra.  Simpler cell types will specify a
//    single seismic velocity that is spacially constant throughout
//    the cell.
//
// Physically speaking, scattering parameters could also be specified
// within a cell, but we take the strategy here of only maintaining a
// link to separate Scatterer and ScatterParam objects, which are to
// be maintained in an independent collection.  We do this because, in
// complex models, a small number of scattering parameter sets are
// likely to be shared over a large number of medium cells, and thus
// maintaining only links avoids storing redundant duplications of the
// scattering parameters and precomputer data (case in point: the
// fully-precomputed scattering functions can take up many kilobytes
// of memory).
//
#ifndef MEDIA_H_
#define MEDIA_H_
//
#include <stdlib.h>
#include "raytype.hpp"
#include "geom.hpp"
#include "elastic.hpp"
#include "array.hpp"

//////
// CLASSES: -- Forward Declarations --
//

class CellFace;   /* Needed by TravelRec class, which for stylistic
                     reasons I prefer to define before the CellFace
                     class */
class MediumCell; /* Needed by CellFace class */

class GCAD_RetVal;  /* Needed by CellFace function GetCircArcDist */

//   FROM OTHER HEADERS:  (Referenced here by pointer only -
//                         no need to include full header.)

class Scatterer;  /* Defined in scatterers.hpp */
class RTCoef;     /* Defined in rtcoef.hpp */
class GridData;   /* Passed by constant reference to Tetra constructor */


//////
// TYPEDEFS:
//

//////
// CLASSES: Definitions
//
// INCLUDING:
//
//   o  class TravelRec
//   o  class CellFace
//   o  class MediumCell
//   o  class RCUCylinder
//   o  class Tetra
//   o  class GCAD_RetVal
//

//////
// CLASS:   ::::  TravelRec  ::::
// ENCAPS:  Travel Record for a Phonon
///@ingroup phonon_support
///@brief
///
///   The TravelRec class describes a single "hop" or "leg" of the
///   journey of a Phonon along a complex trajectory.  For example,
///   from a Phonon's current location to a scatterring event, or to a
///   cell boundary, etc.
///
///   Data members are public.  (This class intended to be used more
///   like a struct than a class.)
///
class TravelRec {
public:
  
  Real     PathLength;  ///<  Distance travelled along (possibly curved)
                        ////  path
  Real     TravelTime;  ///<  Time it took to travel that distance
  R3::XYZ      NewLoc;  ///<  Location of phonon after travelling
  S2::ThetaPhi NewDir;  ///<  Direction of phonon after travelling
  Real    Attenuation;  ///<  Intrinsic Attenutation determined by Q
  CellFace    * pFace;  ///<  Points to CellFace through which the ray
                        ///   exited the MediumCell. (Undefined if N/A.)
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


//////
// CLASS:  MediumCell
//
// PURPOSE:
//
//   To serve as an abstract base class for the category of classes
//   that I will refer to as "MediumCell" classes. Here we declare the
//   virtual methods that will comprise the interface to the whole
//   family of MediumCell classes.
//
//   This arrangement is motivated by a desired ability to divide our
//   "Earth media" into "Cells" in a variety of ways, in order to
//   model both differing geometries and different approaches to
//   modeling seismic propagation.
//
//   Examples of different geometries might include cells that divide
//   up space into half-spaces, layers, or tetrahedra, etc.
//
//   Examples of different seismic modelling might be cells that treat
//   velocities as uniform-throughout, or as linear velocity
//   gradients, etc.
//
//   No actual objects of class MediumCell should ever be
//   instantiated.  (In fact they can't be, because the virtual
//   methods are not implemented.)  Rather, objects of the derived
//   classes should be instantiated.  But those objects can be pointed
//   to by MediumCell pointers.
//
//   There is a slight performace cost in using virtual methods, as
//   the call is routed through an extra level of indirection.  In
//   particular, gcc (without optimizations) seems add three
//   additional 'mov' instructions, each involving a memory access
//   (One to pull the *this pointer, one to pull the vtable pointer,
//   and one to pull the function pointer), to the calling code.  On
//   the whole I'm thinking (hoping) that this performance cost will
//   be minimal, especially with optimizations turned on.
//
//   There is a very definite memory cost, however, as each and every
//   MediumCell-derived object must contain a pointer to the vtable as
//   an additional element. This will add 8-bytes to each MediumCell,
//   or 8-MBytes for every one-million Cells in our model.  This may
//   motivate a return to non-virtual method calling (and thus
//   constraining our code to handling only ONE type of cell) after
//   thorough testing of our most-capable cell type is complete.
//
class MediumCell {
protected:

  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Class-Static Member Variables  (MediumCell Base Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  //
  //          These values need runtime initialization prior to the
  //          construction of any MediumCell-derived objects.  The
  //          responsibility to do this initialization lies with
  //          the Model constructor in model.cpp.  Initialization
  //          is achieved via public property-set methods.
  //
  
  static Real  cmPhononFreq;  // Phonon Frequency (Hertz)

public:

  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Set Methods  (MediumCell Base Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::

  static void SetFrequencyHertz(Real freq) {
    cmPhononFreq = freq;
  }

  static Real GetFrequencyHertz() {
    return cmPhononFreq; 
  }

  // :::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Interface Declarations  (MediumCell Base Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::::
  //
  //                The following pure-virtual function
  //                declarations define a set of methods (an
  //                "interface") that MUST be implemented by
  //                any derived class in order to compile.
  //

  virtual CellFace & Face(CellFace::face_id_e face_id) = 0;
                // Returns a read-write reference to the CellFace
                // object, presumed to be on the surface of the
                // MediumCell, identified by face_id

  virtual Scatterer * GetActiveScatterer() = 0;

  virtual void Link(Scatterer *) = 0;

  virtual TravelRec AdvanceLength(raytype rt, Real len, 
                                  const R3::XYZ & startloc,
                                  const S2::ThetaPhi & startdir) = 0;

  virtual TravelRec GetPathToBoundary(raytype rt,
                                      const R3::XYZ & startloc,
                                      const S2::ThetaPhi & startdir) = 0;

  virtual Real GetVelocAtPoint(const R3::XYZ & loc, 
                               raytype type) const = 0;

  virtual Real GetWavelengthAtPoint(const R3::XYZ & loc,
                                    raytype type) const =0;

  virtual Real GetDensityAtPoint(const R3::XYZ & loc) const =0;

  virtual Real GetQatPoint(const R3::XYZ & loc, raytype type) const =0;

  virtual Real IsPointInside(const R3::XYZ & loc) const =0;


};


//////
// CLASS:  RCUCylinder  ("Right-Cirular Uniform Cylinder")
//
// FROM:   MediumCell
//
// ENCAPS: 
//
//   The RCUCylinder class encapsulates a MediumCell in the shape of a
//   right-circular cylinder that is coaxial with the Z-axis, and with
//   end-caps located at specified Z-depths.  The end-caps can take
//   any orientation in space, and thus can be used to model inclined
//   interfaces between cells, such as a dipping Moho, e.g.  The top
//   and bottom endcaps will be treated as CellFaces that can be
//   marked for collection, reflection, transmission or loss, etc.
//   The curved side will always be treated as a loss face, and will
//   be mapped to a CellFace object whose location and orientation
//   have no particular meaning (since CellFace objects do not handle
//   curvature).
//
//   It will be possible to stack RCUCylinder objects on top of each
//   other to establish different velocity regions.
//
//   Seismic velocity inside the RCUCylinder is a single, spacially
//   uniform quantity, (specified separately for P and S waves).
//
class RCUCylinder : public MediumCell {
protected:

  Real  mVelTop[2];         // Velocities at "top" of cylinder
  Real  mVelBot[2];         // Velocities at "bottom" (ignored)
  Real  mDensity;           // Density (Arbitrary units - so far this
  Real  mQ[2];              // is only used in R/T computations, which
                            // depend on ratio across boundary, so
                            // units don't matter.)

  CellFace       mTopFace;  // The top endcap
  CellFace    mBottomFace;  // The bottom endcap

  Scatterer      * mpScat;  // Scatterer object that operates within
                            // this cell.  Set by Link method.
                            // Responsibility for constructing/
                            // destructing lies elsewhere.

  static Real     cmRange;  // Radius of the Cylinders (class variable)
  static bool  cmRangeSet;  // Range needs to be set prior to path
                            // calculation. As a check, constructor will
                            // throw exception if not set.  Typically,
                            // this will be set in the Model() constructor
                            // via a call to SetCylinderRange()

  static CellFace  cmLossFace;  // The one is where we discard the phonon,
                                // and "represents" a phonon leaving via
                                // the cylinder wall.  All objects will
                                // share this face.

public:

  // ::::::::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (RCUCylinder Class)  :::
  // ::::::::::::::::::::::::::::::::::::::::::

  RCUCylinder(R3::XYZ N1_top, R3::XYZ N2_top, R3::XYZ N3_top,
              R3::XYZ N1_bot, R3::XYZ N2_bot, R3::XYZ N3_bot,
              Elastic::Velocity vpvs_top,
              Elastic::Velocity vpvs_bot,
              Real rho_top,  Real rho_bot,
              Elastic::Q q_top, Elastic::Q q_bot);


  // ::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Set Methods  (RCUCylinder Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::

  virtual void Link(Scatterer * pScat) {mpScat = pScat;}

  static void SetRange(Real range) {
    cmRange = range;
    cmRangeSet = true;
  }


  // :::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Get Methods  (RCUCylinder Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::

  virtual Scatterer * GetActiveScatterer() {return mpScat;}

  virtual Real GetVelocAtPoint(const R3::XYZ &, raytype) const;
  virtual Real GetWavelengthAtPoint(const R3::XYZ &, raytype) const;
  virtual Real GetDensityAtPoint(const R3::XYZ &) const;
  virtual Real GetQatPoint(const R3::XYZ &, raytype) const;

  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Component Reference Methods  (RCUCylinder Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  virtual CellFace & Face(CellFace::face_id_e face_id);
                // Returns a read-write reference to either the
                // top or bottom CellFace object, as determined
                // by the value of face_id


  // ::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Interrogative Methods  (RCUCylinder Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::

  virtual Real IsPointInside(const R3::XYZ & loc) const;
                // Reports whether a given point is "inside"
                // the cell by returning a "mismatch" factor.


  // :::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Do-Something Methods  (RCUCylinder Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::

  virtual TravelRec AdvanceLength(raytype rt, Real length, 
                                  const R3::XYZ & startloc,
                                  const S2::ThetaPhi & startdir);

  virtual TravelRec GetPathToBoundary(raytype rt,
                                      const R3::XYZ & startloc,
                                      const S2::ThetaPhi & startdir);


protected:

  // :::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Class-Specific Helpers  (RCUCylinder Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::::

  static Real HelperAttenuation(Real cycles, Real Q);
                    // Compute amplitude attenuation factor.

};


//////
// CLASS:  Tetra  
//
// FROM:   MediumCell
//
// ENCAPS:
//
//
class Tetra : public MediumCell {
protected:
 
  R3::XYZ   mVelGrad[RAY_NBT];  // Velocity Gradient for P and S
  R3::XYZ   mDensGrad;          // Density Gradient
  Real      mVel0[RAY_NBT];     // Velocity at (0,0,0)
  Real      mDens0;             // Density at (0,0,0)
  Real      mQ[2];              // Average Q of four nodes


  Array::Quad<CellFace> mFaces; // Array consisting of 
                                // four Cell Faces
  


  Scatterer           * mpScat;  // Scatterer object that operates within
                                 // this cell.  Set by Link method.
                                 // Responsibility for constructing/
                                 // destructing lies elsewhere.

  

public:

  // ::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (Tetra Class)  :::
  // ::::::::::::::::::::::::::::::::::::
  
  Tetra (R3::XYZ N1, R3::XYZ N2, R3::XYZ N3, R3::XYZ N4,
         const GridData & dataA, const GridData & dataB,
         const GridData & dataC, const GridData & dataD);


  // :::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Set Methods  (Tetra Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::
  
  virtual void Link(Scatterer * pScat) {mpScat = pScat;}


  // :::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Get Methods  (Tetra Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::
  
  virtual Scatterer * GetActiveScatterer() {return mpScat;}

  virtual Real GetVelocAtPoint(const R3::XYZ &, raytype) const;
  virtual Real GetWavelengthAtPoint(const R3::XYZ &, raytype) const;
  virtual Real GetDensityAtPoint(const R3::XYZ &) const;
  virtual Real GetQatPoint(const R3::XYZ &, raytype) const;


  // ::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Component Reference Methods  (Tetra Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::

  virtual CellFace & Face(CellFace::face_id_e face_id);
                // Returns a read-write reference to either the
                // top or bottom CellFace object, as determined
                // by the value of face_id

  // ::::::::::::::::::::::::::::::::::::::::::::
  // ::: Interrogative Methods  (Tetra Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::

  virtual Real IsPointInside(const R3::XYZ & loc) const;
                // Reports whether a given point is "inside"
                // the cell by returning a "mismatch" factor.


  // :::::::::::::::::::::::::::::::::::::::::::
  // ::: Do-Something Methods  (Tetra Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::

  virtual TravelRec AdvanceLength(raytype rt, Real len, 
                                  const R3::XYZ & startloc,
                                  const S2::ThetaPhi & startdir);
  virtual TravelRec GetPathToBoundary(raytype rt,
                                      const R3::XYZ & startloc,
                                      const S2::ThetaPhi & startdir);
};

class CoordinateTransformation {

private:

  R3::XYZ mLocRotTrans;
  R3::XYZ mTranslation;
  R3::Matrix mRotMatrix;
  Real mRadius;

public:

  CoordinateTransformation(const Real & Vxo, const R3::XYZ & g, const R3::XYZ & loc, const R3::XYZ & t){

    R3::XYZ v2 = g.Cross(t);
    R3::XYZ v1 = v2.Cross(g);
    R3::XYZ v3 = g;
    v1.Normalize();
    v2.Normalize();
    v3.Normalize();

    Real txprime = t.Dot(v1);
    Real tzprime = t.Dot(v3);

    Real s = txprime/(Vxo);
    Real R = 1/(s * g.Mag()); // radius

    R3::Matrix S = R3::RowsMatrix(v1,v2,v3);
    R3::XYZ x0rot = S * loc;
    R3::XYZ translate = R3::XYZ(x0rot.x()+R*tzprime,
                                x0rot.y(),x0rot.z()+(-1)*R*txprime);
   
    R3::XYZ x0prime = x0rot + translate.ScaledBy(-1);

    mLocRotTrans = x0prime;
    mTranslation = translate;
    mRotMatrix = S;
    mRadius = R;

  }


public:

  R3::XYZ PrimeLoc() const { return mLocRotTrans;}
  R3::XYZ Trans() const { return mTranslation;}
  R3::Matrix RotMat() const { return mRotMatrix;}
  Real Radius() const { return mRadius;}


};


//////
// CLASS:  GCAD_RetVal
//
//  This class encapsulates the return value of the GetCircArcDistToFace() member function.

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


///
#endif //#ifndef MEDIA_H_
//
