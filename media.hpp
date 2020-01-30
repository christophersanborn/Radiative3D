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
#include "array.hpp"
#include "raytype.hpp"
#include "elastic.hpp"
#include "raypath.hpp"
#include "media_cellface.hpp"

//////
// CLASSES: -- Forward Declarations --
//
//   FROM OTHER HEADERS:  (Referenced here by pointer only -
//                         no need to include full header.)

class Scatterer;  /* Defined in scatterers.hpp */
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
//   o  class MediumCell
//   o  class RCUCylinder
//   o  class Tetra
//   o  class SphereShell
//   o  class GCAD_RetVal
//

//////
// CLASS:   ::::  TravelRec  ::::
// ENCAPS:  Travel Record for a Phonon
//
//   Encodes one "hop" or one "leg" of the journey of a phonon along a
//   ray.  For example, from a phonon's current location to a
//   scatterring event, or to a cell boundary, etc.
//
//   Data members are public.  (This class intended to be used more
//   like a struct than a class)
//
class TravelRec {
public:

  Real     PathLength;  // Distance travelled along (possibly curved)
                        //  path
  Real     TravelTime;  // Time it took to travel that distance
  R3::XYZ      NewLoc;  // Location of phonon after travelling
  S2::ThetaPhi NewDir;  // Direction of phonon after travelling
  Real    Attenuation;  // Intrinsic Attenutation determined by Q
  const CellFace * pFace;   // Points to CellFace through which the ray
                            // exited the MediumCell. (Undefined if N/A.)
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
  //    These values need runtime initialization prior to the construction
  //    of any MediumCell-derived objects.  The responsibility to do this
  //    initialization lies with the Model constructor in model.cpp.
  //    Initialization is achieved via public property-set methods.
  //

  static Real  cmPhononFreq;  // Phonon Frequency (Hertz)


  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Base Class Member Variables  (MediumCell Base Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  //
  //    Base class manages Scatterer linkage.
  //

  Scatterer      * mpScat;      // Scatterer object that operates within
                                // this cell.  Set by Link method.
                                // Responsibility for constructing/
                                // destructing lies elsewhere.
  MediumCell() :
    mpScat(nullptr){}


public:

  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Class-Wide Property-Set Methods  (MediumCell Base Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  static void SetFrequencyHertz(Real freq) {cmPhononFreq = freq;}
  static Real GetFrequencyHertz() {return cmPhononFreq;}

  // :::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Set Methods  (MediumCell Base Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::

  void Link(Scatterer * pScat) {mpScat = pScat;}

  Scatterer * GetActiveScatterer() const {return mpScat;}


  // ::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Interrogative Methods  (RCUCylinder Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::

  Real IsPointInside(const R3::XYZ & loc) const;
        // Reports whether a point is "inside" the MediumCell by returning a
        // "mismatch" factor. Positive mismatch means the point is OUTSIDE
        // the cell, and negative values mean INSIDE.  Zero implies points
        // is exactly on boundary.  Method makes use of virtualized Face()
        // method to probe all CellFaces for "insideness".

  // :::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Interface Declarations  (MediumCell Base Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::::
  //
  //    The following pure-virtual function declarations define a set of
  //    methods (an "interface") that MUST be implemented by any derived
  //    class in order to compile.
  //

  virtual CellFace & Face(Index face_id) = 0;
        // Returns a read-write reference to the virtual CellFace
        // object, presumed to be on the surface of the MediumCell,
        // identified by face_id

  virtual Count NumFaces() const = 0;
        // Gets number of faces bounding the MediumCell derivative.

  // :::
  // ::  Functions defining a Cell's elastic profile:
  // :

  virtual Real GetVelocAtPoint(const R3::XYZ & loc, raytype type) const = 0;
  virtual Real GetWavelengthAtPoint(const R3::XYZ & loc, raytype type) const =0;
  virtual Real GetDensityAtPoint(const R3::XYZ & loc) const =0;
  virtual Real GetQatPoint(const R3::XYZ & loc, raytype type) const =0;

  // :::
  // ::  Functions that compute travel paths within the Cell:
  // :

  virtual TravelRec AdvanceLength
  (raytype rt, Real len, const R3::XYZ & startloc, const S2::ThetaPhi & startdir) const = 0;
        // This one computes where/when a Phonon will end up if if follows
        // it's cell-specific ray path (which could be straight-line or
        // curved) for a given total path length 'len'.

  virtual TravelRec GetPathToBoundary
  (raytype rt, const R3::XYZ & startloc, const S2::ThetaPhi & startdir) const = 0;
        // This one computes where/when a Phonon will end up if if
        // follows it's cell-specific ray path until it hits a boundary
        // of the cell.

protected:

  // :::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Useful Helper Functions  (MediumCell Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::::
  //
  // Kitchen sink type stuff...
  //

  static Real HelperUniformAttenuation(Real cycles, Real Q);
        // Compute amplitude attenuation factor in uniform-Q cells. (Note: All
        // derived cell types are currently assuming Q is uniform.) Note:
        // Computes amplitude attenuation, NOT energy attenuation.

};


//////
// CLASS:  RCUCylinder  ("Right-Cirular Uniform Cylinder")
//
// FROM:   MediumCell
//
// ENCAPSULATES:
//
//   The RCUCylinder class encapsulates a MediumCell in the shape of a right-
//   circular cylinder that is coaxial with the Z-axis, and with end caps
//   located at specified Z-depths.  The end caps can take any orientation in
//   space, and thus can be used to model inclined interfaces between cells,
//   such as a dipping Moho.  The top and bottom endcaps will be treated as
//   CellFaces that can be marked for collection, reflection, transmission or
//   loss, etc.  The curved side will always be treated as a loss face, and
//   for efficiency a single class-wide CylinderFace object is shared between
//   all RCUCylinder instances to serve as a common loss face.
//
//   It will be possible to stack RCUCylinder objects on top of each other
//   to establish different velocity regions.
//
//   Seismic velocity inside the RCUCylinder is a single, spacially uniform
//   quantity, (specified separately for P and S waves).
//
class RCUCylinder : public MediumCell {
protected:

  Real  mVelTop[2];         // Velocities at "top" of cylinder
  Real  mVelBot[2];         // Velocities at "bottom" (ignored)
  Real  mDensity;           // Density (Arbitrary units - so far this
  Real  mQ[2];              // is only used in R/T computations, which
                            // depend on ratio across boundary, so
                            // units don't matter.)

  PlaneFace      mTopFace;  // The top endcap
  PlaneFace   mBottomFace;  // The bottom endcap

  static bool  cmRangeSet;  // Radius of loss face needs to be set prior to
                            // path calculation. As a check, constructor
                            // will throw exception if not set.  Typically,
                            // this will be set in the Model() constructor
                            // via a call to SetCylinderRange()

  static CylinderFace cmLossFace; // The one is where we discard the phonon,
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


  // :::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Set Methods  (RCUCylinder Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::

  static void SetRange(Real range) {
    cmLossFace.SetRadius(range);
    cmRangeSet = true;
  }

  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Interface-Required Public Methods  (RCUCylinder Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  virtual CellFace & Face(Index face_id) override;
  virtual Count NumFaces() const override {return 3;}

  // :::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Get Methods  (RCUCylinder Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::

  virtual Real GetVelocAtPoint(const R3::XYZ &, raytype) const override;
  virtual Real GetWavelengthAtPoint(const R3::XYZ &, raytype) const override;
  virtual Real GetDensityAtPoint(const R3::XYZ &) const override;
  virtual Real GetQatPoint(const R3::XYZ &, raytype) const override;

  // :::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Do-Something Methods  (RCUCylinder Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::

  virtual TravelRec AdvanceLength(raytype rt, Real length, const R3::XYZ & startloc,
                                  const S2::ThetaPhi & startdir) const override;

  virtual TravelRec GetPathToBoundary(raytype rt, const R3::XYZ & startloc,
                                      const S2::ThetaPhi & startdir) const override;

};


//////
// CLASS:  Tetra
//
// FROM:   MediumCell
//
// ENCAPSULATES:
//
//   A tetrahedral medium cell bounded by four planar PlaneFace CellFaces.
//   The velocites and densities internal to the cell follow linear gradients,
//   that are computed from the values of the four corners of the cell, and
//   thus imply circular curved ray paths.
//
class Tetra : public MediumCell {
protected:

  R3::XYZ   mVelGrad[RAY_NBT];  // Velocity Gradient for P and S
  R3::XYZ   mDensGrad;          // Density Gradient
  Real      mVel0[RAY_NBT];     // Velocity at (0,0,0)
  Real      mDens0;             // Density at (0,0,0)
  Real      mQ[2];              // Average Q of four nodes


  Array::Quad<PlaneFace> mFaces; // Array consisting of
                                 // four Cell Faces

public:

  // ::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (Tetra Class)  :::
  // ::::::::::::::::::::::::::::::::::::

  Tetra (R3::XYZ N1, R3::XYZ N2, R3::XYZ N3, R3::XYZ N4,
         const GridData & dataA, const GridData & dataB,
         const GridData & dataC, const GridData & dataD);

  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Interface-Required Public Methods  (Tetra Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  virtual CellFace & Face(Index face_id) override;
  virtual Count NumFaces() const override {return 4;}

  // :::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Get Methods  (Tetra Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::

  virtual Real GetVelocAtPoint(const R3::XYZ &, raytype) const override;
  virtual Real GetWavelengthAtPoint(const R3::XYZ &, raytype) const override;
  virtual Real GetDensityAtPoint(const R3::XYZ &) const override;
  virtual Real GetQatPoint(const R3::XYZ &, raytype) const override;

  // :::::::::::::::::::::::::::::::::::::::::::
  // ::: Do-Something Methods  (Tetra Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::

  virtual TravelRec AdvanceLength(raytype rt, Real len, const R3::XYZ & startloc,
                                  const S2::ThetaPhi & startdir) const override;

  virtual TravelRec GetPathToBoundary(raytype rt, const R3::XYZ & startloc,
                                      const S2::ThetaPhi & startdir) const override;
};


//////
// CLASS:  SphereShell
//
// FROM:   MediumCell
//
// ENCAPSULATES:
//
//  Spherical shell in which the velocity inside is either a degree-0 or
//  degree-2 function of radius.  (Either straight-line ray paths or circular
//  arc raypaths.)
//
//  Internal representation of velocity is v = a*r^2 + c, and likewise for
//  density.  Q values are uniform and taken from the top-surface grid data.
//
class SphereShell : public MediumCell {
protected:

  // Velocity and Elastic Structure:

  Real mVelCoefA[RAY_NBT];      // Coef on r^2
  Real mVelCoefC[RAY_NBT];      // Constant (velocity at r=0)
  Real mZeroRadius[RAY_NBT];    // Radius r at which v==0
  Real mZeroRadius2[RAY_NBT];   // RadiusSquared r^2 at which v==0
  Real mDensCoefA;
  Real mDensCoefC;
  Real mQ[RAY_NBT];

  // Geometry:

  Array::Pair<SphereFace> mFaces; // Array consisting of
                                  // two cell faces

public:

  // ::::::::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (SphereShell Class)  :::
  // ::::::::::::::::::::::::::::::::::::::::::

  SphereShell(Real RadTop, Real RadBot,
              const GridData & DataTop, const GridData & DataBot);

  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Interface-Required Public Methods  (SphereShell Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  virtual CellFace & Face(Index face_id) override;
  virtual Count NumFaces() const override {return 2;}

  // :::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Get Methods  (SphereShell Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::

  virtual Real GetVelocAtPoint(const R3::XYZ &, raytype) const override;
  virtual Real GetWavelengthAtPoint(const R3::XYZ &, raytype) const override;
  virtual Real GetDensityAtPoint(const R3::XYZ &) const override;
  virtual Real GetQatPoint(const R3::XYZ &, raytype) const override;

  // :::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Do-Something Methods  (SphereShell Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::

  virtual TravelRec GetPathToBoundary(raytype rt, const R3::XYZ & startloc, const S2::ThetaPhi & startdir) const override;
  virtual TravelRec AdvanceLength(raytype rt, Real len, const R3::XYZ & startloc, const S2::ThetaPhi & startdir) const override;

  TravelRec GetPath_Variant_D0(raytype rt, const R3::XYZ & loc, const S2::ThetaPhi & dir) const;
        // GetPath handler for uniform velocity profile.
  TravelRec GetPath_Variant_D2(raytype rt, const R3::XYZ & loc, const S2::ThetaPhi & dir) const;
        // GetPath handler for quadratic velocity profile.
  TravelRec AdvanceLength_Variant_D0(raytype rt, Real len, const R3::XYZ & startloc, const S2::ThetaPhi & startdir) const;
        // AdvanceLength handler for uniform velocity profile.
  TravelRec AdvanceLength_Variant_D2(raytype rt, Real len, const R3::XYZ & startloc, const S2::ThetaPhi & startdir) const;
        // AdvanceLength handler for quadratic velocity profile.

  // ::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Specific Helper Methods  (SphereShell Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::

  RayArcAttributes GetRayArcD2(raytype rt, const R3::XYZ & loc, const S2::ThetaPhi & dir) const;
        // Get ray arc attributes in a radial quadratic velocity profile.

};


//////
// CLASS:  CoordinateTransformation
//
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

///
#endif //#ifndef MEDIA_H_
//
