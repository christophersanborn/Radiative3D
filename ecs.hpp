// ecs.hpp
//
// This file develops an Earth Coordinate System class to provide
// translation between geographic coordinate representation and
// internal coordinate representation for input and output, as well as
// to compute reference directions (eg., "which way is up?") which may
// be position-dependent for some mapping strategies.
//
// A goal of this class is to completely abstract this conversion
// process so that the choice of curvilinear or flat coordinate
// systems, and what the exact mapping parameters will be, becomes a
// simple user choice and not something that needs to be hard-coded
// into the main body of the software.
//
// For the moment, ECS functionality is accessed via a global 'ECS'
// object (extern-declared within this file below).  This choice to
// use a global object for this purpose may change in the future.
//
#ifndef ECS_H_
#define ECS_H_
//
#include <stdexcept>
#include "geom.hpp"
#include "elastic.hpp"  /* Earth-flattening needs elastic properties */

//////
// *** FORWARD DECLARATIONS:
//

//////
// *** CLASSES:
//

//////
// CLASS:  EarthCoords
//
//   Provides coordinate-conversion and local reference frame
//   introspection services for an Earth Coordinate System which
//   establishes run-time selectable and configurable mappings between
//   a variety of "Earth-anchored" coordinate systems and the internal
//   XYZ cartesian coordintate system used during simulation.
//
//   TERMINOLOGY and ABBREVIATIONS:
//
//   A variety of terms and abbreviations are used in the public interface and
//   documentation of this class.  A brief outline follows:
//
//   o. ECS, UCS, GS: - These stand for Earth Coordinate System, User Coordi-
//                      nate System, and Grid Space, respectively.  They are
//                      synonymous and refer to the choice of coordinates used
//                      for input into the program.  I.e., this is the coord-
//                      inate system the user wishes to describe Earth models
//                      with.  The Grid family of classes (grid.hpp) reads and
//                      stores its location data in ECS coordinates.
//
//   o. ICS, MCS, MS: - These stand for Internal Coordinate System, Model Co-
//                      ordinate System, and Model Space, respectively.  This
//                      is the coordinate system used internally by the simu-
//                      lation engine.  The Model class (model.hpp) encapsu-
//                      lates the Earth model as a collection of MediumCell
//                      derived objecs (media.hpp), and they store their loc-
//                      ation data in ICS coordinates.  The ICS is a Cartesian
//                      XYZ-based rectilinear coordinate system.
//
//   o. OCS:          - This stands for Output Coordinate System, and refers
//                      to the coordinate system that the user wishes to re-
//                      cieve program output in (e.g. event reports, location
//                      metadata on Seismometers, etc.).  This may differ from
//                      the choice of Input Coordinate System, as one choice
//                      may be convenient for describing the Earth model,
//                      whereas another may be prefered for post-processing
//                      and visualization.
//
//   o. Mapping:      - This is the process of converting from coordinates in
//                      one category to another category, and in the implemen-
//                      tation here may include one or more distorting trans-
//                      formations in addition to the simple re-expression of
//                      coordinates from one system to another.  For example,
//                      if the user has requested that an Earth-flatteing
//                      transformation be applied (a la Aki and Richards EFT),
//                      then in addition to simple coordinate transformation,
//                      there is also an adjustment applied to depths and
//                      velocities that can be used as a proxy for Earth curv-
//                      ature in layered models that can't or don't include
//                      actual curved surfaces.  As another example, a modeler
//                      that does use curved surfaces might choose an output
//                      OCS that reverses the curvature to facilitate visual-
//                      ization of the output.
//
//   o. Orthogonal:   - This refers to a coordinate system in which surfaces
//                      of constant depth are "flat" in the Model Space.  Or
//                      in other words, in which Earth curvature is not
//                      simulated directly (though is may still be simulated
//                      approximately by an EFT).
//
//   o. Curved:       - This refers to a coordinate system in which surfaces
//                      of constant depth are curved surfaces in the Model
//                      Space.
//
//   o. Flattened:    - This refers to a mapping scheme in which an Orthogonal
//                      coordinate system is used (and therefore surfaces of
//                      constant depth follow flat surfaces in the model
//                      space) but in which depths are adjusted, (in tandem
//                      with velocity adjustments), by an EFT to partially
//                      approximate the effects of real curved surfaces.
//
//   o. Convert:      - This is the term we use for mapping from ECS coord-
//                      inates to ICS coordinates.  A set of overloaded member
//                      functions named Convert() are provided by the class to
//                      handle this conversion of locus data and elastic
//                      properties (for EFTs).
//
//   o. Back-convert: - This refers to the mapping from the ICS back to the
//                      ECS.  A member function called BackConvert() is
//                      provided to handle this mapping.
//
//   o. Out-convert:  - This refers to the mapping from the ICS to the OCS
//                      (Output Coordinate System). It is at this stage that
//                      "Un-Curving" or "Un-Flattening" may occur, depending
//                      on user selections, in order to facilitate visual-
//                      ization or analyisis.  A member function called
//                      OutConvert() is provided to perform this mapping.
//
//
class EarthCoords {
public:


  // :::::::::::::::::::::::::::::::::::::::::
  // ::: Member Types  (EarthCoords Class) :::
  // :::::::::::::::::::::::::::::::::::::::::

  class Generic {     // A Generic, three element coordinate tuple.  Sorta
                      // like R3::XYZ, except with absolutely zero assumptions
  private:;           // about geometric meaning and zero defined geometric
    Real mX1;         // operations. Purpose of this type is to be converted
    Real mX2;         // into an R3::XYZ-based model space by the various
    Real mX3;         // mechanisms of the EarthCoords class, and to take
                      // advantage of the C++ typing system to enforce proper
  public:;            // usage of the class.
    Generic() :
      mX1(0),mX2(0),mX3(0)
    {}

    Generic(Real x1, Real x2, Real x3) :
      mX1(x1),mX2(x2),mX3(x3)
    {}

    void SetTriple(Real x1, Real x2, Real x3) {
      mX1=x1; mX2=x2; mX3=x3;
    }

    Real x1() const { return mX1; }
    Real x2() const { return mX2; }
    Real x3() const { return mX3; }
    Real Radius(const EarthCoords& ecs) const { return ecs.ExtractRadius(*this); }
         // We can ask it its radius if we give it an ECS to interpret it by.
    bool IsNull() const {
      return ((mX1==0.0) && (mX2==0.0) && (mX3==0.0));
    }

  };


  // ::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Custom Exceptions  (EarthCoords Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::

  class Unimplemented : public std::runtime_error {
  public: // what(): "ECS [method]: Unimplemented for this mapping".
    Unimplemented(std::string method_name);
  };

  class UnknownMapCode : public std::invalid_argument {
  public: // what(): "ECS [method]: Unknown Map Code".
    UnknownMapCode(std::string method_name);
  };

  class ECSMappingError : public std::logic_error {
  public:
    ECSMappingError(std::string token);
  };

  // ::::::::::::::::::::::::::::::::::
  // ::: Enums  (EarthCoords Class) :::
  // ::::::::::::::::::::::::::::::::::

  enum earthcoords_e  // Choice of mapping between Grid Space
  {                   // and Model Space:
                      //
    ENU_ORTHO,        // No conversion; Assumes +x,+y,+z --> East,North,Up,
                      //   respectively.  Surface-level agnostic, but by
                      //   convention assumes z=0 is sea-level.
    RAE_ORTHO,        // Range, Azimuth, Elevation.  No curvature is assumed,
                      //   so range and azi get converted to x,y coords, and
                      //   elevation becomes z.  Assumes XYZ-->ENU, and
                      //   surface at z=0.
    RAE_CURVED,       // Range, Azimuth, Elevation.  Modelspace XYZ=0,0,0 is
                      //   the origin, with spherical curvature about the
                      //   point xyz=0,0,-mRadE. At origin, North is in +y
                      //   direction, East is in +x direction, and azimuth is
                      //   measured in degrees east of north.
    RAE_SPHERICAL,    // Range, Azimuth, Elevation. ECS origin is taken to be
                      //   equatorial zero-longitude (i.e. "Null Island" or
                      //   "Anker's Point" in the Gulf of Guinea).  In the ICS,
                      //   we map this point to xyz=0,0,mRadE.  We take the ICS
                      //   origin to be the Earth's radial center, and we make
                      //   +y the polar axis, so that the North Pole is as ICS
                      //   xyz=0,mRadE,0.  Azimuth is measured in degrees east
                      //   of north.
    MAP_UNSUPPORTED   // Will trigger exception if used.
                      //
  };

  enum outcoords_e {  // Choice of coordinate system for output: (Used by the
                      //   OutConvert() method to put coordinate values into a
                      //   system suitable for post-analysis or
                      //   visualization.)
                      //
    OUT_NOTRANSFORM,  // Do not transform on output. (Ie, output uses the
                      //   internal, or "model space" coordinate system.)
    OUT_ECS,          // Transform back to ECS system that user selected for
                      //   input. (Same results as would be achieved from a
                      //   call to BackConvert().)
    OUT_ENU_ORTHO     // Transform to Cartesian XYZ, undoing curvature and
                      //   Earth-Flattening transformations, if they were a
                      //   part of the initial ECS mapping.  (And if they were
                      //   not, then this produces the same results as
                      //   OUT_NOTRANSFORM.)
  };


private:
  ;
  // ::::::::::::::::::::::::::::::::::::::::
  // ::: Member Data  (EarthCoords Class) :::
  // ::::::::::::::::::::::::::::::::::::::::
  //
  //            The member variables code the TYPE of mapping and
  //            PARAMETERS of that mapping.

  earthcoords_e  mMapCode;  // Mapping Scheme to use for input.
  outcoords_e    mOutCode;  // Target coordinates for output.
  bool           mFlatten;  // Whether or not to apply Earth-flattening
                            // transformation (ignored for CURVED mappings)
  Real              mRadE;  // Radius of the Earth

      // Fixed point and fixed direction cache: must update whenever
      // mRadE or mMapCode changes.

  R3::XYZ   mCacheEarthCenter;  // Locations...
  R3::XYZ     mCacheNorthPole;  //
  R3::XYZ    mCacheNullIsland;
  R3::XYZ   mCacheEasternPode;
  R3::XYZ  mCacheSingularEast;  // Directions (unit vectors)...
  R3::XYZ    mCacheSingularUp;
  bool            mCacheValid;  // (False in mappings where these locations are meaningless)

  void RefreshCache();

public:
  ;
  // :::::::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (EarthCoords Class) :::
  // :::::::::::::::::::::::::::::::::::::::::

  EarthCoords() :
    mMapCode ( ENU_ORTHO     ),
    mOutCode ( OUT_ENU_ORTHO ),
    mFlatten ( false         ),
    mRadE    ( 6371.0        )
  { RefreshCache(); }


  // ::::::::::::::::::::::::::::::::::::::::
  // ::: Set-Methods  (EarthCoords Class) :::
  // ::::::::::::::::::::::::::::::::::::::::
  //
  //        These methods form the interface by which the user of the class
  //        selects the desired coordinate and transformation system.
  //

  void SetEarthFlattening(bool truefalse) {mFlatten = truefalse;}
  void SetEarthRadius(Real erad) {mRadE = erad; RefreshCache();}
  void SetMapping(earthcoords_e maptype) {mMapCode = maptype; RefreshCache();}
  void SetOCSMapping(outcoords_e maptype) {mOutCode = maptype;}


  // ::::::::::::::::::::::::::::::::::::::::
  // ::: Get-Methods  (EarthCoords class) :::
  // ::::::::::::::::::::::::::::::::::::::::

  Real GetEarthRadius() const {return mRadE;}
  bool IsEarthFlattening() const {return mFlatten;}
  bool CurvedCoords() const {
    // Returns true if selected mapping involves curvature.
    return (mMapCode == RAE_CURVED || mMapCode == RAE_SPHERICAL);
  }


  // ::::::::::::::::::::::::::::::::::::::::::::
  // ::: Extract Methods  (EarthCoords Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::
  //
  //        These extract named coordinates from the coordinate tuple in the
  //        currently selected ECS coordinate system.  Note: there are no
  //        coordinate-space to model-space conversions here - we stay in
  //        ECS. This is just about getting a particular coordinate from the
  //        tuple.  (Some basic translation may occur, however. E.g.,
  //        extracting "East" or "North" from an RAE tupple.)
  //

  Real ExtractElevation(Generic ecs_loc) const;
  Real ExtractRadius(Generic ecs_loc) const;


  // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Get Reference Location Methods  (EarthCoords Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  //
  //        These return named fixed locations and special fixed reference
  //        directions in the internal coordinate system. "NullIsland"
  //        refers to the equatorial zero (I.e. where LLE=0,0,0), and
  //        "EasternPode" refers to LLE=0,+90,0.  "Singular East" is the
  //        direction used for East when it is otherwise undefined (when on
  //        the polar axis).
  //

  R3::XYZ GetEarthCenter() const {if (mCacheValid) return mCacheEarthCenter; else throw ECSMappingError("EarthCenter");}
  R3::XYZ GetNorthPole() const {if (mCacheValid) return mCacheNorthPole; else throw ECSMappingError("NorthPole");}
  R3::XYZ GetNullIsland() const {if (mCacheValid) return mCacheNullIsland; else throw ECSMappingError("NullIsland");}
  R3::XYZ GetEasternPode() const {if (mCacheValid) return mCacheEasternPode; else throw ECSMappingError("EasternPode");}
  R3::XYZ GetSingularEast() const {if (mCacheValid) return mCacheSingularEast; else throw ECSMappingError("SingularEast");}
  R3::XYZ GetSingularUp() const {if (mCacheValid) return mCacheSingularUp; else throw ECSMappingError("SingularUp");}

  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Get-Local-Direction Methods  (EarthCoords Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  //
  //        These return named directions as vectors in the internal co-
  //        ordinate system.  For curvilinear mappings, these directions
  //        are location-dependant.
  //
  //        Note that the location arguments are assumed to be in the
  //        INTERNAL coordinate system.
  //
  //        The GetRadial() and GetTransverse() methods compute the Radial
  //        and Transverse directions (as in the Radial, Transverse, and Z
  //        directions commonly used as a basis in seismic analysis) defined
  //        at the location 'from' with respect to a reference location at
  //        'ref'.  (Typically 'ref' will be the event location and 'from'
  //        the seismometer location.)  A GetZ() method is not provided as
  //        the direction is definitionally equivalent to the GetUp()
  //        direction and is also not dependent on the 'ref' location.  Note
  //        that in seismology, the 'transverse' coordinate typically
  //        increases clockwise, and so we do the same here, meaning RTZ
  //        will form a LEFT-handed coordinate system.
  //
  //        SINGULARITIES: Directional singluarities are removed as follows:
  //        EAST: Along the polar axis, East is defined as paralell to the
  //        line from the Earth's center to the Eastern Pode.  NORTH: Along
  //        the polar axis, North is defined by Up cross East (same as it is
  //        defined elsewhere).  UP: At the Earth's center, Up is defined as
  //        along the line from the Earth's center to the North Pole.
  //        TRANSVERSE: When 'ref' and 'from' are the same or above/below
  //        each other according to local Up direction, then Transverse is
  //        assigned to point South.
  //

  R3::XYZ GetUp(R3::XYZ from) const;
  R3::XYZ GetNorth(R3::XYZ from) const;
  R3::XYZ GetEast(R3::XYZ from) const;
  R3::XYZ GetDown(R3::XYZ from) const  {return GetUp(from).Negative();}
  R3::XYZ GetSouth(R3::XYZ from) const {return GetNorth(from).Negative();}
  R3::XYZ GetWest(R3::XYZ from) const  {return GetEast(from).Negative();}

  R3::XYZ GetRadial(R3::XYZ ref, R3::XYZ from) const;
  R3::XYZ GetTransverse(R3::XYZ ref, R3::XYZ from) const;


  // :::::::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Coordinate Transformations  (EarthCoords Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::::::::
  //
  //        To/From the internal coordinate system.  These are generally
  //        intended to be called by GridCell objects when queried via their
  //        Get methods, in order to convert on the fly from ECS coordinates
  //        into internal representation.  Thus, the grid will exist in the
  //        user-chosen ECS coordinates, but the model will be built
  //        automatically in internal coordinates.
  //

  R3::XYZ Convert(Generic ecs_loc) const;     // From ECS to internal
  Generic BackConvert(R3::XYZ int_loc) const; // From internal to ECS
  Generic OutConvert(R3::XYZ int_loc) const;  // From internal to prefered
                                              // coordinate system for output.

  Generic OutConvertDirectional(R3::XYZ int_loc, R3::XYZ dir) const;
                // This converts a directional quantity 'dir' at the ICS
                // location 'int_loc' into OCS representation.  (Principle use
                // case is seismometer axes, which are included in the output
                // files, and may need to be rotated to undo Earth curvature.)


  // :::::::::::::::::::::::::::::::::::::::::::::
  // ::: Earth-Flattening  (EarthCoords Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::
  //
  //        These are for performing Earth-flattening algorithms a la Aki
  //        and Richards box 9.9.  Note that the depth transformation will
  //        be called automatically by the Convert() method if the chosen
  //        ECS system includes flattening, but the property transformations
  //        (velocity) are handled via separate functions. To facilitate
  //        property transformations, the Convert() method is overloaded to
  //        take and return an Elastic::HElastic object, and will flatten or
  //        not based on the current ECS flattening settings.  If flattening
  //        is not called for, then the Convert() method returns the
  //        HElastic object unmodified, as the elastic properties are not
  //        otherwise coordinate-dependent.
  //
  //        The Convert() function is intended to be called by the GetData()
  //        method of the GridNode class, so that transparent, on-the-fly
  //        conversion happens before the data is passed to the MediumCell
  //        constructors that depend on the data. (The specific flattening
  //        methods are mostly not needed outside this class itself.)
  //
  //        The property transforms must be supplied with the unconverted/
  //        untransformed (ECS frame) locations at which the properties are
  //        to be transformed, and for the reverse transformations the
  //        internal coords must be supplied.
  //
  //

  Elastic::HElastic Convert(Generic ecs_loc, Elastic::HElastic prop) const; // From ECS to ICS
  Elastic::HElastic BackConvert(R3::XYZ int_loc, Elastic::HElastic prop) const; // From ICS to ECS
  Elastic::HElastic OutConvert(R3::XYZ int_loc, Elastic::HElastic prop) const; // From ICS to OCS

  Elastic::Velocity FlattenVelocity(Generic ecs_loc, Elastic::Velocity v_sph) const;
  Elastic::Velocity UnflattenVelocity(R3::XYZ int_loc, Elastic::Velocity v_fl) const;

  Elastic::Density FlattenDensity(Generic ecs_loc, Elastic::Density rho) const;
  Elastic::Density UnflattenDensity(R3::XYZ int_l, Elastic::Density rho) const;
  Elastic::Q FlattenQ(Generic ecs_loc, Elastic::Q q_sph) const;
  Elastic::Q UnflattenQ(R3::XYZ int_l, Elastic::Q q_fl) const;
  Elastic::HetSpec FlattenHetSpec(Generic ecs_loc, Elastic::HetSpec hs) const;
  Elastic::HetSpec UnflattenHetSpec(R3::XYZ int_l, Elastic::HetSpec hs) const;
  //
  //    Density, Q, and HetSpec: Currently do-nothing placeholders.  Uncertain
  //    whether there *is* any good flattening transform for density.  As for
  //    Q and Heterogeneity Spectrum, would have to look into this further,
  //    but I think the appropriate transformation would include introducing
  //    a directional anisotropy in the Q and HetSpec parameters.  This is a
  //    capability not currently implemented, but which is being considered
  //    for future implementation.
  //

  Real FlattenDepth(Real z) const;    // Note: z is negative-sense, ie more
  Real UnflattenDepth(Real z) const;  // negative means more deeper, whereas
                                      // positive means above surface
                                      // reference level.


  // :::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Get-Rotation Methods  (EarthCoords Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::

  //            These return rotation matrices for rotating vectors or
  //            tensors to some desired orientation in the local Earth
  //            reference frame.
  //
  //            Note that the 'from' argument is assumed to be in
  //            the INTERNAL coordinate system.
  //

  R3::Matrix GetXYZToLocalNEDRotation(R3::XYZ from) const;
                // Produce a rotation matrix which acts on an object
                // oriented to the internal XYZ system and reorients
                // it to align with the local NED system as defined at
                // location 'from'.
                //
                // (Equivalently, this can be thought of as a
                // basis-change matrix which probes an object
                // expressed in the NED reference and gives its
                // expression in the internal XYZ system.)
                //


}; // class EarthCoords
////
extern EarthCoords ECS;  // Global object
///                      // (initialized in global.cpp)

///
#endif //#ifndef ECS_H_
//
