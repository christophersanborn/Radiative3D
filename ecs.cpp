// ecs.cpp
//
#include <cmath>
#include <iostream>     /* for debugging output */
#include "ecs.hpp"

//////
// EXCEPTIONS:  EarthCoords Class
//
EarthCoords::
Unimplemented::Unimplemented(std::string method_name) :
  std::runtime_error(
      std::string("EarthCoords ") + method_name +
      std::string(": Unimplemented for this mapping")
    )
{}////
////
EarthCoords::
UnknownMapCode::UnknownMapCode(std::string method_name) :
  std::invalid_argument(std::string("EarthCoords ") + method_name +
                        std::string(": Unknown Map Code"))
{}////
////
EarthCoords::
ECSMappingError::ECSMappingError(std::string method_name) :
  std::logic_error(
      std::string("Invalid Mapping Operation in current map mode: ") + method_name
    )
{}////
////


//////
// METHOD:  EarthCoords :: ExtractElevation()
//
//   Extracts elevation, defined as distance above surface reference
//   level, from an ECS coord tuple.  The method to extract this
//   coordinate depends on the encoding scheme of the tuple (I.e., the
//   coordinate system chosen).
//
Real EarthCoords::ExtractElevation(Generic ecs_loc) const {

  switch (mMapCode) {

  case ENU_ORTHO:         // All supported mapping thus far pack elevation
  case RAE_ORTHO:         // as third position in tuple.
  case RAE_CURVED:        //
  case RAE_SPHERICAL:
    return ecs_loc.x3();
    break;

  default:
    throw UnknownMapCode("ExtractElevation");
    break;

  }
}


//////
// METHOD:  EarthCoords :: RefreshCache()
//
//  Precompute all the ICS fixed reference locations and directions.
//
void EarthCoords::RefreshCache() {

  mCacheValid = false;

  switch (mMapCode) {

  case ENU_ORTHO:
  case RAE_ORTHO:
    // Do nothing; leave mCacheValid false
    break;

  case RAE_CURVED:
    mCacheEarthCenter = R3::XYZ(0,0,-mRadE);
    mCacheNorthPole = R3::XYZ(0,mRadE,-mRadE);
    mCacheNullIsland = R3::XYZ(0,0,0);
    mCacheEasternPode = R3::XYZ(mRadE,0,-mRadE);
    mCacheSingularEast = R3::XYZ(1,0,0);
    mCacheSingularUp = R3::XYZ(0,1,0);
    mCacheValid = true;
    break;

  case RAE_SPHERICAL:
    mCacheEarthCenter = R3::XYZ(0,0,0);
    mCacheNorthPole = R3::XYZ(0,mRadE,0);
    mCacheNullIsland = R3::XYZ(0,0,mRadE);
    mCacheEasternPode = R3::XYZ(mRadE,0,0);
    mCacheSingularEast = R3::XYZ(1,0,0);
    mCacheSingularUp = R3::XYZ(0,1,0);
    mCacheValid = true;
    break;

  default:
    throw UnknownMapCode("GetEarthCenter");
    break;

  }
}

//////
// METHOD:  EarthCoords :: GetUp()
//
//   Given an ICS coordinate in the current ECS mapping scheme, return
//   an ICS unit vector pointing in the upward (skyward) direction. If
//   result is singular then a fallback value is returned.
//
R3::XYZ EarthCoords::GetUp(const R3::XYZ loc) const {

  switch (mMapCode) {

  case ENU_ORTHO:
  case RAE_ORTHO:
    return R3::XYZ(0,0,+1);
    break;

  case RAE_CURVED:
  case RAE_SPHERICAL:
    return GetEarthCenter().VectorTo(loc).UnitElse(GetSingularUp());
    break;

  default:
    throw UnknownMapCode("GetUp");
    break;

  }

}


//////
// METHOD:  EarthCoords :: GetNorth()
//
//   Given an ICS coordinate in the current ECS mapping scheme, return
//   an ICS vector pointing in the northward direction.
//
R3::XYZ EarthCoords::GetNorth(const R3::XYZ loc) const {

  switch (mMapCode) {

  case ENU_ORTHO:
  case RAE_ORTHO:
    return R3::XYZ(0,+1,0);
    break;

  case RAE_CURVED:
    return GetUp(loc).Cross(GetEast(loc));
    break;

  default:
    throw UnknownMapCode("GetNorth");
    break;

  }

}


//////
// METHOD:  EarthCoords :: GetEast()
//
//   Given an ICS coordinate in the current ECS mapping scheme, return
//   an ICS vector pointing in the eastward direction.
//
R3::XYZ EarthCoords::GetEast(const R3::XYZ loc) const {

  switch (mMapCode) {

  case ENU_ORTHO:
  case RAE_ORTHO:
    return R3::XYZ(+1,0,0);
    break;

  case RAE_CURVED:
  case RAE_SPHERICAL:
  {
    const R3::XYZ ChordNorth = loc.VectorTo(GetNorthPole());
    const R3::XYZ upward = GetEarthCenter().VectorTo(loc);
    const R3::XYZ eastward = ChordNorth.Cross(upward);
    return eastward.UnitElse(GetSingularEast());
  } break;

  default:
    throw UnknownMapCode("GetEast");
    break;

  }

}


//////
// METHOD:  EarthCoords :: GetRadial()
//
//   Returns a vector parallel (in-plane) to the Earth's surface and
//   which points in the direction away from a reference location.  (Ie,
//   it is the direction of increasing "range" from the reference.)  The
//   reference location would typically be an event source (e.g. an
//   earthquake location) and the resultant radial direction would be
//   as-defined at a seismometer location, for example.
//
R3::XYZ EarthCoords::GetRadial(const R3::XYZ ref, const R3::XYZ loc) const {

  switch (mMapCode) {

  case ENU_ORTHO:   // Works the same
  case RAE_ORTHO:   // whether curved
  case RAE_CURVED:  // or not.
  case RAE_SPHERICAL:
    return GetUp(loc).Cross(GetTransverse(ref,loc));
    break;

  default:
    throw UnknownMapCode("GetRadial");
    break;

  }

}


//////
// METHOD:  EarthCoords :: GetTransverse()
//
//   Returns a vector parallel (in-plane) to the Earth's surface and
//   which points in the transverse direction with regards to a ref-
//   erence location.  (I.e., it points in the direction of stationary
//   (unchanging) "range" from the reference.)  The reference location
//   would typically be an event source (e.g. an earthquake location)
//   and the resultant transverse direction would be as-defined at a
//   seismometer location, for example.
//
//   Note that the Transverse direction can be defined to point in
//   either a clockwise or counter-clockwise direction w.r.t. the ref-
//   erence location.  We choose the clockwise direction here, even
//   though as this results in (Radial, Transverse, Z) forming a left-
//   hand basis, as this seems to be the convention in the seismology
//   community.  For our purposes, this choice has little practical
//   relevance, as Radiative3D computes envelopes, not plain seis-
//   mograms, and thus the sign of the trace is not changed by
//   negating the recording axis.
//
R3::XYZ EarthCoords::GetTransverse(const R3::XYZ ref,
                                   const R3::XYZ loc) const {

  switch (mMapCode) {

  case ENU_ORTHO:     // Math works regardless of curvature. (The only part
  case RAE_ORTHO:     // that depends on curvature is the "up" direction,
  case RAE_CURVED:    // and that's handled appropriately by GetUp().)
  case RAE_SPHERICAL: {
    R3::XYZ ChordRadial = ref.VectorTo(loc);
    R3::XYZ UnscaledTransverse = ChordRadial.Cross(GetUp(loc));
    if (UnscaledTransverse.IsSquaredZero()) { // If not well resolved,
      UnscaledTransverse = GetSouth(loc);     // then fall back on South.
    }                                         // (Happens if loc on top of ref)
    return (UnscaledTransverse.Unit());
  } break;

  default:
    throw UnknownMapCode("GetTransverse");
    break;

  }

}


//////
// METHOD:  EarthCoords :: Convert(ecs_loc)
//
//   Responsible for converting coordinates from the user-chosen Earth
//   Coordinate System (ECS) into the internal, or "model space," XYZ
//   cartesian coordinate system.  The conversion may also calls a
//   helper function to perform an Earth-flattening transformation on
//   the depth coordinate if mFlatten is true and the chosen mapping
//   supports Earth flattening.
//
R3::XYZ EarthCoords::Convert(Generic ecs_loc) const {

  Real X,Y,Z;

  switch (mMapCode) {

  case ENU_ORTHO: {               // *** Easting, Northing, Up (Elevation)
    X = ecs_loc.x1();             // ***
    Y = ecs_loc.x2();
    Z = (mFlatten) ? FlattenDepth(ecs_loc.x3())
                   : ecs_loc.x3();
    } break;


  case RAE_ORTHO: {               // *** Range, Azimuth, Elevation
    Real range = ecs_loc.x1();    // ***
    Real phi = Geometry::DtoR * (90.0 - ecs_loc.x2());
    X = range * cos(phi);
    Y = range * sin(phi);
    Z = (mFlatten) ? FlattenDepth(ecs_loc.x3())
                   : ecs_loc.x3();
    } break;


  case RAE_CURVED:
  case RAE_SPHERICAL: {           // *** Range, Azimuth, Elevation
    Real range = ecs_loc.x1();    // ***
    Real theta = range/mRadE;
    Real phi = Geometry::DtoR * (90.0 - ecs_loc.x2());
    Real r = mRadE + ecs_loc.x3();
    X = r * sin(theta) * cos(phi);
    Y = r * sin(theta) * sin(phi);
    Z = (r * cos(theta)) + GetEarthCenter().z();
    } break;


  default:                        // *** Unknown; Throw exception
                                  // ***
    throw UnknownMapCode("Convert");
    break;

  }

  return (R3::XYZ(X,Y,Z));

}


//////
// METHOD:  EarthCoords :: BackConvert(int_loc)
//
//   From internal coordinate system back to the user's choice of
//   Earth Coordinate system.  Basically, the reverse of Convert().
//
EarthCoords::Generic EarthCoords::BackConvert(const R3::XYZ int_loc) const {

  throw Unimplemented("BackConvert");

}


//////
// METHOD:  EarthCoords :: OutConvert(int_loc)
//
//   From internal coordinate system to user's choice of output
//   coordinates.
//
EarthCoords::Generic EarthCoords::OutConvert(const R3::XYZ int_loc) const {

  switch (mOutCode) { // Early return for trivial cases:

    case OUT_NOTRANSFORM:
    return Generic(int_loc.x(), int_loc.y(), int_loc.z());

    case OUT_ECS:
    return BackConvert(int_loc);

    case OUT_ENU_ORTHO:
    break; // Continue after switch

    default:
    throw UnknownMapCode("OutConvert");

  } // end switch (mOutCode)
    // We are now, basically, handling the OUT_ENU_ORTHO case.

  // *** Now we either Un-Curve, Unflatten, or return As-Is:
  // ***

  if (CurvedCoords()) {
    Real x = int_loc.x();               // x wrt Earth center and Model origin
    Real y = int_loc.y();               // y wrt Earth center and Model origin
    Real zeta = int_loc.z() - GetEarthCenter().z();    // z wrt Earth center
    Real r = sqrt(x*x + y*y + zeta*zeta);  // radius from Earth center
    Real rxy = sqrt(x*x + y*y);            // perp. radius from polar axis
    Real ranges = mRadE * acos(zeta/r);    // range at Earth surface
    Real rangesx = (rxy>0)?ranges*(x/rxy):0; // x component of surface range
    Real rangesy = (rxy>0)?ranges*(y/rxy):0; // y component of surface range
    Real elev = r - mRadE;                 // height (elevation) above surface
    return Generic(rangesx, rangesy, elev);   // Curvature removed. (Distorts
                                              // lateral distances at depth)
  } else if (mFlatten) {
    return Generic(int_loc.x(), int_loc.y(), UnflattenDepth(int_loc.z()));
  } else {
    return Generic(int_loc.x(), int_loc.y(), int_loc.z());
  }

}


//////
// METHOD:  EarthCoords :: OutConvertDirectional()
//
EarthCoords::Generic
EarthCoords::OutConvertDirectional(const R3::XYZ loc,
                                   const R3::XYZ dir) const {

  switch (mOutCode) {

  case OUT_NOTRANSFORM:
    return Generic(dir.x(), dir.y(), dir.z());
    break;

  case OUT_ECS:
    throw Unimplemented("OutConvertDirectional");
    break;

  case OUT_ENU_ORTHO:
    return Generic( dir.Dot(GetEast(loc)),  // Return projections in the local
                    dir.Dot(GetNorth(loc)), // (East, North, Up) directions
                    dir.Dot(GetUp(loc)) );  //
    break;

  default:
    throw UnknownMapCode("OutConvertDirectional");
    break;
  }

}


//////
// METHOD:  EarthCoords :: Convert(ecs_loc, properties)
//
//   Responsible for transforming elastic properties, if called for by
//   the current ECS state.  Principally, it is the Earth Flattening
//   transformation that requires this kind of modification, as plain-
//   old coordinate mappings of course have nothing to do with the
//   elastic properties. So here we apply flattening transformation if
//   mFlatten == true.
//
Elastic::HElastic
EarthCoords::Convert(const Generic ecs_loc, 
                     const Elastic::HElastic prop) const {
  if (mFlatten) {
    Elastic::Velocity   flatV = FlattenVelocity(ecs_loc, prop.getV());
    Elastic::Density flatDens = FlattenDensity(ecs_loc, prop.getDens());
    Elastic::Q          flatQ = FlattenQ(ecs_loc, prop.getQ());
    Elastic::HetSpec   flatHS = FlattenHetSpec(ecs_loc, prop.getHS());
    return Elastic::HElastic(flatV, flatDens, flatQ, flatHS);
  }
  return prop;  // else return unmodified prop:
}


//////
// METHODS:  EarthCoords :: BackConvert(int_loc, properties)
//                       :: OutConvert(int_loc, properties)
//
//   Transforms elastic properties from the ICS representation back to
//   the ECS representation or to OCS representation.  Principally,
//   this just undoes the Earth Flattening Transformation, if the EFT
//   was part of the initial Convert process.  (Elastic properties are
//   not otherwise affected by corodinate transforms, EFTs being a
//   special case.)
//
Elastic::HElastic
EarthCoords::BackConvert(const R3::XYZ int_loc, 
                         const Elastic::HElastic prop) const {
  if (mFlatten) {
    Elastic::Velocity   sphV = UnflattenVelocity(int_loc, prop.getV());
    Elastic::Density sphDens = UnflattenDensity(int_loc, prop.getDens());
    Elastic::Q          sphQ = UnflattenQ(int_loc, prop.getQ());
    Elastic::HetSpec   sphHS = UnflattenHetSpec(int_loc, prop.getHS());
    return Elastic::HElastic(sphV, sphDens, sphQ, sphHS);
  }

  return prop;  // else return unmodified prop:

}////
////
Elastic::HElastic
EarthCoords::OutConvert(const R3::XYZ int_loc, 
                        const Elastic::HElastic prop) const {

  if (mOutCode == OUT_NOTRANSFORM) {    // Then return as-is.
    return prop;                        //
  } else {                              // Else undo Earth Flattening
    return BackConvert(int_loc, prop);  // Transform (EFT), if needed.
  }                                     // (No other transforms apply
                                        // for elastic props.)
}


//////
// METHODS:  EarthCoords :: FlattenDepth()
//                       :: UnflattenDepth()
//
//   Applies an Earth-flattening transformation to a the depth value
//   provided in the argument.  Here "depth" is a negative-sense value
//   specifying distance above (positive) or below (negative) the
//   surface reference level of the Earth.
//
Real EarthCoords::FlattenDepth(Real z) const {

  Real r = mRadE + z;
  return (mRadE * log(r/mRadE));

}////
////
Real EarthCoords::UnflattenDepth(Real z_fl) const {

  Real z_sph = mRadE * expm1(z_fl/mRadE); // mRadE * (exp(z_fl/mRadE) - 1)
  return z_sph;

}


//////
// METHOD:  EarthCoords :: FlattenVelocity()
//
//   Applies the Earth-flattening transformation to the pseudo-
//   spherical (unflattened) velocities contained in v_sph.  The
//   transformation is depth-dependent, and so the un-transformed
//   depth is extracted from the untransformed coordinate tuple in
//   ecs_loc.
//
Elastic::Velocity EarthCoords::FlattenVelocity(Generic ecs_loc,
                                     Elastic::Velocity v_sph) const {

  Real z_sph = ExtractElevation(ecs_loc);
  Real r = mRadE + z_sph;   // Radius from Earth's center

  // Compute flattened velocities:
  Real FlatFactor = mRadE/r;
  Elastic::Velocity v_fl = Elastic::VpVs(v_sph.Vp()*FlatFactor,
                                         v_sph.Vs()*FlatFactor);

  return v_fl;

}////
////
Elastic::Velocity
EarthCoords::UnflattenVelocity(R3::XYZ int_loc,
                               Elastic::Velocity v_fl) const {

  Real z_flat = int_loc.z();            // Assumes an ORTHO map scheme, which
  Real z_sph  = UnflattenDepth(z_flat); //  is the only kind that makes sense
  Real r = mRadE + z_sph;               //  for Earth Flattening. TODO: write
                                        //  a GetMSDepth() method to
                                        //  future/idiot-proof.

  Real UnflatFactor = r/mRadE;
  Elastic::Velocity v_sph = Elastic::VpVs(v_fl.Vp()*UnflatFactor,
                                          v_fl.Vs()*UnflatFactor);

  return v_sph;

}


//////
// METHODS:  EarthCoords :: FlattenDensity()
//                       :: UnflattenDensity()
//
//   A do-nothing function for now, just returns its input.  There
//   seems to be a wide variety of differing transformations for
//   density, each with narrow ranges of applicability.  Need to
//   research it further, but for now going to assume that
//   not-transforming is "good enough".  Note that the density
//   transformation has no effect on travel time kinematics or
//   scattering.  It would have effect on reflection/transmission
//   coefficients, however.
//
Elastic::Density
EarthCoords::FlattenDensity(const Generic ecs_loc,
                            const Elastic::Density rho_sph) const {

  return rho_sph;  // TODO: possibly implement a real transformation
                   // factor someday.

}////
////
Elastic::Density
EarthCoords::UnflattenDensity(const R3::XYZ int_loc,
                              const Elastic::Density rho_flat) const {

  return rho_flat; // TODO: possibly implement a real transformation
                   // factor someday.

}


//////
// METHODS:  EarthCoords :: FlattenQ()
//                       :: UnflattenQ()
//
//   Do-nothing transformations for now. (Pending future need.)
//
Elastic::Q
EarthCoords::FlattenQ(const Generic ecs_loc,
                      const Elastic::Q q_sph) const {

  return q_sph;   // TODO: possibly implement a real transformation
                  // factor someday.

}////
////
Elastic::Q
EarthCoords::UnflattenQ(const R3::XYZ int_loc,
                        const Elastic::Q q_flat) const {

  return q_flat;  // TODO: possibly implement a real transformation
                  // factor someday.

}


//////
// METHODS:  EarthCoords :: FlattenHetSpec()
//                       :: UnflattenHetSpec()
//
//   Do-nothing transformations for now. (Pending future need.)
//
Elastic::HetSpec
EarthCoords::FlattenHetSpec(const Generic ecs_loc,
                            const Elastic::HetSpec hs_sph) const {

  return hs_sph;  // TODO: possibly implement a real transformation
                  // factor someday.

}////
////
Elastic::HetSpec
EarthCoords::UnflattenHetSpec(const R3::XYZ int_loc,
                              const Elastic::HetSpec hs_flat) const {

  return hs_flat; // TODO: possibly implement a real transformation
                  // factor someday.

}


//////
// METHOD:  EarthCoords :: GetXYZToLocalNEDRotation()
//
//   Produce a rotation matrix which acts on an object oriented to the
//   internal XYZ system and reorients it to align with the local NED
//   system as defined at location 'from'.
//
//   (Equivalently, this can be thought of as a basis-change matrix
//   which probes an object expressed in the NED reference and gives
//   its expression in the internal XYZ system.)
//
//   E.g., for moment tensor M = Sum{a,b} |Na> Mab <Nb| expressed in
//   the basis {|Na>} = {|North>, |East>, |Down>}, we can get the i,j
//   elements in the {|Xi>} = {|X>, |Y>, |Z>} basis by:
//
//     Mij = <Xi| (Sum{a,b} |Na> Mab <Nb|) |Xj>
//         = Sum{a,b}  <Xi|Na> Mab <Nb|Xj>
//         = Sum{a,b}    Tia   Mab   Ybj         (Y = transpose(T))
//         =                   Mij
//
//  Where the transformation elements are:
//
//     Tia = <Xi|Na>
//
R3::Matrix EarthCoords::GetXYZToLocalNEDRotation(R3::XYZ from) const {
  R3::XYZ N1 = GetNorth(from);
  R3::XYZ N2 = GetEast(from);
  R3::XYZ N3 = GetDown(from);
  R3::XYZ X1(1,0,0);
  R3::XYZ X2(0,1,0);
  R3::XYZ X3(0,0,1);
  return (
          X1.Dot(N1) * X1.Outer(X1) +   //  |X1> <X1|N1> <X1|  +
          X1.Dot(N2) * X1.Outer(X2) +   //  ...
          X1.Dot(N3) * X1.Outer(X3) +   //
          X2.Dot(N1) * X2.Outer(X1) +
          X2.Dot(N2) * X2.Outer(X2) +
          X2.Dot(N3) * X2.Outer(X3) +
          X3.Dot(N1) * X3.Outer(X1) +
          X3.Dot(N2) * X3.Outer(X2) +
          X3.Dot(N3) * X3.Outer(X3)
         );
}
