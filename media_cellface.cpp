// media_cellface.cpp
//
// The CellFace class hierarchy offers a set of classes that define the
// bounding surfaces of MediumCell objects.  CellFace is a pure virtual base
// class defining an interface for interacting with the bounding surfaces of
// various geometries. See full commentary in media_cellface.hpp.
//
#include <cassert>
#include "media_cellface.hpp"
#include "media.hpp"
#include "rtcoef.hpp"

// In this file:
// CLASS IMPLEMENTATIONS FOR:
//
//   o  Class CellFace
//   o  Class GCAD_RetVal
//

//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  CLASS:  CellFace                                            ****
// ****                                                              ****
//
//   Methods Defined Here:
//
//     o   LinkTo()
//     o   VelocityJump()
//     o   GetRTBasis()
//

//////
// METHOD:  CellFace :: LinkTo()
//
//   Called by the model-building routines, this function handles the
//   book-keeping of linking together the two CellFaces (one from each
//   MediumCell of an adjoining pair) that represent the interface
//   between two model cells.
//
//   Takes a reference to the other CellFace, and a boolean argument
//   which signals whether the Grid specified that there should be a
//   first-order discontinuity across this interface.  The caller is
//   responsible for determining whether this is the case, and
//   communicating it via this argument.
//
void CellFace::LinkTo(CellFace & other, bool disc) {
  mpOther = &other;
  mAdjoin = true;
  mGridDiscon = disc;
  other.mpOther = this;
  other.mAdjoin = true;
  other.mGridDiscon = disc;
}//
//   If called with the following form, then assumption is that a
//   grid-indicated discontinuity across the interface has already
//   been flagged by a call to SetDiscontinuous() on one or both of
//   the faces.  (This could happen, e.g., at an earlier stage of the
//   model build-out process.)  Check mGridDiscon on both faces, and
//   if either is true, then explicitly set both to true to ensure
//   both faces agree.
//
void CellFace::LinkTo(CellFace & other) {
  bool disc = (mGridDiscon || other.mGridDiscon); // true if either is true
  mpOther = &other;
  mAdjoin = true;
  mGridDiscon = disc;
  other.mpOther = this;
  other.mAdjoin = true;
  other.mGridDiscon = disc;
}


//////
// METHOD:  CellFace :: VelocityJump()
//
//   Returns a number characterizing the size of the velocity
//   discontinuity across the CellFace (from one MediumCell to
//   another) at a given location.  Computes the number as a
//   fractional velocity change. Returns a non-negative number, and
//   returns the largest of the fractional changes (P or S) at the
//   given location.
//
Real CellFace::VelocityJump(const R3::XYZ & loc) const {

  Real vel1, vel2, dvp, dvs;

  // First for P velocities:
  vel1 = mpCell->GetVelocAtPoint(loc, RAY_P);
  vel2 = mpOther->mpCell->GetVelocAtPoint(loc, RAY_P);
  dvp = std::abs(2*(vel2-vel1)/(vel2+vel1));

  // Then for S velocities:
  vel1 = mpCell->GetVelocAtPoint(loc, RAY_S);
  vel2 = mpOther->mpCell->GetVelocAtPoint(loc, RAY_S);
  dvs = std::abs(2*(vel2-vel1)/(vel2+vel1));

  return (dvp > dvs) ? dvp : dvs;  // Return the greater of dvp, dvs

}


//////
// METHOD:  GetRTBasis()
//
//   Populates an RTCoef object with a set of basis vectors in which
//   reflected or transmitted rays can be expanded; populates the
//   velocity and density parameters for both sides of the interface.
//
//   The information generated by this function and stored in the
//   RTCoef object are subsequently used (by member functions of the
//   RTCoef class or in other locations) to compute reflection and
//   transmission probabilities, and to choose an outcome.  The info
//   generated here is basically just to get the ball rolling on that
//   front.
//
//   The information generated here depends on the incidence angle
//   w.r.t. the CellFace, and needs access to MediumCells of which this
//   CellFace object has exclusive knowledge, and that's why this method
//   appears in the CellFace class, instead of in the RTCoef class or
//   elsewhere.
//
RTCoef CellFace::GetRTBasis(const R3::XYZ & loc,
                            const R3::XYZ & dir) const {

  RTCoef rt(Normal(loc), dir);  // We will return this object

  MediumCell * pCellR = mpCell; // Reflection cell
  MediumCell * pCellT;          // Transmission cell

  rt.DensityR = pCellR->GetDensityAtPoint(loc);
  rt.VelocR[RAY_P] = pCellR->GetVelocAtPoint(loc,RAY_P);
  rt.VelocR[RAY_S] = pCellR->GetVelocAtPoint(loc,RAY_S);

  if (HasNeighbor()) {          // Get values on transmit side
    pCellT = mpOther->mpCell;
    rt.DensityT = pCellT->GetDensityAtPoint(loc);
    rt.VelocT[RAY_P] = pCellT->GetVelocAtPoint(loc,RAY_P);
    rt.VelocT[RAY_S] = pCellT->GetVelocAtPoint(loc,RAY_S);
  }
  else {                        // Else: free surface.
    rt.DensityT = 0.0;          // These values simulate a free
    rt.VelocT[RAY_P] = 1e-12;   // surface to pretty good
    rt.VelocT[RAY_S] = 1e-12;   // approximation.
    rt.NoTransmit = true;       //
  }

  return rt;

}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  CLASS:  PlaneFace  (from CellFace)                          ****
// ****                                                              ****
//
//   Methods Defined Here:
//
//     ...
//

//////
// CONSTRUCTOR:  PlaneFace()
//
PlaneFace::PlaneFace (const R3::XYZ & N1, const R3::XYZ & N2,
                      const R3::XYZ & N3,
                      MediumCell * powner) :
  CellFace(powner),
  mPoint   ( N1     )       // Need a point on the Face. N1 is on
                            // the face (albeit on the edge of it.)
{
  // The three nodes, N1, N2, and N3, define a triangle that
  // represents the cell face.  The order of the nodes determines
  // the orientation of the face, with the "outward" side being the
  // side from which the nodes would be seen to trace a
  // counter-clockwise path.
  R3::XYZ Ray1 = N1.VectorTo(N2);
  R3::XYZ Ray2 = N1.VectorTo(N3);
  mNormal = Ray1.Cross(Ray2);   // Points in outward direction;
  mNormal.Normalize();          // Set magnitude to 1.0;

}


//////
// CONSTRUCTOR:  PlaneFace()
//
PlaneFace::PlaneFace (const R3::XYZ & N1, const R3::XYZ & N2,
                      const R3::XYZ & N3, const R3::XYZ & N4,
                      MediumCell * powner) :
  CellFace(powner),
  mPoint   ( N1     )

{
  // The three nodes, N1, N2, and N3, define a triangle that
  // represents the cell face.  The fourth node N4 is the excluded
  // node which is used to determine the direction of the normal vector.
  R3::XYZ Ray1 = N1.VectorTo(N2);
  R3::XYZ Ray2 = N1.VectorTo(N3);
  mNormal = Ray1.Cross(Ray2);
  mNormal.Normalize();

  if(GetDistanceAboveFace(N4) > 0)    // Flips direction of mNormal if
    mNormal = mNormal.ScaledBy(-1);   // pointing inward

}


//////
// METHOD:  PlaneFace :: GetDistanceAboveFace()
//
//   Returns the direct (shortest) distance from the CellFace plane to
//   the given point, with the return value using the following sign
//   convention:
//
//     > 0:  Point is "above" the CellFace, with "above" being the
//           side pointed to by the mNormal vector.
//     = 0:  Point is "on" the surface.
//     < 0:  Point is "below" the surface.
//
Real PlaneFace::GetDistanceAboveFace(const R3::XYZ & loc) const {
  //
  return mNormal.Dot(mPoint.VectorTo(loc));
}


//////
// METHOD:  PlaneFace :: LinearRayDistToExit()
//
//   Returns the directed (along a ray) distance from the starting
//   location ('loc') to the CellFace, under the presumption that we
//   are crossing the face in the "exiting" direction.  Should it be
//   that the direction of our ray will cause us to cross the face in
//   "entering" direction, or else to not cross the face at all
//   (parallel motion), then special values are returned to signal
//   this condition.
//
// INPUTS:
//
//   o  loc     -   The starting location of the ray
//   o  dir     -   The direction of the ray (unit vector -
//                           responsibility lies with caller
//                           to provide unit-magnitude vector)
// RETURNS:
//
//    +inf  -   Means we will NEVER cross the face in exiting fashion,
//              and our future lies INSIDE the face (ie, we are
//              entering or have entered, or are running parallel and
//              are already inside (including on-surface))
//
//     > 0  -   We will exit the face at some future time.
//
//    == 0  -   We are exiting right now (loc is on surface and dir
//              points outward)
//
//     < 0  -   We have already exited the face, at some point in the
//              past.
//
//    -inf  -   We will NEVER exit the face, because we are already
//              outside, and running parallel to the face.
//
Real PlaneFace::LinearRayDistToExit(const R3::XYZ & loc,
                                    const R3::XYZ & dir) const {

  Real d_sh;        // dist_shortest

  // ::::::
  // :: First get the direct (shortest) distance to the plane.
  // :
  //       d_sh  =  \hat{N} .dot. (\vec{F} - \vec{P})
  //
  //       With N being the outward surface normal from the cell face,
  //       F being a location known to reside on the surface, and P
  //       being the starting location of the ray.
  //
  //       d_sh < 0 means we are outside the surface.  Otherwise we
  //       are inside or on the surface.
  //

  d_sh = mNormal.Dot(loc.VectorTo(mPoint));

  // ::::::
  // :: Next, get the direction factor:
  // :
  //       Sign tells us entering/exiting status, and value becomes a
  //       scale factor to convert direct to directed distance.
  //

  Real d_fact = mNormal.Dot(dir);   // d_fact > 0 means we are moving in
                                    // the exiting direction.  d_fact = 0
                                    // means we are moving parallel, and
                                    // d_fact < zero means we are moving
                                    // in the entering direction.

  // ::::::
  // :: First, handle special cases:
  // :

  if (d_fact < 0) {             // Then we are moving in the entering
    return (1./0.); /* +INF */  // direction. Treat face as invisible
  }                             // (ie, we will NEVER hit it).

  if (d_fact == 0) {            // Then we are running exactly parallel

    if (d_sh < 0) {             // Then we are outside, and will never
      return (-1./0.); /*-INF*/ // enter. Return -inf to signal that
    }                           // we "exited a long time ago".
    else {
      return (1./0.); /*+INF*/  // Else we are inside, and will never
    }                           // exit.

  } //
  ///  If we get here, then d_sh > 0.
  //

  // ::::::
  // :: Finally, calculate and return dist to face:
  // :
  //      dist = d_sh / d_fact
  //

  return d_sh / d_fact;

}



//////
// METHOD:  PlaneFace :: GetCircArcDistToFace()
//
//
//
GCAD_RetVal PlaneFace::GetCircArcDistToFace(const Real & R, const R3::XYZ & C,
                                            const R3::Matrix & S) const{

  GCAD_RetVal retval;
  bool continuous = true;

  // Rotate and Translate Cell face
  R3::XYZ rotNorm = S * mNormal;
  R3::XYZ rotx0 = S * mPoint;
  R3::XYZ x0prime = rotx0 + C.ScaledBy(-1);

  // Equation of cell face (plane):
  // ax + by + cz + d = 0
  // d = -ax0 - by0 - cz0
  Real d = (-1)*rotNorm.Dot(x0prime);

  // Distance from origin of new system,
  // bisecting intersection of face with circular arc
  Real D = -d/sqrt(rotNorm.x()*rotNorm.x()+rotNorm.z()*rotNorm.z());

  // Compute azimuthual angle to bisector
  R3::XYZ norm2D = R3::XYZ(rotNorm.x(),0,rotNorm.z());
  norm2D.Normalize();

  Real ColatAngletoBisector = atan2(norm2D.x(),norm2D.z());
  Real ColatAngletoExit = 0;  // azimuthal angle to exit
  Real ColatAngletoEntry = 0;

  if (D/R < 1 && D/R > -1){
    // In 2D rotated system, determine
    // angle to each intersection point
    // (plane-circulararc intersection)
    Real DblPrimeAngletoIntersection = acos(D/R);

    //////////////////////////////////
    /// Physical Viability of Exit ///
    //////////////////////////////////

    // Test if exit in negative velocity zone.
    // If in negative velocity zone or "behind" current location
    // it is physically impossible to reach exit.
    // +inf indicates future exit unacceptable


    if(ColatAngletoBisector > (-1)*Geometry::Pi90 && ColatAngletoBisector < Geometry::Pi90){
      ColatAngletoEntry = ColatAngletoBisector + DblPrimeAngletoIntersection;
      ColatAngletoExit = ColatAngletoBisector - DblPrimeAngletoIntersection;
      continuous = false;
    }
    else if (ColatAngletoBisector <= (-1)*Geometry::Pi90){
      ColatAngletoEntry = ColatAngletoBisector + DblPrimeAngletoIntersection;
      ColatAngletoExit = ColatAngletoBisector - DblPrimeAngletoIntersection + Geometry::Pi360;
    }
    else if (ColatAngletoBisector >= Geometry::Pi90){
      ColatAngletoEntry = ColatAngletoBisector + DblPrimeAngletoIntersection - Geometry::Pi360;
      ColatAngletoExit = ColatAngletoBisector - DblPrimeAngletoIntersection;
    }
    else{
      std::cerr << "GCAD Bisector nan" << std::endl;
      exit(1);
    }
  }

  if (ColatAngletoBisector >= Geometry::Pi90 || ColatAngletoBisector <= (-1)*Geometry::Pi90)
    ColatAngletoBisector = 1.0/0.0;

  if ( ColatAngletoEntry >= Geometry::Pi90 )
    ColatAngletoEntry = 1.0/0.0;
  if ( ColatAngletoEntry <= (-1)*Geometry::Pi90 )
    ColatAngletoEntry = -1.0/0.0;
  if ( ColatAngletoExit >= Geometry::Pi90 )
    ColatAngletoExit = 1.0/0.0;
  if ( ColatAngletoExit <= (-1)*Geometry::Pi90 )
    ColatAngletoExit = -1.0/0.0;

  if ( D/R >= 1 ) {
    ColatAngletoEntry = -1.0/0.0;
    ColatAngletoExit = 1.0/0.0;
  }
  if ( D/R <= -1 ){
    ColatAngletoEntry = 1.0/0.0;
    ColatAngletoExit = -1.0/0.0;
    ColatAngletoBisector = -1.0/0.0;
    continuous = false;
  }


  retval = GCAD_RetVal(ColatAngletoEntry, ColatAngletoExit, ColatAngletoBisector, continuous);


  return retval;


}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  CLASS:  CylinderFace  (from CellFace)                       ****
// ****                                                              ****
//
//   Methods Defined Here:
//
//     ...
//

//////
// CONSTRUCTOR:  CylinderFace()
//
CylinderFace::CylinderFace (Real radius, MediumCell * powner) :
  CellFace( powner ),
  mRadius ( radius ),
  mRad2   ( radius*radius )
{
  assert(radius >= 0);
  std::cout << "~~~~> Wants CylinderFace at radius " << mRadius << "\n";
}

//////
// METHOD:  CylinderFace :: Normal()
//
//  Returns surface normal direction for the given point.
//
R3::XYZ CylinderFace::Normal(R3::XYZ loc) const {
  R3::XYZ squashed(loc.x(), loc.y(), 0);  // Keep just zonal coords
  return squashed.UnitElse(R3::XYZ(1,0,0));
}

//////
// METHOD:  CylinderFace :: GetDistanceAboveFace()
//
//   Returns the direct (shortest) distance from the face manifold to
//   the given point, with the return value using the following sign
//   convention:
//
//     > 0:  Point is "above" the CellFace, with "above" being the
//           side pointed to by the Normal vector.
//     = 0:  Point is "on" the surface.
//     < 0:  Point is "below" the surface, or "inside" of it.
//
Real CylinderFace::GetDistanceAboveFace(const R3::XYZ & loc) const {
  R3::XYZ squashed(loc.x(), loc.y(), 0);  // Keep just zonal coords
  return squashed.Mag() - mRadius;
}


//////
// METHOD:  PlaneFace :: LinearRayDistToExit()
//
//   Returns the directed (along a ray) distance from the starting location
//   ('loc') to the CellFace, under the presumption that we are crossing the
//   face in the "exiting" direction. (Else return special sentinal values.
//   See corresponding function in RCUCylinder for additional info.)
//
//   Returns the distance (in 3D) along a ray from a starting point to
//   the intersection of a cylinder (treated in 2D as a circle in the
//   XY plane centered on the origin).
//
// INPUTS:
//
//   o  loc   - Starting location;
//   o  dir   - Unit vector in direction of travel (caller responsible
//              for ensuring vector is unit-magnitude)
//
// THEORY OF OPERATION:
//
//   The distance travelled is the solution to a quadratic equation:
//
//       A*d^2 + B*d + C = 0
//
//   with A = D^2, B = 2*L.dot.D, and C = L^2 - R^2.  The vector D is
//   a two-component vector containing the X and Y direcion cosines, L
//   is a two-vector containing the X and Y location coordinates, and
//   R is the radius of the cylinder.
//
//   Note that A is always positive (or zero).  B is positive whenever
//   the direction trends radially outwards, otherwise is zero (when
//   completely tangential) or negative (trending inwards).  And C is
//   negative so long as the starting location is inside the cylinder
//   radius
//
// RETURNS:
//
//   There are as many as two solutions to the quadratic.  The greater
//   (most positive, or least negative) solution represents the one
//   where the directed ray is exiting the cylinder along its travel
//   path (typically this is ahead of the starting location), and the
//   lesser solution is the one where the ray is entering the circle.
//   If the starting location is inside the circle, there will be one
//   positive, and one negative solution. If outside, there may be two
//   positives, two negative, or no solution at all.
//
//   If the solutions exist, the function returns the greater (most
//   positive / least negative) solution.  If no solutions exist, the
//   function returns either +infinity (if the ray will travel forever
//   w/o exiting the cylinder, i.e. parallel and inside the cylinder)
//   or else it will return -infinity, as a signal that the ray never
//   was and never will be inside the cylinder.
//
Real CylinderFace::LinearRayDistToExit(const R3::XYZ & loc, const R3::XYZ & dir) const {

  Real A = dir.x()*dir.x() + dir.y()*dir.y();           // D^2
  Real C = loc.x()*loc.x() + loc.y()*loc.y() - mRad2;   // L^2 - R^2

  if (A==0) {             // First test for no-solution: case when
                          // motion is parallel to cylinder wall
    if (C <= 0) {
      return (1./0.);     // Returns +inf, a signal meaning raypath is
                          // parallel to cylinder and either inside or
    }                     // on-the-surface. (ie, it will never exit)
    else {
      return (-1./0.);    // Returns -inf, signaling that the ray is
                          // outside the cylinder and will never enter
    }                     // it.

  }

  Real B = 2 * (loc.x()*dir.x() + loc.y()*dir.y());   // 2*L.dot.D
  Real urad = B*B - 4*A*C;   // (urad: Under the RADical)

  if (urad < 0) {         // Second test for no-solution: outside
                          // cylinder and was never / will never enter
    return (-1./0.);      // Returns -inf.
  }

  Real margin = sqrt(urad);
  Real dist = (margin - B) / (2*A);   // {-B + sqrt(B^2-4AC)}/2A
                                      // (The "positive" root)
  return dist;

}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  CLASS:  SphereFace  (from CellFace)                         ****
// ****                                                              ****
//
//   Methods Defined Here:
//
//     ...
//

//////
// CONSTRUCTOR:  SphereFace()
//
SphereFace::SphereFace (Real radius, face_id_e toporbottom,
                        MediumCell * powner) :
  CellFace( powner ),
  mRadius ( toporbottom==F_TOP ? radius : -radius ),
  mRad2   ( radius*radius )
{
  assert(radius >= 0);
  assert(toporbottom == F_TOP || toporbottom == F_BOTTOM);
}

//////
// METHOD:  SphereFace :: Normal()
//
//  Returns surface normal direction for the given point using sign of
//  mRadius to determine whether "outward" us "up" or "down".
//
R3::XYZ SphereFace::Normal(R3::XYZ loc) const {
  return IsOutwardNormal() ? loc.UnitElse(R3::XYZ(0,0,1))
                           : loc.UnitElse(R3::XYZ(0,0,1)).Negative();
}

//////
// METHOD:  SphereFace :: GetDistanceAboveFace()
//
//   Returns the direct (shortest) distance from the face plane to
//   the given point, with the return value using the following sign
//   convention:
//
//     > 0:  Point is "above" the CellFace, with "above" being the
//           side pointed to by the Normal vector.
//     = 0:  Point is "on" the surface.
//     < 0:  Point is "below" the surface, or "inside" of it.
//
Real SphereFace::GetDistanceAboveFace(const R3::XYZ & loc) const {
  return (mRadius > 0) ? loc.Mag() - mRadius
                       : -(mRadius + loc.Mag());
}


//////
// METHOD:  SphereFace :: LinearRayDistToExit()
//
//   Return the distance to the point where a straight-line ray will exit (or
//   has exited, if value negative) the interior region of the Spherical
//   CellFace.
//
// RETURNS:
//
//    +inf  -   We are inside the face, and will never CROSS the face in an
//              exiting direction. (Though we may have exited previously and
//              re-entered.  Or we might "touch" the face in a tangent without
//              crossing.  Or we might be on the face but entering permanently)
//              This return value can only happen for inward-normal spherical
//              faces.
//     > 0  -   We will exit the face at some future time.  We are not necessarily
//              inside the face *right now*, but there is a definitive future exit.
//    == 0  -   We are on surface and exiting.
//     < 0  -   We have already exited and are outside the face. (Does not
//              preclude re-entereing, e.g. if we're in the "bubble" inside an
//              inward-normal sphere.)
//    -inf  -   We are outside the face and will never CROSS the face, nor have we
//              ever been definitively inside the face.  (Though it can happen that
//              we touch the face at a tangent point.)
//
// SOLUTION METHOD:
//
//   The return value is a solution to a quadratic equation, with certain
//   sentinal values communicated as infinities.  The quadratic equation has
//   two solutions, d(+) and d(-), though they may be imaginary or degenerate,
//   and these cases require special handling.  Degeneracy and existence are
//   determing by 'urad', the "value under the radical".  Return values depend
//   on whether we are an outward-normal face or an inward-normal face, and are
//   given by the following table:
//
//                               Out-Normal     In-Normal
//           urad < 0               -inf           +inf
//           urad == 0              -inf           +inf
//           urad > 0:
//                midpt <= 0        d(+)           +inf
//                midpt > 0         d(+)           d(-)
//
//   The equation:
//
//     d^2 + 2 loc.dot.dir d + loc^2 = R^2, with solutions:
//     d = -(loc.dot.dir) +/- sqrt( R^2 - loc^2 + (loc.dot.dir)^2 )
//
Real SphereFace::LinearRayDistToExit(const R3::XYZ & loc,
                                     const R3::XYZ & dir) const {

  Real midpt = -loc.Dot(dir);
  Real urad = mRad2 + midpt*midpt - loc.MagSquared();

  if (urad <= 0) { // No intersections or degenerate
    return IsOutwardNormal() ? (-1./0.) : (1./0.);
  }

  Real sqrad = sqrt(urad);
  Real dplus = midpt + sqrad;

  if (IsOutwardNormal()) return dplus;
  if (midpt <= 0) return (1./0.);

  Real dminus = midpt - sqrad;

  return dminus;

}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  CLASS:  GCAD_RetVal                                         ****
// ****                                                              ****
//

GCAD_RetVal::GCAD_RetVal(Real enter, Real exit, Real half, bool continuous) :
  mEntry(enter), mExit(exit), mHalf(half), mContinuous(continuous)
{}


GCAD_RetVal::GCAD_RetVal() :
  mEntry(0), mExit(0), mHalf(0), mContinuous(true)
{}


bool GCAD_RetVal::Inside (const Real & theta) const{
  Real error = 0.0000000001;
  if(mContinuous == true){
    if ( theta <= Exit() && theta >= (Entry()-error) )
      return true;
  }

  if(mContinuous == false){
    if ( (theta >= (-1)*Geometry::Pi90 && theta <= Exit() ) ||
         (theta >= (Entry()-error) && theta <= Geometry::Pi90) )
      return true;
  }

  return false;

}


bool GCAD_RetVal::IsProper (const GCAD_RetVal & a, const GCAD_RetVal & b, const GCAD_RetVal & c) const{

  if ( a.Inside(Exit()) == true &&
       b.Inside(Exit()) == true &&
       c.Inside(Exit()) == true )
    return true;
  else
    return false;

}
