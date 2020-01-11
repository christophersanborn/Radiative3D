// media_cellface.cpp
//
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
//     o   CellFace()
//     o   LinkTo()
//     o   VelocityJump()
//     o   GetRTBasis()
//     o   DistToPoint()
//     o   DistToExitFace()
//     o   GetCircArcDistToFace()
//

//////
// CONSTRUCTOR:  CellFace()
//
CellFace::CellFace (const R3::XYZ & N1, const R3::XYZ & N2,
                    const R3::XYZ & N3,
                    MediumCell * powner) :
  mCollect ( false  ),
  mReflect ( false  ),
  mAdjoin  ( false  ),
  mGridDiscon ( false ),
  mpOther  ( 0      ),
  mpCell   ( powner ),
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
// CONSTRUCTOR:  CellFace()
//
CellFace::CellFace (const R3::XYZ & N1, const R3::XYZ & N2,
                    const R3::XYZ & N3, const R3::XYZ & N4,
                    MediumCell * powner) :
  mCollect ( false  ),
  mReflect ( false  ),
  mAdjoin  ( false  ),
  mGridDiscon ( false ),
  mpOther  ( 0      ),
  mpCell   ( powner ),
  mPoint   ( N1     )

{
  // The three nodes, N1, N2, and N3, define a triangle that
  // represents the cell face.  The fourth node N4 is the excluded
  // node which is used to determine the direction of the normal vector.
  R3::XYZ Ray1 = N1.VectorTo(N2);
  R3::XYZ Ray2 = N1.VectorTo(N3);
  mNormal = Ray1.Cross(Ray2);
  mNormal.Normalize();

  if(DistToPoint(N4) > 0)             // Flips direction of mNormal if
    mNormal = mNormal.ScaledBy(-1);   // pointing inward

}


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

  RTCoef rt(Normal(), dir);     // We will return this object

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


//////
// METHOD:  CellFace :: DistToPoint()
//
//   Returns the direct (shortest) distance from the CellFace plane to
//   the given point, with the return value using the following sign
//   convention:
//
//     > 0:  Point is "above" the CellFace, with "above" being the
//           side pointed to by the mNormal vector.
//
//     = 0:  Point is "on" the surface.
//
//     < 0:  Point is "below" the surface.
//
//   Note: This function uses the opposite sign convention as the
//   related function, DirectedDistToPlane()
//
Real CellFace::DistToPoint(const R3::XYZ & loc) const {
  //
  return mNormal.Dot(mPoint.VectorTo(loc));
}


//////
// METHOD:  CellFace :: DistToExitFace()
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
Real CellFace::DistToExitFace(const R3::XYZ & loc,
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
// METHOD:  CellFace :: GetCircArcDistToFace()
//
//
//
GCAD_RetVal CellFace::GetCircArcDistToFace(const Real & R, const R3::XYZ & C, 
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
