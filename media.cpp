// media.cpp
//
#include <cmath>      /* sqrt, etc. */
#include <iostream>   /* std::cerr; remove when no longer needed */
#include <stdio.h> //remove
#include <algorithm>
#include "media.hpp"
#include "grid.hpp"
#include "rtcoef.hpp"
#include <iomanip>
// In this file:
// CLASS IMPLEMENTATIONS FOR:
//
//   o  Class CellFace
//   o  Class MediumCell
//   o  Class RCUCylinder
//   o  Class Tetra
//
// Search on "&&&&" to jump between class implementations in this
// file.
//

// HELPER FUNCTION PROTOTYPES:
//
Real RayIntersectCircleXY(R3::XYZ, R3::XYZ, Real);
                      // Used in RCUCylinder::GetPathToBoundary() 


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
///
///   Called by the model-building routines, this function handles the
///   book-keeping of linking together the two CellFaces (one from each
///   MediumCell of an adjoining pair) that represent the interface
///   between two model cells.
///
///   Takes a reference to the other CellFace, and a boolean argument
///   which signals whether the Grid specified that there should be a
///   first-order discontinuity across this interface.  The caller is
///   responsible for determining whether this is the case, and
///   communicating it via this argument.
///
///   @callergraph
///
void CellFace::LinkTo(CellFace & other, bool disc) {
  mpOther = &other;
  mAdjoin = true;
  mGridDiscon = disc;
  other.mpOther = this;
  other.mAdjoin = true;
  other.mGridDiscon = disc;
}//
///   If called with the following form, then assumption is that a
///   grid-indicated discontinuity across the interface has already
///   been flagged by a call to SetDiscontinuous() on one or both of
///   the faces.  (This could happen, e.g., at an earlier stage of the
///   model build-out process.)  Check mGridDiscon on both faces, and
///   if either is true, then explicitly set both to true to ensure
///   both faces agree.
///
///   @callergraph
///
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
///
///   Returns a number characterizing the size of the velocity
///   discontinuity across the CellFace (from one MediumCell to
///   another) at a given location.  Computes the number as a
///   fractional velocity change. Returns a non-negative number, and
///   returns the largest of the fractional changes (P or S) at the
///   given location.
///
///   @callergraph
///
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
///
///   Populates an RTCoef object with a set of basis vectors in which
///   reflected or transmitted rays can be expanded; populates the
///   velocity and density parameters for both sides of the interface.
///
///   The information generated by this function and stored in the
///   RTCoef object are subsequently used (by member functions of the
///   RTCoef class or in other locations) to compute reflection and
///   transmission probabilities, and to choose an outcome.  The info
///   generated here is basically just to get the ball rolling on that
///   front.
///
///   The information generated here depends on the incidence angle
///   w.r.t. the CellFace, and needs access to MediumCells of which this
///   CellFace object has exclusive knowledge, and that's why this method
///   appears in the CellFace class, instead of in the RTCoef class or
///   elsewhere.
///
///   @callergraph
///
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
///
///   Returns the direct (shortest) distance from the CellFace plane to
///   the given point, with the return value using the following sign
///   convention:
///
///     > 0:  Point is "above" the CellFace, with "above" being the
///           side pointed to by the mNormal vector.
///
///     = 0:  Point is "on" the surface.
///
///     < 0:  Point is "below" the surface.
///
///   Note: This function uses the opposite sign convention as the
///   related function, DirectedDistToPlane()
///
///   @callergraph
///
Real CellFace::DistToPoint(const R3::XYZ & loc) const {
  //
  return mNormal.Dot(mPoint.VectorTo(loc));
}


//////
// METHOD:  CellFace :: DistToExitFace()
///
///   Returns the directed (along a ray) distance from the starting
///   location ('loc') to the CellFace, under the presumption that we
///   are crossing the face in the "exiting" direction.  Should it be
///   that the direction of our ray will cause us to cross the face in
///   "entering" direction, or else to not cross the face at all
///   (parallel motion), then special values are returned to signal
///   this condition.
///
/// INPUTS:
///
///   o  loc     -   The starting location of the ray
///   o  dir     -   The direction of the ray (unit vector -
///                           responsibility lies with caller
///                           to provide unit-magnitude vector)
/// RETURNS:
///
///    +inf  -   Means we will NEVER cross the face in exiting fashion,
///              and our future lies INSIDE the face (ie, we are
///              entering or have entered, or are running parallel and
///              are already inside (including on-surface))
///
///     > 0  -   We will exit the face at some future time.
///
///    == 0  -   We are exiting right now (loc is on surface and dir
///              points outward)
///
///     < 0  -   We have already exited the face, at some point in the
///              past.
///
///    -inf  -   We will NEVER exit the face, because we are already
///              outside, and running parallel to the face.
///
///   @callergraph
///
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
///
///   @callergraph
///   @callgraph
///  
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
// ****  CLASS:  MediumCell                                          ****
// ****                                                              ****
//

//
// Static Member Initialization:  (MediumCell Class)
//

Real MediumCell::cmPhononFreq = 1.0;    // Responsibility to set at
                                        // runtime lies with Model
                                        // constructor


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  CLASS:  RCUCylinder                                         ****
// ****                                                              ****
//
//   Methods Defined Here:
//
//     o   RCUCylinder()
//     o   GetVelocAtPoint()
//     o   GetWavelengthAtPoint()
//     o   GetDensityAtPoint()
//     o   GetQatPoint()
//     o   Face()
//     o   IsPointInside()
//     o   AdvanceLength()
//     o   GetPathToBoundary()
//

//
// Static Member Initialization:  (RCUCylinder Class)
//

Real RCUCylinder::cmRange = 1200.0;     // Responsibility to set at
bool RCUCylinder::cmRangeSet = false;   // runtime lies with Model
                                        // constructor.
CellFace RCUCylinder::cmLossFace(       // The CellFace through which
         R3::XYZ(0,0,1e10),             // phonons get lost through
         R3::XYZ(1,0,1e10),             // the cylinder sidewall.  Face
         R3::XYZ(0,1,1e10), 0           // isn't neede to correspond to
         );                             // actual geometry, just needed
                                        // for bookkeeping purposes.

//////
// CONSTRUCTOR:  RCUCylinder()
//
RCUCylinder::RCUCylinder(
              R3::XYZ N1_top, R3::XYZ N2_top, R3::XYZ N3_top,
              R3::XYZ N1_bot, R3::XYZ N2_bot, R3::XYZ N3_bot,
              Elastic::Velocity vpvs_top,
              Elastic::Velocity vpvs_bot,
              Real rho_top,  Real rho_bot,
              Elastic::Q q_top, Elastic::Q q_bot)
:
  mTopFace    ( N1_top, N2_top, N3_top, this ), // Makes a face pointing up.
  mBottomFace ( N1_bot, N3_bot, N2_bot, this )  // Makes a face pointing down.
{

  if (cmRangeSet==false)  // (Ensure static vars have been initialized)
    {throw(std::logic_error("RCUCylinder Range unset."));} 

  mVelTop[RAY_P] = vpvs_top.Vp();   // Seismic Velocities
  mVelTop[RAY_S] = vpvs_top.Vs();
  mVelBot[RAY_P] = vpvs_bot.Vp();
  mVelBot[RAY_S] = vpvs_bot.Vs();
  mDensity = rho_top;               // Density, arbitrary units
  mQ[RAY_P] = q_top.Qp(vpvs_top);
  mQ[RAY_S] = q_top.Qs(vpvs_top);

}


//////
// METHOD:   RCUCylinder :: GetVelocAtPoint()
///
///  Gets elastic velocity at location within model cell.
///
Real RCUCylinder::GetVelocAtPoint(const R3::XYZ & loc, 
                                  raytype type) const {
  return mVelTop[type]; // (Pretty simple in a constant-velocity
                        //  medium...)
}


//////
// METHOD:   RCUCylinder :: GetWavelengthAtPoint()
//
Real RCUCylinder::GetWavelengthAtPoint(const R3::XYZ & loc, 
                                       raytype type) const {
  return mVelTop[type] / cmPhononFreq;

}


//////
// METHOD:   RCUCylinder :: GetDensityAtPoint()
//
Real RCUCylinder::GetDensityAtPoint(const R3::XYZ & loc) const {

  return mDensity;

}

/////
// METHOD:   RCUCylinder :: GetQatPoint()
//
Real RCUCylinder::GetQatPoint(const R3::XYZ & loc, raytype type) const {

  return mQ[type];

}


//////
// METHOD:   RCUCylinder :: Face()
///
///   Returns a read-write reference to either the top or bottom
///   CellFace object, as determined by the value of face_id
///
CellFace & RCUCylinder::Face(CellFace::face_id_e face_id) {
  switch (face_id) {
  case CellFace::F_TOP:
    return mTopFace;
  case CellFace::F_BOTTOM:
    return mBottomFace;
  default:
    /* TODO: Need Proper error-handling here */ 
    *((int*)0) = 1; // force segfault, kindof a brutal exit
    return *((CellFace*)0); // Bogus return - suppresses compile warning
  }
}


//////
// METHOD:   RCUCylinder :: IsPointInside()
//
//   Reports whether a given point is "inside" the Cylinder cell by
//   returning a "mismatch" factor.
//
//   If the return value is:
//
//     <= 0:  then the point is definitively inside or on-the-boundary
//            of the cell.
//
//      > 0:  then the point is outside the cell.  The numeric value
//            represents the direct-distance to the nearest cell-wall.
//
// THEORY of operation:
//
//   A point is only inside the cell if it is inside ALL of the cell's
//   bounding surfaces.  Thus for each bounding surface (two planes
//   and a cylinder) we compute the distance "above" the surface (such
//   that a positive value means "outside" the surface and negative
//   means inside), and thus, among these three distances we report
//   the "greatest" value as the mismatch, which roughly equates to
//   the distance inside (if negative) or outside (if positive) the
//   cell.
//
Real RCUCylinder::IsPointInside(const R3::XYZ & loc) const {

  Real max_dist;

  Real dist                     // Distance outside the 
    = sqrt(loc.x()*loc.x()      // cylinder wall
           + loc.y()*loc.y())   //
      - cmRange;                //

  max_dist = dist;              // (So far the greatest)

  dist = mTopFace.DistToPoint(loc);     // Get dist outside top end-cap
  if (dist > max_dist)                  //
    max_dist = dist;                    //
 
  dist = mBottomFace.DistToPoint(loc);  // Get dist outside bottom cap
  if (dist > max_dist)                  //
    max_dist = dist;                    //

  return max_dist;              // Return the max

}


//////
// METHOD:   RCUCylinder :: AdvanceLength()
//
// WHAT-IT-DOES:
///
///   Computes a travel record for the path travelled through the
///   MediumCell assuming the path length is known in advance.
///
///
TravelRec
RCUCylinder::AdvanceLength(raytype rt, Real len,
                           const R3::XYZ & startloc,
                           const S2::ThetaPhi & startdir) {
  TravelRec rec;

  rec.PathLength = len;
  rec.TravelTime = len / mVelTop[rt];
  rec.NewLoc     = startloc + R3::XYZ(startdir).ScaledBy(len);
  rec.NewDir     = startdir;
  rec.Attenuation = HelperAttenuation(rec.TravelTime * cmPhononFreq, mQ[rt]);

  return rec;

}


//////
// METHOD:   RCUCylinder :: GetPathToBoundary()
//
// WHAT-IT-DOES:
//
///   Computes a travel record for a path travelled through the
///   MediumCell from the starting location/direction until a cell
///   boundary is hit.  Travel record encodes the length/time
///   travelled, and location of boundary intersection and direction of
///   raypath tangent at that point.
//
// RETURNS:
//
//   o  Distance travelled, as the actual return value
//   o  The TravelRec, (received by pointer and modified)
//
TravelRec 
RCUCylinder::GetPathToBoundary(raytype rt,
                               const R3::XYZ & loc, 
                               const S2::ThetaPhi & dir) {
  TravelRec rec;
  enum exitface_e {TOP, BOTTOM, LOSS, NUM_EF};
  Real dists[NUM_EF];

  // :::
  // :: Get distances to cylinder sidewall and top/bottom faces:
  // :
  //                Note: negative distances imply that we have
  //                already exited through the surface. (They most
  //                likely result from slight numberical errors.) We
  //                squash these negative numbers to zero to avoid
  //                engaging in retrograde motion.  The zero should
  //                result in immediate handover (at the caller level)
  //                to the neighboring cell, which presumably will be
  //                the definitive home of the point in question (or
  //                at least stands a good chance of it).
  //

  dists[LOSS]   = RayIntersectCircleXY(loc, dir, cmRange);
  if (dists[LOSS] < 0) dists[LOSS] = 0; // (exclude negatives)
                    //
                    //   >0  : Implies an ordinary path to the cylinder.
                    //
                    //  ==0  : Implies on the surface and exiting or
                    //         outside surface and not entering
                    //
                    //   <0  : (We should never see this value)
                    //
                    //  +inf : Running parallel to cylinder and not
                    //         outside it.
                    //

  dists[TOP]    = mTopFace.DistToExitFace(loc, dir);
  dists[BOTTOM] = mBottomFace.DistToExitFace(loc, dir);
  if (dists[TOP] < 0)  dists[TOP] = 0;       // Squash negatives to 
  if (dists[BOTTOM] < 0)  dists[BOTTOM] = 0; // zero.

                    //   >0  : Implies positive (ordinary) path to
                    //         face.
                    //
                    //  ==0  : On face and exiting or else ouside face
                    //         and not entering
                    //
                    //   <0  : (We should never see this)
                    //
                    //  +inf : Entering through face, or else running
                    //         parallel either on or inside face.
                    //         Interpret as "We'll hit another face
                    //         before we hit this one."
                    //

  // ::::::
  // :: Determine exit face: (Phase 1)
  // :
  //                We default to Cylinder wall (LOSS) exit.  But
  //                check whether Top or Bottom face have shorter
  //                distances.
  //

  exitface_e  exf = LOSS;         // The default
  Real shortest = dists[LOSS];    //
  if (dists[TOP] < shortest) {
    exf      = TOP;
    shortest = dists[TOP];
  }
  if (dists[BOTTOM] < shortest) {
    exf      = BOTTOM;
    shortest = dists[BOTTOM];
  }

  
  // ::::::
  // :: Count the Zeros: (Phase 2)
  // :
  //                If >1 dists are EXACTLY zero, then we are exiting
  //                through a corner. Ditch the phonon rather than
  //                deal with the complexity of figuring out into
  //                which cell to inject it.
  //

  unsigned zero_count = 0;
  for (unsigned i = 0; i < NUM_EF; i++) {
    if (dists[i]==0) zero_count++;
  }

  if (zero_count > 1) {
    exf      = LOSS;
    shortest = dists[LOSS];
  }
  
  // ::::::
  // :: Populate the TravelRec and return it to the caller:
  // :

  rec.PathLength  = shortest;
  rec.TravelTime  = shortest / mVelTop[rt];
  rec.NewLoc      = loc + R3::XYZ(dir).ScaledBy(shortest);
  rec.NewDir      = dir;
  rec.Attenuation = HelperAttenuation(rec.TravelTime * cmPhononFreq, mQ[rt]);

  switch (exf) {
  case TOP:
    rec.pFace = &mTopFace;
    break;
  case BOTTOM:
    rec.pFace = &mBottomFace;
    break;
  default:
    rec.pFace = &cmLossFace;
  }

  return rec;

}


//////
// METHOD:   RCUCylinder :: HelperAttenuation()
//
//   Computes the amplitude attenuation factor given the number of
//   wave cycles (traveltime * frequency) and the Q factor in effect
//   throughout the travel arc (assumes Q is a uniform quantity).
//
//   Returns the exponential factor in the amplitude decay eqation:
//
//     A(t) = A0 exp( -omega*t / 2Q )
//
//   Note: (1) The factor of 2 in the denominator is because we are
//   computing amplitude decay. When amplitude is squared to get
//   energy, the factor of 2 goes away.  (2) Omega = 2*pi*frequency,
//   and since we're given cycles, the factor of 2*pi must appear in
//   the numerator of our exponent. The 2 in the numerator and the 2
//   in the denominator cancel out, however.
//
Real RCUCylinder::HelperAttenuation(Real cycles, Real Q) {

  return exp( ((-1)*Geometry::Pi * cycles) /
              (Q)
            );

}


//////
// HELPER FUNCTION:  RayIntersectCircleXY()
//
//   Returns the distance (in 3D) along a ray from a starting point to
//   the intersection of a cylinder (treated in 2D as a circle in the
//   XY plane centered on the origin).
//
//
// INPUTS:
//
//   o  loc   - Starting location;
//   o  dir   - Unit vector in direction of travel (caller responsible 
//              for ensuring vector is unit-magnitude)
//   o  rad   - Radius of the circle/cylinder
//
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
Real RayIntersectCircleXY(R3::XYZ loc, R3::XYZ dir, Real rad) {

  Real A = dir.x()*dir.x() + dir.y()*dir.y();           // D^2
  Real C = loc.x()*loc.x() + loc.y()*loc.y() - rad*rad; // L^2 - R^2

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
// ****  CLASS:  Tetra                                               ****
// ****                                                              ****
//
//   Methods Defined Here:
//
//     o   Tetra()
//     o   GetVelocAtPoint()
//     o   GetWavelengthAtPoint()
//     o   GetDensityAtPoint()
//     o   GetQatPoint()
//     o   Face()
//     o   IsPointInside()
//     o   AdvanceLength()
//     o   GetPathToBoundary()
//

//////
// CONSTRUCTOR:  Tetrahedral()
//
Tetra::Tetra(R3::XYZ N1, R3::XYZ N2, R3::XYZ N3, R3::XYZ N4, 
             const GridData & dataA, const GridData & dataB, 
             const GridData & dataC, const GridData & dataD): 

  mFaces (CellFace( N2, N3, N4,    N1, this ),
          CellFace( N3, N4, N1,    N2, this ),
          CellFace( N4, N1, N2,    N3, this ),
          CellFace( N1, N2, N3,    N4, this ))

  {
    //
    // Temporary Matrix and Columns 
    // used to compute gradients
    //
    R4::Rows LocMatrix = R4::Rows(N1, 1, N2, 1, N3, 1, N4, 1);

    R4::Column DensCol = R4::Column(dataA.Rho(), dataB.Rho(), 
                                    dataC.Rho(), dataD.Rho());

    R4::Column VelColP = R4::Column(dataA.Vp(), dataB.Vp(), 
                                    dataC.Vp(), dataD.Vp());

    R4::Column VelColS = R4::Column(dataA.Vs(), dataB.Vs(), 
                                    dataC.Vs(), dataD.Vs());

    R4::Column VelCoefP = LocMatrix.SolveAXB(VelColP);
    R4::Column VelCoefS = LocMatrix.SolveAXB(VelColS);
    R4::Column DensCoef = LocMatrix.SolveAXB(DensCol);
    
    /////////////////////////
    // Tetra Class Members //
    /////////////////////////
    mVelGrad[RAY_P] = VelCoefP.TruncXYZ();  // P Velocity Gradient
    mVelGrad[RAY_S] = VelCoefS.TruncXYZ();  // S Velocity Gradient
    mDensGrad       = DensCoef.TruncXYZ();  // Density Velocity Gradient
    mVel0[RAY_P]    = VelCoefP.x(3);        // P Velocity at (0,0,0)
    mVel0[RAY_S]    = VelCoefS.x(3);        // S Velocity at (0,0,0)
    mDens0          = DensCoef.x(3);        // Density at (0,0,0)
    mQ[RAY_P]       = (dataA.Qp() + dataB.Qp()        // TODO: needs 
                        + dataC.Qp() + dataD.Qp())/4; // geometric center weights
    mQ[RAY_S]       = (dataA.Qs() + dataB.Qs()        // TODO: needs 
                        + dataC.Qs() + dataD.Qs())/4; // geometric center weights
  }


//////
// METHOD:   Tetra :: GetVelocAtPoint()
//
Real Tetra::GetVelocAtPoint(const R3::XYZ & loc, 
                            raytype type) const {

  return  loc.Dot(mVelGrad[type]) + mVel0[type];
  
}


//////
// METHOD:   Tetra :: GetWavelengthAtPoint()
//
Real Tetra::GetWavelengthAtPoint(const R3::XYZ & loc, 
                                 raytype type) const {
  
  
  return  (GetVelocAtPoint(loc,type)/ cmPhononFreq);
}


//////
// METHOD:   Tetra :: GetDensityAtPoint()
//
Real Tetra::GetDensityAtPoint(const R3::XYZ & loc) const {

  return  loc.Dot(mDensGrad) + mDens0;

}

/////
// METHOD:   Tetra :: GetQatPoint()
//
Real Tetra::GetQatPoint(const R3::XYZ & loc, raytype type) const {

  return mQ[type];

}


//////
// METHOD:   Tetra :: Face()
//
//   Returns a read-write reference to either the top or bottom
//   CellFace object, as determined by the value of face_id
//
CellFace & Tetra::Face(CellFace::face_id_e face_id) {
  return mFaces[face_id];
}


//////
// METHOD:   Tetra :: IsPointInside()
//
//   Reports whether a given point is "inside" the Tetrahedral cell by
//   returning a "mismatch" factor.
//
//   If the return value is:
//
//     <= 0:  then the point is definitively inside or on-the-boundary
//            of the cell.
//
//      > 0:  then the point is outside the cell.  The numeric value
//            represents the direct-distance to the nearest cell-wall.
//
// THEORY of operation:
//
//   A point is only inside the cell if it is inside ALL of the cell's
//   bounding surfaces.  Thus for each bounding surface we compute the
//   distance "above" the surface (such that a positive value means
//   "outside" the surface and negative means inside), and thus, among
//   these three distances we report the "greatest" value as the
//   mismatch, which roughly equates to the distance inside (if
//   negative) or outside (if positive) the cell.
//
Real Tetra::IsPointInside(const R3::XYZ & loc) const {

  Real distance = mFaces[0].DistToPoint(loc);
  
  for (int i = 1; i < 4; i++){
   
    if(mFaces[i].DistToPoint(loc) > distance){
      
      distance = mFaces[i].DistToPoint(loc);}
   
  }

  return distance;
    

}


//////
// METHOD:   Tetra :: AdvanceLength()
//
// WHAT-IT-DOES:
//
//   Computes a travel record for the path travelled through the
//   MediumCell assuming the path length is known in advance.
//
//
// THEORY of operation:
//
//  A phonon travels along a circular arc in the t-g plane
//  where t is the tangential direction of the phonon and g is
//  the of the cell. Using a coordinate transformation and a given
//  length, the new location and new direction are calculated in this
//  2D plane.
TravelRec
Tetra::AdvanceLength(raytype rt, Real len,
                     const R3::XYZ & startloc,
                     const S2::ThetaPhi & startdir) {
 
  CoordinateTransformation CT =  CoordinateTransformation(GetVelocAtPoint(startloc, rt),mVelGrad[rt],
                                                          startloc, startdir.XYZ());
  
  Real theta = len/CT.Radius();  //Angle separating 2D current location
                                 //and 2D exit.

  //////////////////////////////
  /// Determine New Location ///
  //////////////////////////////

  //2D new location in rotated system
  R3::XYZ newLocRot2D = R3::XYZ(CT.Radius()*sin(theta/2),0,CT.Radius()*cos(theta/2));

  //rotation angle to return to prime 2D system
  Real angletoX0 = atan2(CT.PrimeLoc().x(),CT.PrimeLoc().z());
  Real rotAngle = angletoX0 + (theta/2);
  rotAngle = (rotAngle > Geometry::Pi360) ? rotAngle-=Geometry::Pi360 : rotAngle;
 
  //2D new location
  R3::XYZ newLoc2D = R3::XYZ(cos(rotAngle)*newLocRot2D.x() + sin(rotAngle)*newLocRot2D.z() , 0,
                            -sin(rotAngle)*newLocRot2D.x() + cos(rotAngle)*newLocRot2D.z());
  
  //Transform the new location from 2D to 3D space
  R3::XYZ newLoc = CT.RotMat().T() * (newLoc2D + CT.Trans());
 
 
  if(IsPointInside(newLoc) > 0.1){
     std::cerr << "new loc " << IsPointInside(newLoc) << " km outside cell\n";
     std::cerr << "length: " << len << std::endl;
     std::cerr << "( " << newLoc.x() << "   " << newLoc.y() << "   " << newLoc.z() << "   )" << std::endl;
     exit(1);
  }

  ///////////////////////////////
  /// Determine New Direction ///
  ///////////////////////////////
  
  //New Direction in 2D plane is the 
  //tangential direction at the new location
  Real AngleNewLoc2D = atan2(newLoc2D.x(),newLoc2D.z());
  R3::XYZ newDir2D = R3::XYZ(cos(AngleNewLoc2D),0,(-1)*sin(AngleNewLoc2D));

  //Transform the new direction from 2D to 3D space
  R3::XYZ newDir = CT.RotMat().T()*newDir2D;
  newDir.Normalize();
  
  
  Real travelTime =(1/mVelGrad[rt].Mag())*( log(fabs(tan((AngleNewLoc2D/2 + Geometry::Pi45))))
                                            - log(fabs(tan((angletoX0/2 + Geometry::Pi45)))));
  
  //Generate Travel Rec using new direction and new location
  TravelRec rec;
 
  rec.PathLength = len;
  rec.TravelTime = travelTime;
  rec.NewLoc = newLoc;
  rec.NewDir = S2::Node(newDir.x(),newDir.y(),newDir.z());
  rec.Attenuation = exp( ((-1) * Geometry::Pi * rec.TravelTime * cmPhononFreq) / mQ[rt] );
  
  return rec;
}


//////
// METHOD:   Tetra :: GetPathToBoundary()
//
// WHAT-IT-DOES:
//
//   Computes a travel record for a path travelled through the
//   MediumCell from the starting location/direction until a cell
//   boundary is hit.  Travel record encodes the length/time
//   travelled, and location of boundary intersection and direction of
//   raypath tangent at that point.
//
// RETURNS:
//
//   o  Distance travelled, as the actual return value
//   o  The TravelRec, (received by pointer and modified)
//
TravelRec 
Tetra::GetPathToBoundary(raytype rt,
                         const R3::XYZ & loc, 
                         const S2::ThetaPhi & dir){
  if(IsPointInside(loc) > 0.1){
    std::cerr << "IsPointInside: " << IsPointInside(loc) << std::endl;
    exit(1);
  }
 
  CoordinateTransformation CT =  CoordinateTransformation(GetVelocAtPoint(loc, rt),
                                                          mVelGrad[rt],loc,dir.XYZ());

  //Determine azimuthal angle to current location
  //in transformed 2D system
  Real ColatAngletoX0 = atan2(CT.PrimeLoc().x(),CT.PrimeLoc().z());  

  //Compute entry and exit angles for each cell face.
  //Angles stored in an XYZ object:
  // X = Entry into face
  // Y = Exit out of face
  // Z = Bisector (+/- Pi to orient with cell face normal) 
  GCAD_RetVal retvals[4];

  for(int i = 0; i < 4; i++)
    {
      retvals[i] = mFaces[i].GetCircArcDistToFace(CT.Radius(),CT.Trans(),CT.RotMat());
    }

  Real len = 1.0/0.0;
  int faceID = 0;
 
  for (int i=0; i < 4; i++){
    if( retvals[i].IsProper(retvals[(i+1)%4],retvals[(i+2)%4],retvals[(i+3)%4]) == true ){
      
      Real newlen = (retvals[i].Exit() - ColatAngletoX0) * CT.Radius();
      
      if ( newlen < 0 && (ColatAngletoX0 > retvals[i].Half()) ){
        newlen = len; //dismiss exit
      }
      
      if( newlen < len){
        len = newlen;
        faceID = i;
      }
    }
  }
  /*
  std::cerr << faceID << std::endl;
  std::cerr << len << std::endl;

  std::cerr << std::setprecision(15);
  
  std::cerr << "xo: " << ColatAngletoX0 << std::endl;
  for(int i=0; i<4; i++){
  std::cerr << "entry: " << retvals[i].Entry() << "  exit: " << retvals[i].Exit() 
            << " isproper: " <<  retvals[i].IsProper(retvals[(i+1)%4],retvals[(i+2)%4],retvals[(i+3)%4])
            << " exit-xo: " << retvals[i].Exit() - ColatAngletoX0 << std::endl;
}
  */
  


  //Generate Travel Rec using function Advance Length
  TravelRec rec;
  rec = Tetra::AdvanceLength(rt,len,loc,dir);
  rec.pFace = &mFaces[faceID];

  return rec;
}


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
 
