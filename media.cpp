// media.cpp
//
#include <cmath>      /* sqrt, etc. */
#include <iostream>   /* std::cerr; remove when no longer needed */
#include <stdio.h> //remove
#include <algorithm>
#include "media.hpp"
#include "grid.hpp"
#include <iomanip>
// In this file:
// CLASS IMPLEMENTATIONS FOR:
//
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
PlaneFace RCUCylinder::cmLossFace(      // The CellFace through which
          R3::XYZ(0,0,1e10),            // phonons get lost through
          R3::XYZ(1,0,1e10),            // the cylinder sidewall.  Face
          R3::XYZ(0,1,1e10), 0          // isn't neede to correspond to
          );                            // actual geometry, just needed
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
//
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
//
//   Returns a read-write reference to either the top or bottom
//   CellFace object, as determined by the value of face_id
//
CellFace & RCUCylinder::Face(CellFace::face_id_e face_id) {
  switch (face_id) {
  case CellFace::F_TOP:
    return mTopFace;
  case CellFace::F_BOTTOM:
    return mBottomFace;
  default:
    std::cerr << "Invalid CellFace id.\n";
    throw std::exception();
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

  dist = mTopFace.GetDistanceAboveFace(loc);    // Get dist outside top end-cap
  if (dist > max_dist)                          //
    max_dist = dist;                            //

  dist = mBottomFace.GetDistanceAboveFace(loc); // Get dist outside bottom cap
  if (dist > max_dist)                          //
    max_dist = dist;                            //

  return max_dist;              // Return the max

}


//////
// METHOD:   RCUCylinder :: AdvanceLength()
//
// WHAT-IT-DOES:
//
//   Computes a travel record for the path travelled through the
//   MediumCell assuming the path length is known in advance.
//
//
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

  mFaces (PlaneFace( N2, N3, N4,    N1, this ),
          PlaneFace( N3, N4, N1,    N2, this ),
          PlaneFace( N4, N1, N2,    N3, this ),
          PlaneFace( N1, N2, N3,    N4, this ))

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

  Real distance = mFaces[0].GetDistanceAboveFace(loc);
  for (int i = 1; i < 4; i++){
    Real maybedistance = mFaces[i].GetDistanceAboveFace(loc);
    if(maybedistance > distance){
      distance = maybedistance;
    }
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
  rotAngle = (rotAngle > Geometry::Pi360) ? rotAngle-Geometry::Pi360 : rotAngle;

  //2D new location
  R3::XYZ newLoc2D = R3::XYZ(cos(rotAngle)*newLocRot2D.x() + sin(rotAngle)*newLocRot2D.z() , 0,
                            -sin(rotAngle)*newLocRot2D.x() + cos(rotAngle)*newLocRot2D.z());

  //Transform the new location from 2D to 3D space
  R3::XYZ newLoc = CT.RotMat().T() * (newLoc2D + CT.Trans());

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

  //Generate Travel Rec using function Advance Length
  TravelRec rec;
  rec = Tetra::AdvanceLength(rt,len,loc,dir);
  rec.pFace = &mFaces[faceID];

  return rec;
}


//////
// CONSTRUCTOR:  SphereShellD2()
//
SphereShellD2::SphereShellD2(Real RadTop, Real RadBot,
                             const GridData & DataTop, const GridData & DataBot) :
  mFaces( SphereFace(RadTop, CellFace::F_TOP, this),
          SphereFace(RadBot, CellFace::F_BOTTOM, this))
{

  std::cerr << "~~> Wants to make SphereShell > Inner: " << RadBot << "  Outer: " << RadTop << "\n";


}//
//

//////
// METHOD:   SphereShellD2 :: Face()
//
//   Returns a read-write reference to either the top or bottom
//   CellFace object, as determined by the value of face_id
//
CellFace & SphereShellD2::Face(CellFace::face_id_e face_id) {
  return mFaces[face_id];
}

//////
// METHOD:   SphereShellD2 :: GetVelocAtPoint()
// METHOD:   SphereShellD2 :: GetWavelengthAtPoint()
// METHOD:   SphereShellD2 :: GetDensityAtPoint()
// METHOD:   SphereShellD2 :: GetQatPoint()
//
Real SphereShellD2::GetVelocAtPoint(const R3::XYZ & loc, raytype type) const {
  throw std::runtime_error("UnimpSphereShellD2_GetVelocAtPoint");
  return  0;//loc.Dot(mVelGrad[type]) + mVel0[type];
}
//
Real SphereShellD2::GetWavelengthAtPoint(const R3::XYZ & loc, raytype type) const {
  throw std::runtime_error("UnimpSphereShellD2_GetWavelengthAtPoint");
  return  0;//(GetVelocAtPoint(loc,type)/ cmPhononFreq);
}
//
Real SphereShellD2::GetDensityAtPoint(const R3::XYZ & loc) const {
  throw std::runtime_error("UnimpSphereShellD2_GetDensityAtPoint");
  return  0;//loc.Dot(mDensGrad) + mDens0;
}
//
Real SphereShellD2::GetQatPoint(const R3::XYZ & loc, raytype type) const {
  throw std::runtime_error("UnimpSphereShellD2_GetQAtPoint");
  return 0;//mQ[type];
}

//////
// METHOD:   SphereShellD2 :: GetPathToBoundary()
//
TravelRec
SphereShellD2::GetPathToBoundary(raytype rt, const R3::XYZ & loc, const S2::ThetaPhi & dir) {
  throw std::runtime_error("UnimpSphereShellD2_GetPathToBoundary");
  return TravelRec(); // **FIXME**
}

//////
// METHOD:   SphereShellD2 :: AdvanceLength()
//
TravelRec
SphereShellD2::AdvanceLength(raytype rt, Real len, const R3::XYZ & startloc,
                     const S2::ThetaPhi & startdir) {
  throw std::runtime_error("UnimpSphereShellD2_AdvanceLength");
  return TravelRec(); // **FIXME**
}

//////
// METHOD:   SphereShellD2 :: IsPointInside()
//
Real SphereShellD2::IsPointInside(const R3::XYZ & loc) const {

  Real distance = mFaces[0].GetDistanceAboveFace(loc);
  for (int i = 1; i < 2; i++){
    Real maybedistance = mFaces[i].GetDistanceAboveFace(loc);
    if(maybedistance > distance){
      distance = maybedistance;
    }
  }
  return distance;

}
