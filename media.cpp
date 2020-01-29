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
//   o  Class SphereShell
//
// Search on "&&&&" to jump between class implementations in this
// file.
//

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


//////
// METHOD:   MediumCell :: IsPointInside()
//
//   Reports whether a given point is "inside" the Cylinder cell by returning
//   a "mismatch" factor.
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
//   A point is only inside the cell if it is inside ALL of the cell's bounding
//   surfaces.  Thus for each bounding surface (two planes and a cylinder) we
//   compute the distance "above" the surface (such that a positive value means
//   "outside" the surface and negative means inside), and thus, among these
//   three distances we report the "greatest" value as the mismatch, which
//   roughly equates to the distance inside (if negative) or outside (if
//   positive) the cell.
//
Real MediumCell::IsPointInside(const R3::XYZ & loc) const {

  auto x_this = (MediumCell *)this; // Promise not modify...

  Count nFaces = NumFaces();
  Real distance = x_this->Face(0).GetDistanceAboveFace(loc);

  for (Index i = 1; i < nFaces; i++){
    Real maybe_distance = x_this->Face(i).GetDistanceAboveFace(loc);
    if(maybe_distance > distance){
      distance = maybe_distance;
    }
  }

  return distance;

}


//////
// METHOD:   MediumCell :: HelperUniformAttenuation()    (static)
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
Real MediumCell::HelperUniformAttenuation(Real cycles, Real Q) {
  return exp( ((-1)*Geometry::Pi * cycles) / (Q) );
}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  CLASS:  RCUCylinder                                         ****
// ****                                                              ****
//
//   Methods Defined Here:
//
//     o   RCUCylinder()
//     o   Face()
//     o   GetVelocAtPoint()
//     o   GetWavelengthAtPoint()
//     o   GetDensityAtPoint()
//     o   GetQatPoint()
//     o   AdvanceLength()
//     o   GetPathToBoundary()
//

//
// Static Member Initialization:  (RCUCylinder Class)
//

CylinderFace RCUCylinder::cmLossFace  // The CellFace through which phonons get
              (1, nullptr);           // lost through the cylinder sidewall. We
                                      // defer setting reasonable radius, and no
                                      // MediumCell needs to be associated.
bool RCUCylinder::cmRangeSet = false; // Responsibility to set radius at runtime
                                      // lies with Model constructor via call to
                                      // RCUCylinder::SetRange().

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
// METHOD:   RCUCylinder :: Face()
//
//   Returns a read-write reference to either the top or bottom
//   CellFace object, as determined by the value of face_id
//
CellFace & RCUCylinder::Face(Index face_id) {
  switch (face_id) {
  case CellFace::F_TOP:
    return mTopFace;
  case CellFace::F_BOTTOM:
    return mBottomFace;
  case CellFace::F_SIDE:
    return cmLossFace;
  default:
    throw Invalid("Invalid CellFace ID for RCUCylinder medium cell.");
  }
}

//////
// METHODS:  RCUCylinder :: GetVelocAtPoint()
//           RCUCylinder :: GetWavelengthAtPoint()
//           RCUCylinder :: GetDensityAtPoint()
//           RCUCylinder :: GetQatPoint()
//
Real RCUCylinder::GetVelocAtPoint(const R3::XYZ & loc, raytype type) const {
  return mVelTop[type];
}
Real RCUCylinder::GetWavelengthAtPoint(const R3::XYZ & loc, raytype type) const {
  return mVelTop[type] / cmPhononFreq;
}
Real RCUCylinder::GetDensityAtPoint(const R3::XYZ & loc) const {
  return mDensity;
}
Real RCUCylinder::GetQatPoint(const R3::XYZ & loc, raytype type) const {
  return mQ[type];
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
  rec.Attenuation = HelperUniformAttenuation(rec.TravelTime * cmPhononFreq, mQ[rt]);

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
TravelRec RCUCylinder::GetPathToBoundary(raytype rt, const R3::XYZ & loc,
                                                     const S2::ThetaPhi & dir) {
  TravelRec rec;
  enum exitface_e {TOP=CellFace::F_TOP,
                   BOTTOM=CellFace::F_BOTTOM,
                   LOSS=CellFace::F_SIDE,
                   NUM_FACES=3};
  Real dists[NUM_FACES];

  // :::
  // :: Get distances to cylinder sidewall and top/bottom faces:
  // :
  //        Note: negative distances imply that we have already exited
  //        through the surface. (They most likely result from slight
  //        numerical errors.) We squash these negative numbers to zero to
  //        avoid engaging in retrograde motion.  The zero should result in
  //        immediate handover (at the caller level) to the neighboring cell,
  //        which presumably will be the definitive home of the point in
  //        question (or at least stands a good chance of it).
  //

  dists[LOSS]   = cmLossFace.LinearRayDistToExit(loc, dir);
  dists[TOP]    = mTopFace.LinearRayDistToExit(loc, dir);
  dists[BOTTOM] = mBottomFace.LinearRayDistToExit(loc, dir);
  if (dists[LOSS]   < 0) dists[LOSS] = 0;      // Squash negatives to
  if (dists[TOP]    < 0) dists[TOP] = 0;       // zero.
  if (dists[BOTTOM] < 0) dists[BOTTOM] = 0;    //
        //
        //   >0  : Implies ordinary forward path to the cylinder.
        //  ==0  : Implies on the surface and exiting or outside
        //         surface and not entering
        //   <0  : (We should never see this value due to squash)
        //  +inf : Entering through face, or else running parallel either on
        //         or inside face.  Interpret as "We'll hit another face
        //         before we hit this one."

  // ::::::
  // :: Determine exit face: (Phase 1)
  // :
  //        We default to Cylinder wall (LOSS) exit.  But check
  //        whether Top or Bottom face have shorter distances.
  //

  exitface_e  exf = LOSS;         // The default,
  Real shortest = dists[LOSS];    //
  if (dists[TOP] < shortest) {    // But maybe we hit these
    exf = TOP;                    // sooner...
    shortest = dists[TOP];
  }
  if (dists[BOTTOM] < shortest) {
    exf      = BOTTOM;
    shortest = dists[BOTTOM];
  }

  // ::::::
  // :: Count the Zeros: (Phase 2)
  // :
  //        If >1 dists are EXACTLY zero, then we are exiting through a
  //        corner. Ditch the phonon rather than deal with the complexity of
  //        figuring out into which cell to inject it.
  //
  //        UPDATE: This situation is so incredibly improbable in responsibly-
  //        designed models that the check isn't worth it.  If we've hit a
  //        corner with the loss face, it'll be prefered anyway. If between
  //        top and bottom... I can think of little real harm of letting it
  //        prefer the top exit, given the already existing ambiguity implied
  //        by a model where layer interfaces cross within the cylinder bounds.
  //

  // (Corner-check code removed.)

  // ::::::
  // :: Populate the TravelRec and return it to the caller:
  // :

  rec.PathLength  = shortest;
  rec.TravelTime  = shortest / mVelTop[rt];
  rec.NewLoc      = loc + R3::XYZ(dir).ScaledBy(shortest);
  rec.NewDir      = dir;
  rec.Attenuation = HelperUniformAttenuation(rec.TravelTime * cmPhononFreq, mQ[rt]);

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


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  CLASS:  Tetra                                               ****
// ****                                                              ****
//
//   Methods Defined Here:
//
//     o   Tetra()
//     o   Face()
//     o   GetVelocAtPoint()
//     o   GetWavelengthAtPoint()
//     o   GetDensityAtPoint()
//     o   GetQatPoint()
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
// METHOD:   Tetra :: Face()
//
CellFace & Tetra::Face(Index face_id) {
  return mFaces[face_id];
}


//////
// METHODS:  Tetra :: GetVelocAtPoint()
//           Tetra :: GetWavelengthAtPoint()
//           Tetra :: GetDensityAtPoint()
//           Tetra :: GetQatPoint()
//
Real Tetra::GetVelocAtPoint(const R3::XYZ & loc, raytype type) const {
  return  loc.Dot(mVelGrad[type]) + mVel0[type];
}
Real Tetra::GetWavelengthAtPoint(const R3::XYZ & loc, raytype type) const {
  return  (GetVelocAtPoint(loc,type)/ cmPhononFreq);
}
Real Tetra::GetDensityAtPoint(const R3::XYZ & loc) const {
  return  loc.Dot(mDensGrad) + mDens0;
}
Real Tetra::GetQatPoint(const R3::XYZ & loc, raytype type) const {
  return mQ[type];
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
  rec.Attenuation = HelperUniformAttenuation(rec.TravelTime * cmPhononFreq, mQ[rt]);

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
// CONSTRUCTOR:  SphereShell()
//
// V(r) = A*r*r + C
// V(rT) = A*rT2 + C = vT
// V(rB) = A*rB2 + C = vB
//         A(rT2-rB2) = (vT-vB)
//
SphereShell::SphereShell(Real RadTop, Real RadBot,
                           const GridData & DataTop, const GridData & DataBot) :
  mFaces( SphereFace(RadTop, CellFace::F_TOP, this),
          SphereFace(RadBot, CellFace::F_BOTTOM, this))
{
  assert(RadTop > RadBot);
  assert(RadBot >= 0);
  assert(DataBot.Vp() >= DataTop.Vp()); // Disallow reverse gradients; Hope to relax this
  assert(DataBot.Vs() >= DataTop.Vs()); // restriction once I figure out how to handle it.

  // Coefficient 'A' on the quadratic term:  ( v = A*r^2 + C )
  mVelCoefA[RAY_P] = (DataTop.Vp()-DataBot.Vp()) / (RadTop*RadTop - RadBot*RadBot);
  mVelCoefA[RAY_S] = (DataTop.Vs()-DataBot.Vs()) / (RadTop*RadTop - RadBot*RadBot);
  mDensCoefA = (DataTop.Rho()-DataBot.Rho()) / (RadTop*RadTop - RadBot*RadBot);

  // Coefficient 'C', the constant term:
  mVelCoefC[RAY_P] = DataTop.Vp() - (mVelCoefA[RAY_P]*RadTop*RadTop);
  mVelCoefC[RAY_S] = DataTop.Vs() - (mVelCoefA[RAY_S]*RadTop*RadTop);
  mDensCoefC = DataTop.Rho() - (mDensCoefA*RadTop*RadTop);

  // Zero-surface for velocities:
  Real uradP = -mVelCoefC[RAY_P] / mVelCoefA[RAY_P];
  Real uradS = -mVelCoefC[RAY_S] / mVelCoefA[RAY_S];
  mZeroRadius2[RAY_P] = uradP;
  mZeroRadius2[RAY_S] = uradS;
  mZeroRadius[RAY_P] = (uradP > 0) ? sqrt(uradP) : (1./0.);
  mZeroRadius[RAY_S] = (uradS > 0) ? sqrt(uradS) : (1./0.);

  // Q Values: (Take top value and apply to whole cell)
  mQ[RAY_P] = DataTop.Qp();
  mQ[RAY_S] = DataTop.Qs();

  std::cout << "~~> SphereShell: Vp = " << mVelCoefA[RAY_P] << "*r^2 + " << mVelCoefC[RAY_P] << "; ZeroSurf at r = " << mZeroRadius[RAY_P] << "\n";
  std::cout << "                 Vs = " << mVelCoefA[RAY_S] << "*r^2 + " << mVelCoefC[RAY_S] << "; ZeroSurf at r = " << mZeroRadius[RAY_S] << "\n";

  assert(mZeroRadius[RAY_P] >= 0);
  assert(mZeroRadius[RAY_S] >= 0);
  assert(mZeroRadius[RAY_P] > RadTop);  // May be able to relax this one in future...
  assert(mZeroRadius[RAY_S] > RadTop);  //

}//
//

//////
// METHOD:   SphereShell :: Face()
//
//   Returns a read-write reference to either the top or bottom
//   CellFace object, as determined by the value of face_id
//
CellFace & SphereShell::Face(Index face_id) {
  return mFaces[face_id];
}

//////
// METHODS:  SphereShell :: GetVelocAtPoint()
//           SphereShell :: GetWavelengthAtPoint()
//           SphereShell :: GetDensityAtPoint()
//           SphereShell :: GetQatPoint()
//
Real SphereShell::GetVelocAtPoint(const R3::XYZ & loc, raytype type) const {
  return mVelCoefC[type] + mVelCoefA[type]*loc.MagSquared();
}
Real SphereShell::GetWavelengthAtPoint(const R3::XYZ & loc, raytype type) const {
  return  (GetVelocAtPoint(loc,type)/ cmPhononFreq);
}
Real SphereShell::GetDensityAtPoint(const R3::XYZ & loc) const {
  return mDensCoefC + mDensCoefA*loc.MagSquared();
}
Real SphereShell::GetQatPoint(const R3::XYZ & loc, raytype type) const {
  return mQ[type];
}


//////
// METHOD:   SphereShell :: GetPathToBoundary()
//
//   Compute a ray path that extends until it hits a boundary of the cell.
//
//   We actually dispatch here based on whether cell codes a Degree 0 (uniform
//   velocity) or Degree 2 (quadratic velocity) elastic profile.  The former
//   assumes straight-line rays, the latter assumes circular arc rays.
//
TravelRec SphereShell::
GetPathToBoundary(raytype rt, const R3::XYZ & loc, const S2::ThetaPhi & dir) {
  if (mVelCoefA[rt] < 0) {
    return GetPath_Variant_D2(rt, loc, dir);  // Bottoming ray arcs
  } else if (mVelCoefA[rt] == 0) {
    return GetPath_Variant_D0(rt, loc, dir);  // Straight-line rays
  } else {
    // TODO: Eventual handler for cresting ray arcs...
    throw std::runtime_error("SphereShell: No handler for inverted radial velocity profiles.");
  }
}

//////
// METHOD:   SphereShell :: GetPath_Variant_D0()
//
//   Handler for GetPathToBoundary() wherin we assume uniform velocity, and
//   corresponding straight-line ray paths.
//
TravelRec SphereShell::
GetPath_Variant_D0(raytype rt, const R3::XYZ & loc, const S2::ThetaPhi & dir) {

  TravelRec rec;
  enum exitface_e {TOP=CellFace::F_TOP,
                   BOTTOM=CellFace::F_BOTTOM,
                   NUM_FACES=2};
  Real dists[NUM_FACES];

  dists[TOP] = mFaces[TOP].LinearRayDistToExit(loc, dir);
  dists[BOTTOM] = mFaces[BOTTOM].LinearRayDistToExit(loc, dir);

  exitface_e ef = (dists[TOP] < dists[BOTTOM]) ? TOP : BOTTOM;
  if (dists[ef] < 0) dists[ef] = 0; // Squash retrograde motion

  // ::::::
  // :: Populate the TravelRec and return it to the caller:
  // :

  rec.PathLength  = dists[ef];
  rec.TravelTime  = dists[ef] / mVelCoefC[rt];
  rec.NewLoc      = loc + R3::XYZ(dir).ScaledBy(dists[ef]);
  rec.NewDir      = dir;
  rec.Attenuation = HelperUniformAttenuation(rec.TravelTime * cmPhononFreq, mQ[rt]);
  rec.pFace = &mFaces[ef];

  return rec;

}

//////
// METHOD:   SphereShell :: GetPath_Variant_D2()
//
//   Handler for GetPathToBoundary() wherin we assume quadratic radial
//   velocity, and ray paths that are circular arcs.  Note in this version we
//   take for granted a velocity profile that increases with depth, and that we
//   are below the "zero surface" (i.e. that velocities are positive).  I.e. we
//   assume that ray arcs are "bottoming" ray arcs, as opposed to "cresting"
//   ray arcs.
//
//   The procedure for computing the ray path begins with finding the center
//   point of the ray arc, which is a function of the current location and
//   direction, and of the velocity profille in the cell.
//
TravelRec SphereShell::
GetPath_Variant_D2(raytype rt, const R3::XYZ & loc, const S2::ThetaPhi & dir) {

  enum exitface_e {TOP=CellFace::F_TOP,
                   BOTTOM=CellFace::F_BOTTOM,
                   NUM_FACES=2};

  Real dists[NUM_FACES];
  RayArcAttributes RayArc = GetRayArcD2(rt, loc, dir);

  dists[TOP] = mFaces[TOP].CircularArcDistToExit(loc, dir, RayArc);
  dists[BOTTOM] = mFaces[BOTTOM].CircularArcDistToExit(loc, dir, RayArc);

  exitface_e ef = (dists[TOP] < dists[BOTTOM]) ? TOP : BOTTOM;
  if (dists[ef] < 0) dists[ef] = 0; // Squash retrograde motion

  std::cout << "PathToBoundary: Two faces considered: dists[TOP,BOT] = (" << dists[TOP] << "," << dists[BOTTOM] << ");   (Chose: Face " << ef << " at dist " << dists[ef] <<")\n";

  //FAKE:
  TravelRec rec;
  rec.PathLength  = dists[ef];
  rec.TravelTime  = 20;
  rec.NewLoc      = loc;
  rec.NewDir      = dir;
  rec.Attenuation = 1.0;
  rec.pFace = &mFaces[ef];
  return rec;

  throw std::runtime_error("UnimpSphereShell_GetPath_D2");
  return TravelRec(); // **FIXME**
}

//////
// METHOD:   SphereShell :: GetRayArcD2()
//
//  Computes Ray Arc attributes in radial quadratic velocity gradients.
//
//  Method of Solution: Bottoming radius:
//
//  Assume ray path follows { sin i = G/r * (A*r^2 + C) }, where the paren-
//  thetical is the velocity profile and G is a boundary value constant. We
//  bottom when sin i == 1, which yields: { G*A*r^2 - r + G*C = 0 }, and
//  solution:
//              r_bottom = ( 1 - sqrt(1-4*G*G*A*C) ) / 2*G*A,
//
//  where we take the negative root on the grounds that we are assuming A
//  negative and C positive.  The parameter G comes from assuming that our
//  initial equation must be satisfied for a known initial location, take-off
//  angle, and velocity, by:
//
//              G = (sin i0) * r0 / v0
//
//  Method of Solution: Arc radius:
//
//  We assume that the arc center is at the intersection of two tangent-lines
//  from the v=0 isosurface (since the ray arc must approach perpendicular as
//  it hits the isosurface.)  Thus the center is at the corner of a right-
//  triangle with side-lengths r_zero, r_arc, and hypotenuse (r_bottom +
//  r_arc), allowing solution from: (r_zero and r_bottom are known)
//
//              r_zero^2 + r_arc^2 = (r_bottom + r_arc)^2
//
//  Method of Solution: Arc center:
//
//  Once the arc radius is known, we simply move perpendicular to 'dir' (the
//  current ray tangent) a distance r_arc from the 'loc', and we find the
//  center.
//
RayArcAttributes SphereShell::
GetRayArcD2(raytype rt, const R3::XYZ & loc, const S2::ThetaPhi & dir) const {

  RayArcAttributes Ray;

  // REFERENCE BASIS anchored on CURRENT LOCATION:
  // These will be orthonormal unless dir is straight up or down, then v3=v2=0.
  R3::XYZ v3 = ECS.GetDown(loc);  // Downwards; Also direction of veloc gradient.
  R3::XYZ v2 = v3.Cross(dir).UnitElse(R3::XYZ(0,0,0));  // Out-of-plane.
  R3::XYZ v1 = v2.Cross(v3);      // Tangent to iso-surface, oriented "forward"
                                  //                 w.r.t. direction of travel.

  // Determine incidence angle-sine into isosurface:
  Real sini = v1.Dot(dir);          // Positive or zero.
  if (sini > 1.0) sini = 1.0;       // (very unlikely, numerical error only.)
  Real cosi = v3.Dot(dir);          // [-1, 1]
  // (NOTE: If 'dir' is pure vertical then we'll have sini==0, cosi==+/-1.)

  // Get BOTTOMING Radial Coordinate:
  const Real G = sini * loc.Mag() / GetVelocAtPoint(loc, rt);
  const Real TwoGA = 2. * G * mVelCoefA[rt];
  const Real urad = 1. - (2. * TwoGA * G * mVelCoefC[rt]);
        // urad >= 1, assuming downward-gradient velocity;
        // urad == 1 implies singular straight up/down 'dir' direction.
        // urad  < 1 implies assumptions violated.
  Ray.Bottom = (urad>1) ? (1. - sqrt(urad)) / TwoGA : 0;

  // ARC RADIUS:
  Ray.Radius = (mZeroRadius2[rt]/Ray.Bottom - Ray.Bottom) / 2.0;

  // (NOTE: If 'dir' is vertical then we'll have Ray.Bottom == 0.)
  // (NOTE: If 'dir' is vertical then we'll have Ray.Radius == +inf.)
  // (NOTE: We check 'urad' again below for another detection of vertical up/down.)

  // VECTOR from 'loc' to ARC CENTER:
  R3::XYZ LocToArcCenter = v1.ScaledBy(Ray.Radius*cosi) + v3.ScaledBy(-Ray.Radius*sini);
  Ray.Center = loc + LocToArcCenter;

  // REFERENCE BASIS anchored on ARC CENTER:
  Ray.u3 = ECS.GetDown(Ray.Center); // Points towards center-of-Earth
  Ray.u2 = v2;                      // Out-of-plane
  Ray.u1 = Ray.u2.Cross(Ray.u3);    // Forward-tangent to arc at bottom.

  // Catch straight up/down:
  if (urad <= 1) {
    Ray.Center = R3::XYZ(0,0,0);
    Ray.u3 = Ray.u2 = R3::XYZ(0,0,0);
    Ray.u1 =  dir;
  }

  /*
  std::cout << "GRAD2: OrNormal check: [" << v3.MagSquared() << ", " << v1.MagSquared() << ", " << v1.Cross(v2).MagSquared() << ", " << v2.Cross(v3).MagSquared() << "]"
            << ", Sin(i) is: " << sini << " Cos(i) is: " << cosi << ";  Uptrending-at-loc: " << -(v3.Dot(dir)) << "\n";
  std::cout << "  BottomRad: " << Ray.Bottom << " ArcRadius: " << Ray.Radius << " Sum: " << (Ray.Bottom+Ray.Radius) << " (ArcCenterMag: " << Ray.Center.Mag() << ") Dir.Theta: " << dir.Theta()
            << "  {{ G: " << G << " urad: " << urad << " 1-sqrad: " << (1.-sqrt(urad)) << " }}\n";
  std::cout << "  ArcNormal check: " << Ray.u1.Cross(Ray.u2).Mag() << ", " << Ray.u2.Cross(Ray.u3).Mag() << ", " << Ray.u1.Mag() << "  u1 Parallality check: " << Ray.u1.Dot(dir) << "  Uptrending check: " << -(Ray.u3.Dot(dir)) << "  Out-of-plane check: " << Ray.u2.Dot(dir) << "\n";
  std::cout << "  Location now is:  " << loc.str() << "  direction: " << R3::XYZ(dir).str() << "  v1: " << v1.str() << " v3: " << v3.str() << "\n";
  std::cout << "  RayAttributes is: " << Ray.str() << "\n";
  */

  return Ray;

}

//////
// METHOD:   SphereShell :: AdvanceLength()
//
TravelRec
SphereShell::AdvanceLength(raytype rt, Real len, const R3::XYZ & startloc,
                           const S2::ThetaPhi & startdir) {
  if (mVelCoefA[rt] == 0) {
    return AdvanceLength_Variant_D0(rt, len, startloc, startdir);
  } else {
    return AdvanceLength_Variant_D2(rt, len, startloc, startdir);
  }
}


//////
// METHOD:   SphereShell :: AdvanceLength_Variant_D0()
//
TravelRec SphereShell::
AdvanceLength_Variant_D0(raytype rt, Real len, const R3::XYZ & startloc, const S2::ThetaPhi & startdir) {

  TravelRec rec;

  rec.PathLength = len;
  rec.TravelTime = len / mVelCoefC[rt];
  rec.NewLoc     = startloc + R3::XYZ(startdir).ScaledBy(len);
  rec.NewDir     = startdir;
  rec.Attenuation = HelperUniformAttenuation(rec.TravelTime * cmPhononFreq, mQ[rt]);

  return rec;

}


//////
// METHOD:   SphereShell :: AdvanceLength_Variant_D2()
//
TravelRec SphereShell::
AdvanceLength_Variant_D2(raytype rt, Real len, const R3::XYZ & startloc, const S2::ThetaPhi & startdir) {

  //FAKE:
  TravelRec rec;
  rec.PathLength  = 1e10;
  rec.TravelTime  = 1e10;
  rec.NewLoc      = startloc;
  rec.NewDir      = startdir;
  rec.Attenuation = 1.0;
  rec.pFace = &mFaces[1];
  return rec;

  throw std::runtime_error("UnimpSphereShell_AdvanceLength_V_D2");
  return TravelRec(); // **FIXME**
}
