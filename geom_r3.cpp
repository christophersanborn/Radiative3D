// geom_r3.cpp
//
#include "geom_r3.hpp"
#include "geom_s2.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cerrno>     /* Error checking in math functions */
using namespace std;

// IN THIS FILE:
// Class Implementations for:
//
//   o  Classes in the R3:: namespace
//
//   Search on "&&&&" to jump to the next class implementation
//   block.
//
//

//////////////////////////////////////////////////////////////////////////
//************************************************************************
//* %%%%                                                             #%%#
//* #%%#    IMPLEMENTATIONS for classes in the R3:: namespace        #%%#
//* #%%#                                                             #%%#
//*         Including:
//*
//*         o  class XYZ
//*         o  class Matrix
//*         o  class OrthoAxes
//*
//************************************************************************


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  CLASS:  XYZ                                                 ****
// ****                                                              ****
//

R3::XYZ::XYZ(const S2::ThetaPhi & thph) :
  mX  ( sin(thph.Theta()) * cos(thph.Phi()) ),
  mY  ( sin(thph.Theta()) * sin(thph.Phi()) ),
  mZ  ( cos(thph.Theta())                   )  {
}


//////
// METHOD:  XYZ :: str()
//
std::string R3::XYZ::str(int fwidth) const {
  std::ostringstream s;
  s << std::setprecision(4);
  s << "(" << mX << "," << mY << "," << mZ << ")";
  return s.str();
}

//////
// METHOD:  XYZ :: Outer()
//
//   Compute and return the outer product between two vectors.  This
//   is the tensor produced by multiplying V1 * V2^T.
//
//   Take (*this) to be V1 and 'other' to be V2, such that the outer
//   product OP can be generated via a call like this:
//
//   OP = V1.Outer(V2);
//
const R3::Matrix R3::XYZ::Outer(const R3::XYZ & other) const {

  return R3::Matrix(mX*other.mX, mX*other.mY, mX*other.mZ,
                    mY*other.mX, mY*other.mY, mY*other.mZ,
                    mZ*other.mX, mZ*other.mY, mZ*other.mZ);

}


//////
// METHOD:  XYZ :: ThetaHat()                        [R3:: Namespace]
//
//   Returns a unit-vector in the theta^hat (direction of increasing
//   co-latitude) direction relative to the direction of the ThetaPhi
//   object.
//
const R3::XYZ R3::XYZ::ThetaHat() const {

  Real theta = Theta();     // Current theta and phi components
  Real phi = Phi();         //

  Real rth;                 // Theta and Phi components of 
  Real rph;                 // return value

  if (theta < Geometry::Pi90) {    // Simple case

    rth = Geometry::Pi90 + theta;
    rph = phi;

  }                                 // Else need to flip phi in order to
  else {                            // keep our values in the favored
                                    // domain
    rth = Geometry::Pi270 - theta;
    rph = (phi < Geometry::Pi180) ? phi + Geometry::Pi180
      : phi - Geometry::Pi180;

  }

  return R3::XYZ( sin(rth) * cos(rph),
                  sin(rth) * sin(rph),
                  cos(rth)             );

}


//////
// METHOD:  XYZ :: PhiHat()                          [R3:: Namespace]
//
//   Returns a unit-vector in the phi^hat (direction of increasing
//   azimuth) direction relative to the direction of the XYZ
//   object.
//
const R3::XYZ R3::XYZ::PhiHat() const {

  Real rph = Phi() + Geometry::Pi90;
  return R3::XYZ( cos(rph), sin(rph), 0 );

}


//////
// METHOD:  XYZ :: GetInPlaneUnitPerpendicular()
//
//   Returns a vector perpendicular to (*this) that is coplanar with
//   (*this) and 'other'.  (There are two possible solutions.  This
//   function returns the "positive" solution - the one that satisfies
//   solution.Dot(other) > 0.)
//
//   This function as currently written is kindof a hack, and could
//   benefit from a more thorough analysis of the numerical
//   instabilities, and could also probably be written in more
//   efficient way.
//
//   TODO: Make this function better, clearer, faster.  Or at a
//   minimum, at least prove to myself that it is numerically stable.
//   Fingers crossed for now.
//
const R3::XYZ
R3::XYZ::GetInPlaneUnitPerpendicular(const R3::XYZ & other) const {
  // Get mutually perpendicular vector:
  R3::XYZ mutualperp;
  mutualperp = (*this).Cross(other);
  if (mutualperp.IsZero()) {    // Happens if 'other' is || *this 
    mutualperp = (*this).Cross(R3::XYZ(1,0,0));   // Try another other
    if (mutualperp.IsZero()) {  // Happens if 'other' is || (1,0,0)
      mutualperp = (*this).Cross(R3::XYZ(0,1,0)); // and another
    }
  } // If it wasn't before, it's perpendicular now.
  mutualperp.Normalize();
  // At this point, mutualperp is normalized, and *probably*
  // perpendicular to (*this), but if 'other' and (*this) were not
  // well separated, then we may have just normalized a numerically
  // unstable quantity, which could make us less than fully
  // orthogonal.  However, at this point mutualperp and (*this) ARE
  // well-separated (I think), and thus one more cross-product should
  // give us a good perpendicular.  Furthermore, if 'other' and
  // (*this) were reasonably separated, the result of this
  // cross-product will be coplanar with 'other' and (*this).
  R3::XYZ result;
  result = mutualperp.Cross(*this);
  result.Normalize();
  return result;
}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  CLASS:  Matrix                                              ****
// ****                                                              ****
//

//////
// METHOD:  Matrix :: OutputContents()
//
//   Print matrix contents, with decorations, to stdout.
//
void R3::Matrix::OutputContents() const {
  Real mag = Mag();
  Real tr  = Trace();
  Real iso = tr*tr / (3.0*mag*mag);
  if (tr<0) {iso*=-1;}
  std::cout << "The contents of this matrix are:" << std::endl 
            << "\tx\ty\tz\t" << std::endl
            << "    +---------------------------" << std::endl
            << "  x |\t" << mxx << "\t" << mxy << "\t" << mxz 
            << "\t\tMagnitude: " << mag << "\n"
            << "  y |\t" << myx << "\t" << myy << "\t" << myz
            << "\t\tTrace:     " << tr << "\n"
            << "  z |\t" << mzx << "\t" << mzy << "\t" << mzz
            << "\t\tIso Frac:  " << iso << "\n"
            << std::endl;
}


//////////////////////////////////////////////////////////////////////////
// &&&&                                                              ****
// ****  CLASS:  OrthoAxes                                           ****
// ****                                                              ****
//

//////
// CONSTRUCTOR:   OrthoAxes()
//
R3::OrthoAxes::OrthoAxes(Real the, Real phi, Real rot) :
  mTheta(the), mPhi(phi), mRot(rot) {
  // Construct axes sets based on Theta, Phi, and Rotation.
  Real costhe = cos(the);
  Real sinthe = sin(the);
  Real cosphi = cos(phi);
  Real sinphi = sin(phi);
  Real cosrot = cos(rot);
  Real sinrot = sin(rot);
  //
  mE3.SetXYZ(  sinthe*cosphi,  sinthe*sinphi,  costhe );
  mE1.SetXYZ(  costhe*cosphi,  costhe*sinphi, -sinthe );
  mE2.SetXYZ( -sinphi,         cosphi,            0   );
  //
  mS1.SetXYZ(  cosrot*costhe*cosphi - sinrot*sinphi,
               cosrot*costhe*sinphi + sinrot*cosphi,
              -cosrot*sinthe                          );
  //
  mS2.SetXYZ( -sinrot*costhe*cosphi - cosrot*sinphi,
              -sinrot*costhe*sinphi + cosrot*cosphi,
               sinrot*sinthe                          );
}


//////
// METHOD:   OrthoAxes :: Express()
//
//   Does a basis-transformation of the OrthoAxes object.
//
const R3::OrthoAxes 
      R3::OrthoAxes::Express(const R3::OrthoAxes & axs) 
      const {
  // Assume that 'axs' is an OrthoAxes object which considers THIS
  // OrthoAxes object to be its basis. Returns an OrthoAxes object
  // which expresses 'axs' in the parent basis of THIS OrthoAxes
  // object.  (I.e, does a basis transform on 'axs').

  using namespace R3;
  OrthoAxes ret;                  // This will hold the return value

  ret.mS1 = Express(axs.mS1);     // Directly transform the S vectors
  ret.mS2 = Express(axs.mS2);     //
  ret.mE3 = Express(axs.mE3);     // (E3 is identical with S3)

                                  // We won't need 'axs' anymore. All
                                  // that we need is now in 'ret'.
  
  errno = 0;                      // Clear errno (#include <cerrno>)

  Real costhe = ret.mE3.z();      // Get angles, sines, and cosines
  Real the = acos(costhe);        // (errno -> EDOM if arg not [-1,1])
  Real phi = atan2(ret.mE3.y(),   // ('' '' ''   if y,x both zero)
                     ret.mE3.x());
  Real sinthe = sin(the);         //
  Real cosphi = cos(phi);         //
  Real sinphi = sin(phi);         //

  ret.mTheta = the;               // Store theta
  ret.mPhi = phi;                 // Store phi

  ret.mE1 = XYZ(  costhe*cosphi,  // Theta^Hat direction
                  costhe*sinphi,
                 -sinthe      );

  ret.mE2 = XYZ( -sinphi,         // Phi^Hat direction
                  cosphi,
                  0    );

  Real rot_x = ret.mS1.Dot(ret.mE1);
  Real rot_y = ret.mS1.Dot(ret.mE2);

  ret.mRot = atan2(rot_y,rot_x);  // Angle S1 makes in the E1,E2 plane
                                  // (the final euler angle defining
                                  // the orientation of the coordinate
                                  // system)

  if (errno != 0) {
    std::cerr << "Domain Error in OrthoAxes basis transformation:\n"
              << " costhe = " << costhe << "\n"
              << " S3.x, S3.y = " << ret.mE3.x() << ", "
              << ret.mE3.y() << "\n"
              << " rot_x, rot_y = " << rot_x << ", " << rot_y << "\n"
              << " errno = " << errno << "\n";
    errno = 0;
    // Do not abort; Although perhaps we should.
  }

  // NOTES on NUMERICAL INSTABILITIES: I haven't put any rigorous
  // analysis into the kinds of instabilities that might occur in this
  // subroutine, but the obvious possibilities fall into two
  // categories: (1.) The possibility of orthogonality not being
  // preserved by the the multiplicative transformations of the S
  // vectors, due to round-off errors, and (2.) The possibility of
  // normality not being preserved.
  //
  // I don't do any checks for problem (1).  As for problem (2), the
  // most dangerous affect of this that I see is the possibility that
  // the determination of Theta from an Arc Cosine could suffer a
  // domain error.  At the moment, I don't do any direct domain
  // checking.  Instead I just check 'errno' at the end of everything,
  // which will let me know if SOMETHING went wrong in either the acos
  // or one of the atan2's.  These are reported on stderr, to give me
  // the chance to decide later if something further should be done.
  // It might be prudent, however, based on the notion that if
  // 'costhe' is out-of-domain it will likely be by only the slightest
  // amount, to simply do a quick check where we set 'costhe' to
  // either 1 ot -1 if it is above or below bounds, and simply
  // proceed.  This would likely be a very pragmatic solution,
  // provided OrthoAxes objects aren't subjected to lots of repeated
  // transformations that could have cumulative detrimental effects.
  //
  // So, should I worry about it?  Probably not, and here's why: At
  // least for the moment, the only use case of this OrthoAxes class
  // is in the Phonon class where it is used to rotate phonon
  // trajectories based on scattering angles.  The new trajectories DO
  // NOT get stored as OrthoAxes objects.  They get stored as Theta,
  // Phi, and Rot angles.  This means that there is only ONE basis
  // transformation to worry about before the direction gets re-coded
  // in a way that is definitionally ortho-normal.  On the next
  // scattering event, the conversion into an OrthoAxes representation
  // will begin anew from unambiguous starting parameters.  There are
  // no cumulative iterations, and thus we never have a chance to
  // wander too far from orthonormality.  So I'm not worried.  At
  // least, so far as I'm not forgetting anything...
  //

  return ret;

}
