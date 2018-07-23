// tensors.hpp
//
// This file develops representations of moment tensors for describing
// seismic source events.
//
// Note: this file used to define basic functionality of tensors
// generally, but that functionality has been removed to the R3
// namespace in geom_r3.hpp.  Tensor::Tensor is now a typedef alias to
// R3::Matrix, and this file develops derived classes whose purpose is
// primarily to provide specific interfaces to constructing
// seismic/geophysical moment tensors, while the fundamentals of
// tensor algebra are developed in the Matrix class.
//
#ifndef TENSORS_H_
#define TENSORS_H_
//
#include <stdexcept>
#include <cmath>
#include "geom_r3.hpp"

// For the moment, there is no corresponding tensors.cpp.  Everything
// here so far is defined inline.  (This may change with growing
// complexity, of course.)
//

//////
// *** FORWARD DECLARATIONS:
//


//__________________________________________________________________________
//**************************************************************************
// NAMESPACE: Tensor
// PURPOSE:
//
//   To provide classes that represent Rank-2 Tensors in three dimensions.
//
//__________________________________________________________________________
//**************************************************************************

namespace Tensor {

//////
// TYPEDEFS:
//
typedef R3::Matrix Tensor;  // Internally, our Tensor "base class" is
                            // just a 3x3 matrix.  In the future, we
                            // may wish for some extended funtionality
                            // that may motivate us to derive from
                            // R3::Matrix and expand it, but for now a
                            // typedef alias does the trick.


//////
// CLASSES: Forward Declarations
//

  class Symmetric;    // Symmetric: Specify 6 elements
  class Sym_Eigs;     // Symmetric: Specify 3 eigenvalues and three
                      // euler rotation angles
  class USGS;         // Like Symmetric, but uses Harvard/USGS sign
                      // convention (as opposed to Aki and Richards
                      // which we use in the others.)
  class Euler;        // Rotation
  class SDR;          // Creates Moment Tensor using strike, dip, and rake.

//////
// CLASSES:
//

//////
// Notes from old Tensor base class (before we removed the
// functionality to R3::Matrix and aliased via a typedef):
//
//   This is a base class and generic representation of a Rank-2
//   Dimension-3 tensor (more-or-less a 3x3 matrix).  It provides much
//   of the basic interface to tensor objects or other 3x3 things like
//   rotation matrices, but it provides minimal options for
//   constructor overloads.  For variations on how to "specify" the
//   tensor (e.g. as Strike Dip and Rake, Euler angles, USGS
//   convention, etc.), look at the derived classes, the purpose of
//   which is basically to provide names to the different ways of
//   specifying the elements.  (I.e., derived classes should provide
//   custom constructors, but not add data members or addtional
//   functionality. (Barring a reason to do so not anticipated here.))
//
//   The tensor elements are indexed in the standard row,column named
//   index manner using x,y,z as the index symbols. These are accessed
//   by name via calls to access function xx(), xy(), ... zz().
//
//   Numeric indicial access to the elements may be provided in the
//   future (but hasn't been needed yet) whereby access would be
//   something like element(i,j) with i,j being integers indexing the
//   rox,column element, with 0->x, 1->y, 2->z.
//
//   Note that the tensor class does not care about the choice of
//   basis that x,y,z are intended to represent, other than assumeing
//   a standard right-handed orthogonal coordinate system,... but
//   "which way is up" is not something that the class knows about or
//   cares about.  However, for REFERENCE (provided for the benefit of
//   someone looking here rather than elsewhere), note that the USER
//   INTERFACE code in OTHER modules, when interpreting MOMENT TENSORS
//   provided by the USER, will interpret using the NED convention, or
//   x,y,z <--> North,East,Down, which is favored by Aki and Richards.
//   This convention might not be the one used for internal modelling,
//   and the code that initializes the model and event sources is
//   tasked with making any necessary ROTATIONS to the provided moment
//   tensors so that they are oriented correctly in the actual Earth
//   model.  But do note... I just put this note here for convenience.
//   For definitive documentation of how moment tensors are
//   interpreted, look in the appropriate other modules.
//
//  


//////
// CLASS:  Tensor::Symmetric
//
//   A symmetric 3x3 tensor object, specified by the diagonal elements
//   and then the three off-diagonals in the upper triangle.
//
class Symmetric : public Tensor {
public:
  Symmetric(Real xx, Real yy, Real zz,
            Real xy = 0, Real xz = 0, Real yz = 0) :
    Tensor(xx, xy, xz, xy, yy, yz, xz, yz, zz)  {
  }
  ~Symmetric() {}
};


//////
// CLASS:  Tensor::USGS
//
//   Similar to Tensor::Symmetric, but parameters specified following
//   the USE axes convention of Harvard/USGS earthquake moment tensor
//   notation (C.f. the NED convention by which we interpret the
//   general moment tensor).  The element ordering in the constructor
//   args is is the order that the elements are listed in the CMT
//   database catalog.
//
//   Here the r,t,p axes notation corresponds to r, theta, phi in an
//   Earth-centered polar coordinate system, such that r->Up,
//   t->South, and p->East as measured on the Earth's surface (e.g. at
//   the epicenter).
//
class USGS : public Tensor {
public:
  USGS(Real rr, Real tt, Real pp,
       Real rt = 0, Real rp = 0, Real tp = 0) :
    Tensor( tt, -tp,  rt, 
           -tp,  pp, -rp, 
            rt, -rp,  rr) {}
  ~USGS() {}
};


//////
// CLASS:  Tensor::EulerSDR
//
//   Produces a rotation matrix from a set of Euler angles. This
//   matrix can then be used to rotate, say, a moment tensor (e.g. to
//   achieve desired Strike, Dip, and Rake angles) by multiplying like
//   so: M = E*D*E.T(), where M is the desired moment tensor, D is
//   some default moment tensor (e.g., with S,D,R = 0,0,0), E is the
//   Euler matrix, and E.T() is its transpose.
//
//   Note that there are multiple ways to define a set of Euler angles
//   and that they do not all produce the same result.  The convention
//   used here follows that of the software package MoPaD, which is a
//   tool for plotting moment tensor beach balls and analyzing moment
//   tensor properties. (http://www.mopad.org/) The MoPaD Euler angles
//   have a direct mapping between Strike, Dip, and Rake, in which:
//
//   alpha is the dip angle,
//   beta is the strike angle, and
//   gamma is the negative rake angle.
//
//   All angles are assumed to have been converted from degrees to
//   radians if necessary.
//
class EulerSDR : public Tensor {
public:
  EulerSDR(Real alpha, Real beta, Real gamma) {

    Real ca = cos(alpha);
    Real cb = cos(beta);
    Real cg = cos(gamma);
    Real sa = sin(alpha);
    Real sb = sin(beta);
    Real sg = sin(gamma);

    mxx = cb*cg-ca*sb*sg;       // (Note that this is the transpose of
    mxy = -cb*sg-ca*sb*cg;      // the way these are defined in MoPaD
    mxz = sa*sb;                // because I prefer to left multiply
    myx = sb*cg+ca*cb*sg;       // by the UN-transposed Euler matrix,
    myy = -sb*sg+ca*cb*cg;      // e.g., M=E*D*E.T(), as opposed to
    myz = -sa*cb;               // M=E.T()*D*E which is how MoPaD does
    mzx = sa*sg;                // it.)
    mzy = sa*cg;
    mzz = ca;

  };
};

//////
// CLASS:  Tensor::SDR
//
//   Takes a set of Strike, Dip, and Rake angles and produces a
//   double-couple moment tensor of the desired orientation.  Angles
//   are specified according to the normal conventions (I think) and
//   should be specified in DEGREES.
//
//   Strike in expected in range [-180, 360] and represents degrees
//   east of north.  Dip is expected in range [0, 90].  Rake is
//   expected in range [-180, 180].  But I'm pretty sure nothing
//   breaks if numbers are outside that range (it just becomes hard to
//   visualize the resultant orientation) so no bounds checking is
//   done.  Might code some non-terminating warnings at some point if
//   it seems worthwhile.
//
//   Optionally takes a moment value to set the magnitude of the
//   resultant moment tensor.  The magnitude of the moment tensor has
//   no bearing on the operation of Radiative3D (at present, anyway),
//   but may be of interest to the user as it can scale the displayed
//   values of the tensor when printed on the output stream. A
//   negative value negates the moment tenor.
//
//   Optionally takes an "isotropic fraction", (a value in the range
//   [-1.0, 1.0]), that determines what fraction of the squared
//   magnitude of the resultant MT comes from the isotropic vs. the
//   deviatoric part.  Default is zero, which means pure deviatoric.
//   Negative values imply implosive isotropic part.
//
//   Note: Another way to specify the relative contribution of
//   isotropic and deviatoric moment is via an isotropic angle, as
//   defined (rather intuitively, in my opinion) in Bukchin et al 2001
//   - Isotropic and Nonisotropic Components.  A static public method
//   is provided to convert from iso angle to iso fraction, and can be
//   used by the calling code if it prefers to specify in this way.
//
class SDR : public Tensor {
public:
  SDR(Real strike, Real dip, Real rake, Real iso=0.0, Real moment=1.0) :
    Tensor(0,  0, -1,   // Default (SDR=0,0,0) tensor
           0,  0,  0,   //
          -1,  0,  0)   //
  {

    // Bounds Checking:
    if (iso > 1.0 || iso < -1.0) throw(std::domain_error(
          "Tensor::SDR: iso not in [-1.0, 1.0]"));
    if (moment == 0) throw(std::domain_error(
          "Tensor::SDR: moment must be non-zero"));

    // Rotate basic Double-Couple to desired SDR:
    strike *= Geometry::DtoR;     // Degrees to radians
    dip *= Geometry::DtoR;        //
    rake *= Geometry::DtoR;       //
    EulerSDR rot(dip,strike,-rake);  // Rotation matrix
    Transform(rot);               // Apply to self

    Real I2 = std::fabs(iso);     // Isotrpic and Deviatoric 
    Real D2 = 1.0-I2;             // squared magnitudes, contrived to
                                  // put resultant matrix on unit
                                  // sphere.

    SetSquaredMag(D2);            // (*this) contains the deviatoric
                                  // part. Set squared mag to D2.

    iso = (iso >= 0) ? 1 : -1;    // Pick positive or negative isotropic.
    Tensor ISO = Symmetric(iso, iso, iso);  // Diagonal matrix; using
    ISO.SetSquaredMag(I2);                  // 'iso' on the diagonal
                                            // gets the sign right.
    
    (*this) += ISO;               // Combine iso and devi parts.
    (*this) *= moment;            // Sets desired total moment.
                                  // (Scales unit matrix by scalar
                                  // moment.)
  }


public:
  ;
  // :::::::::::::::::::::::::::::::::::::::::
  // ::: Auxiliary UI Methods  (SDR Class) :::
  // :::::::::::::::::::::::::::::::::::::::::

 
  //////
  // METHOD:  IsoFracFromIsoAngle()     [static]
  //
  // Conversion between representing isotropic component as an angle and
  // representing it as an energy fraction.  As an angle, we define phi
  // as:
  //
  //   tan phi = I / D
  //
  // where I is the isotropic moment and D is the deviatoric moment.  As
  // a fraction, we take the total moment M to be given by:
  //
  //   M^2 = I^2 + D^2
  //
  // and the isotropic fraction is then:
  //
  //   isofrac = I^2 / M^2
  //
  // (This definitiion has the feature that the isofrac and
  // correspondingly defined devifrac sum to unity.)  From these
  // definitions, the conversion from phi to isofrac is:
  //
  //   isofrac = (tan phi)^2 / (1 + (tan phi)^2)
  //
  // With a little bit of care, we can preserve the sign on phi so that
  // explosions (isofrac > 0) and implosions (isofrac < 0) can be
  // distinguised.
  //
  // INPUT:
  //
  // Since the strike, dip, and rake angles of teh SDR class are all
  // teken in in degrees, we take the isoangle here to also be in
  // degrees. Valid values range from -90, representing pure implosion,
  // to +90, representing pure explosion, with 0 representing pure
  // deviatoric (no isotropic component).
  //
  static Real IsoFracFromIsoAngle(Real isoangle) {

    // Bounds Checking:
    if (isoangle > 90.0 || isoangle < -90.0) throw(std::domain_error(
          "Tensor::SDR: isoangle not in [-90.0, 90.0]")); 
    if (isoangle <= -89.999) return(-1.0);  // Avoid infinities from
    if (isoangle >=  89.999) return(1.0);   // tangent function

    Real tiso = tan(Geometry::DtoR * isoangle);
    Real tiso2 = tiso * tiso;

    Real isofrac = tiso2 / (1 + tiso2);

    if (isoangle < 0) {isofrac *= -1;}      // Negative angles for
                                            // implosive events
    return(isofrac);

  }//
  //

};// End class SDR
///
  
};// End namespace Tensor
///
//

#endif //#ifndef TENSORS_H_
