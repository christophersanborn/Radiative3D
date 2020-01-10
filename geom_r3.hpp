// geom_r3.hpp
//
// The geom library defines namespaces containing classes and other
// constructs related to the geometry of spaces in various dimensions.
// This file defines:
//
// o Namespace: R3
//
//    For point and other geometric primatives in a flat three-space.
//    Besically, defines a three-component point vector class.
//
//
#ifndef GEOM_R3_H_
#define GEOM_R3_H_
//
#include <stdexcept>
#include <cmath>
#include "geom_base.hpp"


//////
// *** FORWARD DECLARATIONS:
//
namespace S2{       /* Definitions in geom_s2.hpp */

  class ThetaPhi;

};
namespace R3{       /* Definitions below */

  class Matrix;

};


//////
// *** CLASSES:
//

//__________________________________________________________________________
//**************************************************************************
// NAMESPACE: R3
// PURPOSE:
//
//   To provide a set of classes, types, and constants to represent
//   geometry on a flat (cartesian) three-space.
//__________________________________________________________________________
//**************************************************************************

namespace R3 {


//////
// CLASS:   ::::  XYZ  ::::
//
//   Encapsulates Cartesian coordinate tripples in R3 space, and
//   provides various geometric computation methods involving vectors
//   in the space.
//
class XYZ {
protected:

  // ::::::::::::::::::::::::::::::::
  // ::: Member Data  (XYZ Class) :::
  // ::::::::::::::::::::::::::::::::

  Real mX;
  Real mY;
  Real mZ;


public:

  // :::::::::::::::::::::::::::::::::
  // ::: Constructors  (XYZ Class) :::
  // :::::::::::::::::::::::::::::::::

  XYZ() : mX(0), mY(0), mZ(0) {}        // Default to zero-vector

  XYZ(Real x, Real y, Real z) :         // Construct from coord tripple
    mX(x), mY(y), mZ(z) {}              //

  XYZ(const S2::ThetaPhi & thph);       // Construct a unit-vector from
                                        // a ThetaPhi object


public:

  // :::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Set Methods  (XYZ Class) :::
  // :::::::::::::::::::::::::::::::::::::::::

  void SetXYZ(Real x, Real y, Real z) {
    mX = x;
    mY = y;
    mZ = z;
  }


  // :::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Get Methods  (XYZ Class) :::
  // :::::::::::::::::::::::::::::::::::::::::
  //
  //                              These methods return information on
  //                              the internal state of the object.
  //

  Real x() const {return mX;}
  Real y() const {return mY;}
  Real z() const {return mZ;}

  Real Theta() const {
    return IsSquaredZero() 
             ? (Real)0
             : std::acos( mZ / Mag() );
  }

  Real Phi() const {
    return std::atan2(mY, mX);          // (std:: guarantees correct over-
  }                                     // loaded version for type Real.)

  Real MagSquared() const {
    return (mX*mX + mY*mY + mZ*mZ);
  }

  Real Mag() const {
    return std::sqrt(MagSquared());     // (std:: guarantees correct over-
  }                                     // loaded version for type Real.)

  bool IsZero() const {                 // Tests for "identically zero."
    return ((mX==(Real)0.0)             // Should not be used to predict
            && (mY==(Real)0.0)          // whether Mag() is zero, as small-
            && (mZ==(Real)0.0));        // enough vectors will still truncate
  }                                     // to zero on squaring.


  bool IsSquaredZero() const {          // If true, Mag() will return zero,
    return (MagSquared() == (Real)0);   // even if vector is not identically
  }                                     // zero (tiny^2 --> zero). Conversely,
                                        // false should be a pretty good guar-
                                        // antee that Mag() will be non-zero.

  const XYZ Unit() const {              // TODO: Handle sub-IsSquaredZero() 
    Real maginv = ((Real)1.0)/Mag();    // cases.
    return XYZ(mX*maginv,mY*maginv,mZ*maginv);
  }

  const XYZ UnitElse(const XYZ & fallback) const {
    Real mag = Mag();                   // Return unit vector unless we are a
    if (mag == (Real)0) return fallback;// zero-length vector, in which case
    Real maginv = ((Real)1.0)/mag;      // return fallback.
    return XYZ(mX*maginv,mY*maginv,mZ*maginv);
  }

  const XYZ Negative() const {
    return XYZ(-mX, -mY, -mZ);
  }

  const XYZ ThetaHat() const;
            // Returns a unit-vector in the theta^hat (direction of
            // increasing co-latitude) direction relative to the
            // direction of the XYZ object.

  const XYZ PhiHat() const;
            // Returns a unit-vector in the phi^hat (direction of
            // increasing co-latitude) direction relative to the
            // direction of the XYZ object.


  // ::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Object-Manipulation Methods  (XYZ Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::
  //
  //                              These methods modify the internal
  //                              state of the object.
  //

  void Normalize() {
    Real norm = 1.0 / Mag();
    mX *= norm;
    mY *= norm;
    mZ *= norm;
  }


  // ::::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Geometric Computation Methods  (XYZ Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::::::
  // 
  //                              These methods are all CONSTANT (they
  //                              do NOT modifiy object internal
  //                              state) and return whole-object
  //                              temporaries.
  //
  //                              They implement various geometric
  //                              operations involving three-vectors.
  //

  Real Dot(const XYZ & other) const {
    return (other.mX*mX + other.mY*mY + other.mZ*mZ);
  }

  const XYZ VectorTo(const XYZ & other) const {
    return XYZ(other.mX - mX,
               other.mY - mY,
               other.mZ - mZ);
  }

  Real DistFrom(const XYZ & other) const {
    return XYZ(other.mX - mX,
               other.mY - mY,
               other.mZ - mZ).Mag();
  }

  const XYZ Cross(const XYZ & other) const {
    return XYZ(mY * other.mZ - mZ * other.mY,
               mZ * other.mX - mX * other.mZ,
               mX * other.mY - mY * other.mX);
  }

  const Matrix Outer(const XYZ & other) const;

  const XYZ ScaledBy(const Real scale) const {
    return XYZ(scale * mX,
               scale * mY,
               scale * mZ);
  }

  const XYZ operator + (const XYZ & other) const {
    return XYZ(mX + other.mX,
               mY + other.mY,
               mZ + other.mZ);
  }

  const XYZ GetInPlaneUnitPerpendicular(const XYZ & other) const;
  //
  //      Returns a unit-vector perpendicular to (*this) in the plane
  //      defined by (*this) and 'other'.  If 'other' is parallel (or
  //      very nearly) to (*this) then direction of result is
  //      arbitrary, but guaranteed to be perpendicular and
  //      unit-magnitude.
  //


};// END CLASS XYZ
///


//////
// CLASS:  R3 :: Matrix
//
//   Encapsulates 3x3 matrices and the various algebraic manipulations
//   and properties thereof.
//
class Matrix {
protected:

  Real mxx, mxy, mxz;       // Matrix elements
  Real myx, myy, myz;
  Real mzx, mzy, mzz;


public:

  // ::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (Matrix Class) :::
  // ::::::::::::::::::::::::::::::::::::

  Matrix(Real xx, Real xy, Real xz,
         Real yx, Real yy, Real yz,
         Real zx, Real zy, Real zz) :
    mxx(xx), mxy(xy), mxz(xz),
    myx(yx), myy(yy), myz(yz),
    mzx(zx), mzy(zy), mzz(zz)  {
  }

  Matrix() :
    mxx(0), mxy(0), mxz(0),
    myx(0), myy(0), myz(0),
    mzx(0), mzy(0), mzz(0)  {
  }

  ~Matrix(){}


  // :::::::::::::::::::::::::::::::::::::
  // ::: Member-Access  (Matrix Class) :::
  // :::::::::::::::::::::::::::::::::::::

  Real xx() const {return mxx;}
  Real xy() const {return mxy;}
  Real xz() const {return mxz;}
  Real yx() const {return myx;}
  Real yy() const {return myy;}
  Real yz() const {return myz;}
  Real zx() const {return mzx;}
  Real zy() const {return mzy;}
  Real zz() const {return mzz;}

  XYZ Row(Index i) const {      // Return a row as an XYZ.
    switch (i) {                // Index i in set {0, 1, 2}
    case 0: return XYZ(mxx, mxy, mxz); break;
    case 1: return XYZ(myx, myy, myz); break;
    case 2: return XYZ(mzx, mzy, mzz); break;
    default: throw(std::out_of_range("Matrix::Row: Bad Index")); break;
    }
  }
  XYZ Column(Index i) const {   // Return a column as an XYZ.
    switch (i) {                // Index i in set {0, 1, 2}
    case 0: return XYZ(mxx, myx, mzx); break;
    case 1: return XYZ(mxy, myy, mzy); break;
    case 2: return XYZ(mxz, myz, mzz); break;
    default: throw(std::out_of_range("Matrix::Column: Bad Index")); break;
    }
  }

  void OutputContents() const;


  // ::::::::::::::::::::::::::::::::::
  // ::: Generators  (Matrix Class) :::
  // ::::::::::::::::::::::::::::::::::

  Matrix T() const {        // Returns the transpose of matrix
    return Matrix(mxx, myx, mzx,
                  mxy, myy, mzy,
                  mxz, myz, mzz);
  };

  Real Trace() const {
    return (mxx + myy + mzz);
  }

  Real Frobenius(const Matrix & other) const {
    return (mxx*other.mxx + mxy*other.mxy + mxz*other.mxz +
            myx*other.myx + myy*other.myy + myz*other.myz +
            mzx*other.mzx + mzy*other.mzy + mzz*other.mzz);
  }      // Returns the Frobenius double-contraction
         // (generalization of a dot product to matrices)
         // between this and the other matrix.

  Real Mag() const {
    return (sqrt(Frobenius(*this)));
  }      // Returns the magnitude of the matrix computed as the
         // Frobenius norm.

  Real Mag2() const {
    return (Frobenius(*this));
  }      // Returns the squared magnitude of the matrix computed as
         // the Frobenius self-inner-product.

  Matrix ScaledBy(Real scale) const {
    return Matrix(scale*mxx, scale*mxy, scale*mxz,
                  scale*myx, scale*myy, scale*myz,
                  scale*mzx, scale*mzy, scale*mzz);
  }

  XYZ operator*(XYZ rhs) const {  // Matrix Vector multiplication
    return XYZ( (mxx * rhs.x()) + (mxy * rhs.y()) + (mxz * rhs.z()),
                (myx * rhs.x()) + (myy * rhs.y()) + (myz * rhs.z()),
                (mzx * rhs.x()) + (mzy * rhs.y()) + (mzz * rhs.z()));
  }

  // ::::::::::::::::::::::::::::::::::::
  // ::: Manipulators  (Matrix Class) :::
  // ::::::::::::::::::::::::::::::::::::

  void ScaleBy(Real scale) {
    mxx*=scale; mxy*=scale; mxz*=scale;
    myx*=scale; myy*=scale; myz*=scale;
    mzx*=scale; mzy*=scale; mzz*=scale; 
  } // multiply by a scalar

  void Normalize(Real norm = 1.0) {
    ScaleBy(norm/Mag());
  } // Set magnitude

  void SetSquaredMag(Real n2) {
    if (n2<0) throw(std::domain_error("Matrix::SetSquaredMag: negative arg"));
    ScaleBy(sqrt(n2/Mag2()));
  } // Set squared magnitude

  void Transform(Matrix M) {
                    // Transform matrix by left-multiplying by M
                    // and right-multiplying by M.T().  Use with,
                    // e.g., a rotation matrix, to rotate a
                    // tensor/matrix.
    (*this) *= M.T();
    M *= (*this);
    (*this) = M;
  }

  Matrix & operator*=(const Matrix & rhs) {   // Matrix Multiplication
    Matrix temp;
    temp.mxx = (mxx * rhs.mxx) + (mxy * rhs.myx) + (mxz * rhs.mzx);
    temp.mxy = (mxx * rhs.mxy) + (mxy * rhs.myy) + (mxz * rhs.mzy);
    temp.mxz = (mxx * rhs.mxz) + (mxy * rhs.myz) + (mxz * rhs.mzz);
    temp.myx = (myx * rhs.mxx) + (myy * rhs.myx) + (myz * rhs.mzx);
    temp.myy = (myx * rhs.mxy) + (myy * rhs.myy) + (myz * rhs.mzy);
    temp.myz = (myx * rhs.mxz) + (myy * rhs.myz) + (myz * rhs.mzz);
    temp.mzx = (mzx * rhs.mxx) + (mzy * rhs.myx) + (mzz * rhs.mzx);
    temp.mzy = (mzx * rhs.mxy) + (mzy * rhs.myy) + (mzz * rhs.mzy);
    temp.mzz = (mzx * rhs.mxz) + (mzy * rhs.myz) + (mzz * rhs.mzz);
    (*this) = temp;
    return (*this);
  }

  Matrix & operator*=(Real rhs) {   // Scalar Multiplication
    mxx *= rhs;
    mxy *= rhs;
    mxz *= rhs;
    myx *= rhs;
    myy *= rhs;
    myz *= rhs;
    mzx *= rhs;
    mzy *= rhs;
    mzz *= rhs;
    return (*this);
  }

  Matrix & operator+=(const Matrix & rhs) {   // Matrix Addition
    mxx += rhs.mxx;
    mxy += rhs.mxy;
    mxz += rhs.mxz;
    myx += rhs.myx;
    myy += rhs.myy;
    myz += rhs.myz;
    mzx += rhs.mzx;
    mzy += rhs.mzy;
    mzz += rhs.mzz;
    return (*this);
  }

};// END CLASS Matrix
///
// Arithmetic overloads involving Matrices:
//
inline Matrix operator*(Matrix lhs, const Matrix & rhs) {
  return lhs *= rhs;    // Matrix multiplication
}
inline Matrix operator*(Matrix lhs, Real rhs) {
  return lhs *= rhs;    // Mult by scalar on right
}    
inline Matrix operator*(Real lhs, Matrix rhs) {
  return rhs *= lhs;    // Mult by scalar on left
}
inline Matrix operator+(Matrix lhs, const Matrix & rhs) {
  return lhs += rhs;    // Matrix addition
}
//
// Classes providing custom constructor interfaces to Matrix:
//
// CLASS: RowsMatrix -- Construct a Matrix object from three XYZ's
//                      specifying the rows of the matrix.
//
class RowsMatrix : public Matrix {
public:
  RowsMatrix(XYZ r1, XYZ r2, XYZ r3) :
    Matrix(r1.x(), r1.y(), r1.z(),
           r2.x(), r2.y(), r2.z(),
           r3.x(), r3.y(), r3.z())
  {}
};//
///
// CLASS: ColumnsMatrix -- Construct a Matrix object from three XYZ's
//                         specifying the columns of the matrix.
//
class ColumnsMatrix : public Matrix {
public:
  ColumnsMatrix(XYZ c1, XYZ c2, XYZ c3) :
    Matrix(c1.x(), c2.x(), c3.x(),
           c1.y(), c2.y(), c3.y(),
           c1.z(), c2.z(), c3.z())
  {}
};// END Supporting Classes to class Matrix
///


//////
// CLASS: R3::OrthoAxes
//
// ENCAPS: A right-handed set of three orthogonal axes in cartesian
//   3-space suitable to serve as the basis for a coordinate system.
//   Rotational orientation of the axes can be specified on
//   construction by the provision of three angular parameters, here
//   denoted Theta, Phi, and Rot.  Theta and Phi determine the
//   direction of the S3 axis, and Rot determines the orientation of
//   the other two.
//
//   The class actually maintains TWO sets orthogonal unit vectors,
//   which we will denote set {E} and set {S}. The first set, {E},
//   orients the E1 and E2 vectors in the Theta^Hat and Phi^Hat
//   directions, respectively, and serve as a well-defined reference
//   with which to define the orientation of the S1 and S2 axes.  We
//   interpret Rot in such a way that S1 will align along E1 when
//   Rot==0, and will align along E2 when Rot==Pi/2.
//   
class OrthoAxes {
protected:

  // ::::::::::::::::::::::::::::::::::::::
  // ::: Member Data  (OrthoAxes Class) :::
  // ::::::::::::::::::::::::::::::::::::::

  Real mTheta;
  Real mPhi;
  Real mRot;
  XYZ mE1;  // Theta^Hat Direction ("vertical" transverse direction)
  XYZ mE2;  // Phi^Hat Direction ("horizontal" transverse direction)
  XYZ mE3;  // R^Hat Direction ("outwards" or "longitudinal" direction)
  XYZ mS1;  // New vertical transverse direction after rotation
  XYZ mS2;  // New horizontal transverse direction after rotation


public:

  // :::::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (OrthoAxes Class) :::
  // :::::::::::::::::::::::::::::::::::::::

  OrthoAxes(Real the, Real phi, Real rot);
                // Construct from a set of angles

  OrthoAxes() :
    mTheta(0), mPhi(0), mRot(0),    // If no args, give standard
    mE1( XYZ(1,0,0) ),              // X,Y,Z axes set. ("Lab Frame")
    mE2( XYZ(0,1,0) ),
    mE3( XYZ(0,0,1) ),
    mS1( XYZ(1,0,0) ),
    mS2( XYZ(0,1,0) ) {/* NOP */}


  // :::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Property-Get Methods  (OrthoAxes Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::

  const Real Theta() const {return mTheta;}
  const Real Phi() const {return mPhi;}
  const Real Rot() const {return mRot;}
  const XYZ E1() const {return mE1;}  // E set is the reference set
  const XYZ E2() const {return mE2;}
  const XYZ E3() const {return mE3;}
  const XYZ S1() const {return mS1;}  // S set is fully rotated
  const XYZ S2() const {return mS2;}
  const XYZ S3() const {return mE3;}  // (outward direction unaffected
                                      // by rotation .: S3=E3.)
  

  // ::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Computation Methods  (OrthoAxes Class) :::
  // ::::::::::::::::::::::::::::::::::::::::::::::

  const XYZ Express             // Assume that 'vec' is a vector using
            (const XYZ & vec)   // this OrthoAxes as a basis. Express
            const               // this vector in the parent basis of
                                // this OrthoAxes. (I.e., do a basis
  {                             // transform on vec.)
    return XYZ( 
      vec.x()*mS1.x() + vec.y()*mS2.x() + vec.z()*mE3.x(),  // X-component
      vec.x()*mS1.y() + vec.y()*mS2.y() + vec.z()*mE3.y(),  // Y-component
      vec.x()*mS1.z() + vec.y()*mS2.z() + vec.z()*mE3.z()); // Z-component
  }

  const OrthoAxes Express(const OrthoAxes & axs) const;
                                // Does a basis transform
                                // on an OrthoAxes

};// END CLASS OrthoAxes
///

};// end namespace: R3
///
//

#endif //#ifndef GEOM_R3_H_
