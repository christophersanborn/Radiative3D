// geom_r4.hpp
//
// The geom library defines namespaces containing classes and other
// constructs related to the geometry of spaces in various dimensions.
// This file defines:
//
// o Namespace: R4
//
//    For four-element vectors and 4x4 matrices.  Includes a matrix
//    inverse method for 4x4 matrix.
//
//
#ifndef GEOM_R4_H_
#define GEOM_R4_H_
//
#include "geom_r3.hpp"

//////
// *** FORWARD DECLARATIONS:
//

//////
// *** CLASSES:
//

namespace R4 {


//////
// CLASS:    R4 :: Column
//
//   Encapsulates a column vector in R4.
//
class Column {
private:

  // :::::::::::::::::::::::::::::::::::
  // ::: Member Data  (Column Class) :::
  // :::::::::::::::::::::::::::::::::::

  Real mX[4];


public:

  // ::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (Column Class) :::
  // ::::::::::::::::::::::::::::::::::::

  Column() {        // Construct all zeros
    mX[0] = 0;
    mX[1] = 0;
    mX[2] = 0;
    mX[3] = 0;
  }

  Column(Real x1, Real x2, Real x3, Real x4) {
    mX[0] = x1;                     // Takes all four elements explicitly
    mX[1] = x2;
    mX[2] = x3;
    mX[3] = x4;
  }    

  Column(R3::XYZ xyz, Real x4) {    // Takes an XYZ for the first
    mX[0] = xyz.x();                // three elements and a Real for
    mX[1] = xyz.y();                // the last.
    mX[2] = xyz.z();
    mX[3] = x4;
  }    
    

  // ::::::::::::::::::::::::::::::
  // ::: Access  (Column Class) :::
  // ::::::::::::::::::::::::::::::

  Real x1() const {return mX[0];}       // NOTE: Named access uses
  Real x2() const {return mX[1];}       // index base 1, but...
  Real x3() const {return mX[2];}       //
  Real x4() const {return mX[3];}
  Real x(Index i) const {return mX[i];} // ...but indicial access uses
                                        // index base 0.
  R3::XYZ TruncXYZ() const
  {return R3::XYZ(mX[0],mX[1],mX[2]);}
  // Returns truncated (first three elements) vector
  // as an XYZ vector.
  //


};


//////
// CLASS:    R4 :: Matrix
//
//   Encapsulates 4x4 matrices.
//
class Matrix {

  // :::::::::::::::::::::::::::::::::::
  // ::: Member Data  (Matrix Class) :::
  // :::::::::::::::::::::::::::::::::::

  Real mX[4][4];        // First (left) index for rows, second index
                        // for columns

public:

  // ::::::::::::::::::::::::::::::::::::
  // ::: Constructors  (Matrix Class) :::
  // ::::::::::::::::::::::::::::::::::::

  Matrix(Real x11, Real x12, Real x13, Real x14,
         Real x21, Real x22, Real x23, Real x24,
         Real x31, Real x32, Real x33, Real x34,
         Real x41, Real x42, Real x43, Real x44) {
    mX[0][0] = x11;  mX[0][1] = x12;  mX[0][2] = x13;  mX[0][3] = x14;
    mX[1][0] = x21;  mX[1][1] = x22;  mX[1][2] = x23;  mX[1][3] = x24;
    mX[2][0] = x31;  mX[2][1] = x32;  mX[2][2] = x33;  mX[2][3] = x34;
    mX[3][0] = x41;  mX[3][1] = x42;  mX[3][2] = x43;  mX[3][3] = x44;
  }



  // ::::::::::::::::::::::::::::::
  // ::: Access  (Matrix Class) :::
  // ::::::::::::::::::::::::::::::
  Real x(int i, int j) const {return mX[i][j];}
  Column GetRow(int i) const {return Column(mX[i][0],mX[i][1],mX[i][2],mX[i][3]);}

  // ::::::::::::::::::::::::::::::::::::
  // ::: Result Return  (Matrix Class) ::
  // ::::::::::::::::::::::::::::::::::::

  Column SolveAXB(Column B) const;
  //   Solves Ax=b, where A is the this matrix (*this), B is
  //   provided as an argument, and x is what gets returned.
  //

};


//////
// CLASS:    R4 :: Rows  -from- R4::Matrix
//
//   Provides a custom constructor for the Matrix class allowing you
//   to specify elements "row by row", eg:
//
//       R4::Matrix M = R4::Rows(r1,r2,r3,r4);
//
//   where r1, ... r4 are R4::Column vectors, which become the rows of
//   the matrix.
//
//   Another use is to fill the rows with augmented R3 vectors, like
//   so:
//
//       R4::Matrix M = R4::Rows(v1, c1, v2, c2, v3, c3, v4, c4);
//
//   where v1, ... v4 are R3::XYZ vectors, and c1, ... c4 are scalar
//   constants, which fill out the rightmost column of the matrix.
//
class Rows : public Matrix {

public:

  Rows(Column V1, Column V2, Column V3, Column V4) :
    Matrix(V1.x1(), V1.x2(), V1.x3(), V1.x4(),
           V2.x1(), V2.x2(), V2.x3(), V2.x4(),
           V3.x1(), V3.x2(), V3.x3(), V3.x4(),
           V4.x1(), V4.x2(), V4.x3(), V4.x4())
  {}

  Rows(R3::XYZ V1, Real C1,
       R3::XYZ V2, Real C2,
       R3::XYZ V3, Real C3,
       R3::XYZ V4, Real C4) :
    Matrix(V1.x(), V1.y(), V1.z(), C1,
           V2.x(), V2.y(), V2.z(), C2,
           V3.x(), V3.y(), V3.z(), C3,
           V4.x(), V4.y(), V4.z(), C4)
  {}

};

class AugmentedRow {
public:

  Real mX[5];

public:

  AugmentedRow(){
    mX[0] = 0;
    mX[1] = 0;
    mX[2] = 0;
    mX[3] = 0;
    mX[4] = 0;
  }
  AugmentedRow(Real x1, Real x2, Real x3, Real x4, Real x5) {
    mX[0] = x1;
    mX[1] = x2;
    mX[2] = x3;
    mX[3] = x4;
    mX[4] = x5;
  }
  AugmentedRow(Column col, Real val){
    mX[0] = col.x(0);
    mX[1] = col.x(1);
    mX[2] = col.x(2);
    mX[3] = col.x(3);
    mX[4] = val;
  }

    
  AugmentedRow & operator*=(Real scalar){
    mX[0] *= scalar;
    mX[1] *= scalar;
    mX[2] *= scalar;
    mX[3] *= scalar;
    mX[4] *= scalar;
    return (*this);
  }

  AugmentedRow & operator-=(const AugmentedRow & other ){
    mX[0] -= other.mX[0];
    mX[1] -= other.mX[1];
    mX[2] -= other.mX[2];
    mX[3] -= other.mX[3];
    mX[4] -= other.mX[4];
    return (*this);
  }

  // ::::::::::::::::::::::::::::::::::::
  // ::: Access  (AugmentedRow Class) :::
  // ::::::::::::::::::::::::::::::::::::
  Real x(int i) const {return mX[i];}

  // :::::::::::::::::::::::::::::::::::::::::::
  // ::: Result Return  (AugmentedRow Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::
  const AugmentedRow ScaledBy(const Real scalar) const {
    return AugmentedRow(scalar * mX[0], scalar * mX[1],
                        scalar * mX[2], scalar * mX[3],
                        scalar * mX[4]);
    
  }


};

};// end namespace: R4
///
//

#endif //#ifndef GEOM_R4_H_
