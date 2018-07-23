//  takeoff.hpp
/// @file
///
/// Defines classes related to Take-off Angles (discretizations of the
/// initial propagation directions of Phonon objects from PhononSource
/// objects). Class TakeoffAngle represents a single, discrete take-off
/// angle.  Class TakeoffArray encapsulates a set of take-off angles
/// covering a unit sphere.
///
#ifndef TAKEOFF_H_
#define TAKEOFF_H_
//
#include <vector>
#include "geom.hpp"

//////
// CLASS: TakeoffAngle
///@brief
///
/// Encapsulates a single take-off angle, representing a possibile
/// inital propagation direction for a ray or Phonon object departing
/// from a PhononSource object.  Includes a weighting factor that can
/// be used to account for variation in the size of the angular
/// element represented by each individual ToA.
///
/// Note: the exact meaning of the weigting factor is to be determined
/// and enforced by the TakeoffArray class, which will decide whether
/// the set of all weights will sum to 4*pi, 1.0, or an arbitrary sum.
///
class TakeoffAngle {
protected:
;
  S2::ThetaPhi  mDir;     ///<  Direction of take-off
  Real          mWeight;  ///<  Weight factor

public:
;
  TakeoffAngle(S2::ThetaPhi toa, Real weight) 
    : mDir(toa), 
      mWeight(validate::positive(weight)) {}

  Real Weight() const {return mWeight;}
  void SetWeight(Real weight) {mWeight=validate::positive(weight);}

};//
///

//////
// CLASS: TakeoffArray
///@brief
///
/// Keeps a vector of TakeoffAngles and associated solid angular
/// elements covering a unit sphere. Represents a discretization of
/// initial directions for point sources.
///
class TakeoffArray {
protected:
  ;
  Count mN;
  std::vector<TakeoffAngle> mvTOA;

public:
  ;
  TakeoffArray(Count size);

  int size() const {return mvTOA.size();}

protected:
  ;

};//
///

//////
// CLASS: NewTesselSphere
///@brief
///experimental
class NTesselSphere : TakeoffArray {
public:
  ;

  NTesselSphere(int degree);

  static int NumAnglesFromDegree(int);

};//
///

///
#endif //#ifndef TAKEOFF_H_
//
