// rtcoef.hpp
///@file
///
/// This file develops the RTCoef class, which encapsulates
/// Reflection/Transmission coefficients across a first-order
/// discontinuous locally flat interface, and the probabilities of
/// transmission and reflection into P and S raytypes.
///
/// The theoretical basis for this is developed in Aki and Richards
/// 1980, "Quantitative Seismology", Chapter 5.
///
/// AUTHOR NOTE: While the code in this class is entirely my own work,
///              I benefitted tremendously from being able to inspect
///              the source code to Peter Shearer's PSPhonon, which
///              implements the Aki and Richards treatment of
///              Reflection/Transmission in much the same way as is
///              done here.  ~C. Sanborn
///
///
#ifndef RTCOEF_H_
#define RTCOEF_H_
//
#include "geom.hpp"
#include "raytype.hpp"
#include "complex.hpp"

//////
// CLASSES: Definitions
//
// INCLUDING:
//
//   o  class RTCoef
//

//////
// CLASS:   ::::  RTCoef  ::::
///
///   Encodes a set of probabilities for transmitting or reflecting
///   through an interface, and conversions between P/S polarizations.
///
///   The structure has elements for all possible conversions, but only
///   the ones relevant to a given input type should be non-zero.
///
/// USAGE:
///
///   To use this class, a user must:
///
///   1.)  Construct an RTCoef object, providing an interface normal
///        and a ray direction.
///
///   2.)  Set the public velocity and density members, which quantify
///        the discontinuity across the interface.
///
///   2a.) To model a lossless free-surface: Set tranmission-side
///        density to zero, and transmission-side velocities to
///        very-small numbers (to avoid NaNs). This actually seems to
///        result in zero transmission probability, but, to be certain
///        (in case of numerical inacuracies), also set the public
///        NoTransmit member to 'true'.
///
///   3.)  Call GetCoefs(), providing a raytype (P, SV, or SH), which
///        will compute a set of outcome amplitudes and probabilities
///        based on the provided incident raytype.  (If necessary to
///        the use case, the user may call ChooseSPolType() to make a
///        selection between SV and SH based on the incident
///        polarization relative to the interface.)
///
///   4.)  Call Choose(), which will choose an outcome at random based
///        on the computed relative probabilities.
///
///   5.)  Call the various GetChosen...() methods to retrieve the
///        details of the chosen ray outcome.
///
///
class RTCoef {
private:

  // :::::::::::::::::::::::::::::
  // ::: Enums  (RTCoef Class) :::
  // :::::::::::::::::::::::::::::

  enum rt_result_e {  // Identifies outcomes of a R/T event
    R_P,              //   Reflected P
    R_SV,             //   Reflected SV
    R_SH,             //   Reflected SH
    T_P,              //   Transmitted P
    T_SV,             //   Transmitted SV
    T_SH,             //   Transmitted SH
    RT_NUM            // Number of outcome types (six)
  };

  // ::::::::::::::::::::::::::::::::::::::
  // ::: Member Classes  (RTCoef Class) :::
  // ::::::::::::::::::::::::::::::::::::::

  struct AkiParamsPSV {
    Real      p;    // Horizontal slowness

    Real      a;    // These ones depend on the rho's, beta's and p
    Real      b;    //
    Real      c;    //    (See Aki and Richards 1980 page 149)
    Real      d;    //

    Complex   E;    // Cosine-dependent terms:
    Complex   F;    //
    Complex   G;    //    (See Aki and Richards 1980 page 149)
    Complex   H;    //

    Complex   D;    // Determinant (See A&R p. 149)
  };


public:

  // :::::::::::::::::::::::::::::::::::
  // ::: Member Data  (RTCoef Class) :::
  // :::::::::::::::::::::::::::::::::::
  //
  //        Access level: Public  (this class used as a
  //                               semi-smart structure)
  //
  //        These members must be set externally before the "smart"
  //        features of this class can be used.  Think of them almost
  //        as parameters that you would pass to a constructor, just
  //        without the formalism of a constructor.  Responsibility
  //        for the validity of these data lies with the external
  //        user. No internal validity checks are performed on these
  //        public data members.
  //
  //        (In the current incarnation of things, these members are
  //        set by the CellFace::GetRTBasis() method.  Look there for
  //        more details.)
  //

  Real  VelocR[RAY_NBT];  // Velocities on the Reflection side
  Real  VelocT[RAY_NBT];  // Velocities on the Transmission side
  Real         DensityR;  // Densities on the Reflection Side
  Real         DensityT;  // Densities on the Transmission Side

  bool       NoTransmit;  // If true, transmission is prohibited by a
                          // post-choice check that flips transmit
                          // results to reflection results.

private:

  // :::::::::::::::::::::::::::::::::::::::::::
  // ::: Private Member Data  (RTCoef Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::

  // :
  // ::             Basis Vectors:
  // :

  R3::XYZ       fnorm;  // Outward normal to interface
  R3::XYZ       fpara;  // Parallel to interface, in incidence plane
  R3::XYZ     fparash;  // Parallel to iface, perp to incidence plane.

                        //     The above three vectors form the basis
                        //     in which the resultant ray direction
                        //     will be expanded, with fnorm pointing
                        //     towards the transmission side of the
                        //     interface, and fpara pointing along the
                        //     parallel projection of the incoming
                        //     (incident) ray.
                        //   


  // :
  // ::             Sines and Cosines:
  // :

  Real              mSini;  // Sine of incident angle

  Real      mSino[RT_NUM];  // Sine of outgoing angle (c.f. sini)
                            //      Can be > 1.0, in which case
                            //      corresponding cosine will be
                            //      imaginary.

  Complex   mCoso[RT_NUM];  // Cosine of outgoing angle
                            //      Will be real when sine <= 1.0,
                            //      imaginary otherwise.
                            //      Public access function Coso()
                            //      returns real part, suitable for
                            //      construcing ray vectors.

  // :
  // ::             Intermediate Parameters:
  // :

  AkiParamsPSV       mAki;  // Aki and Richards parameters for the
                            // P-SV interaction.

  // :
  // ::             Amplitudes and Probabilities:
  // :

  Complex    mAmp[RT_NUM];  // Relative amplitude of outgoing wave.
                            //      Relative to incoming wave.
                            //      Complex phase encodes phase
                            //      advance from interaction with
                            //      interface.

  Real      mProb[RT_NUM];  // Probability of each resultant wave.
                            //      Only allowed result types are
                            //      non-zero: Those that are not
                            //      relevant to the incoming type, or
                            //      which are precluded by
                            //      post-critical incidence, are zero.
                            //      Total probability is not
                            //      necessarily normalized to 1.0, and
                            //      it doesn't need to be, as the
                            //      random chooser will effectively
                            //      determine the norm when it
                            //      integrates the probabilites to
                            //      make a choose map.
                            // 

  // :
  // ::             Choice Registers:  (Record the choice made by the 
  // :                                  chooser function.)

  rt_result_e     mChoice;  // Chosen outcome type.  This member will
                            // be set by Choose() and interpreted by
                            // the various reporting methods that
                            // inform the outside world about the
                            // result.

  rt_result_e  mDefChoice;  // Default choice. If probabilities result
                            // in an undefined choice (because TotalP
                            // is either zero or NaN), then this is
                            // the result that we return.  Set by the
                            // GetCoefs() method.

  R3::XYZ   mChosenRayDir;  // Direction of the chosen output ray,
                            // as computed and stored here by the
                            // GetChosenRayDirection() method.


public:

  // :::::::::::::::::::::::::::::::::::
  // ::: Constructor  (RTCoef Class) :::
  // :::::::::::::::::::::::::::::::::::

  RTCoef(const R3::XYZ & fnorm_, const R3::XYZ & phdir);
                // Sets up the basis vectors and the incident sine.
                // The public members (velocities and densities) still
                // need to be set by the user before any of the
                // features of the class can be used.


  // :::::::::::::::::::::::::::::::::::::::::::
  // ::: Computation Methods  (RTCoef Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::

  void GetCoefs(raytype intype);  // Dispatcher for _SH/_PSV below


private:

  void GetCoefs_PSV(raytype);   // Compute output coefficients and
  void GetCoefs_SH();           //  probabilities for all cases
                                //  allowed by the input raytype

  void GetSinesCosinesAndParamsPSV(raytype);  // Helper function to compute
                                              // params needed by the P and
                                              // SV interactions.

  void GetProbabilitiesFromAmplitudesPSV();   // Helper function to compute
                                              // the probabilities from the
                                              // amplitudes in the P-SV
                                              // case.


public:

  // ::::::::::::::::::::::::::::::::::::::::
  // ::: Choose Functions  (RTCoef Class) :::
  // ::::::::::::::::::::::::::::::::::::::::
  // 
  //                These functions make statistical decisions from
  //                probability distributions and random-number
  //                generation.
  //

  raytype ChooseSPolType(const R3::XYZ & pdom) const;
                    // Given a polarization direction, choose SH or SV
                    // based on fractional representation.

  void Choose();    // Select an outcome based on the relative
                    // probabilities coded in mProb[].

  raytype GetChosenRaytype() const;   // Reveal the raytype of the
                                      // choice made in Choose().

  bool DidRayTransmit() const;        // Reveals whether chosen outcome
                                      // was a transmission (as opposed
                                      // to a reflection).

  R3::XYZ GetChosenRayDirection();    // Reveals the direction of the ray
                                      // outcome chosen in Choose().

  R3::XYZ GetChosenParticleDOM();     // Reveals particle motion
                                      // (polarization direction) of the
                                      // ray outcome chosen in Choose()
                                      // (for S-waves).

  Real GetChosenRayPhase() const;     // Reveals phase-shift of the
                                      // chosen outgoing ray.


  // :::::::::::::::::::::::::::::::::::::::::::::::::
  // ::: Diag and Output Functions  (RTCoef Class) :::
  // :::::::::::::::::::::::::::::::::::::::::::::::::

  void PrintChosenRaytype() const;    // Prints probabilities and
                                      // outcome to stdout. (Used for
                                      // diagnostic purposes.)

  static void RunRTCoefTest(Count nSini, 
                            Real rho1, Real alpha1, Real beta1,
                            Real rho2, Real alpha2, Real beta2);
                // Run through nSini input angles (sines) with the
                // given densities and velocities and output the
                // resultant probabilities.
                            

};


///
#endif //#ifndef RTCOEF_H_
//
