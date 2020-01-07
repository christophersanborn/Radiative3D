// rtcoef.cpp
//
#include <iostream>
#include <iomanip>
#include <cstdlib>      /* rand(), RAND_MAX */
#include "rtcoef.hpp"

//////
// CONSTRUCTOR:  RTCoef Class
//
//   Sets up the basis coordinate system in which the outgoing rays
//   resulting from an R/T interface interaction will be expanded. The
//   basis vectors are (fpara, fparash, fnorm), which form an RHS
//   system isomorphic to (x, y, z) with z pointing perpendicularly
//   from the interface in the direction of the transmission side. The
//   outgoing ray direction will be a linear combination of fpara and
//   fnorm, whereas the outgoing PDOM (particle direction of motion,
//   aka polarization) will be a linear combination of fpara and fnorm
//   for P and SV outcomes, and will be fparash for an outcome of SH.
//
//   Inputs to the constructor are two unit vectors, (user is
//   responsible for ensuring unit-magnitude):
//
//     o  fnorm:  "face normal" - Unit vector normal to interface,
//                                pointing towards transmission side.
//
//     o  phdir:  "phonon direction" - Unit vector pointing in the ray
//                                propagation direction.
//
RTCoef::RTCoef(const R3::XYZ & fnorm_, const R3::XYZ & phdir)
  : NoTransmit(false), fnorm(fnorm_) {

  fpara = fnorm.GetInPlaneUnitPerpendicular(phdir);
                // Gets a perpendicular to fnorm that is in the
                // incidence plane and projects in same direction as
                // dir. (Basically, parallel to the interface, and
                // co-directional with incoming ray.)

  fparash = fnorm.Cross(fpara);
                // And get a mutual perpendicular. This is a basis
                // vector for SH-polarized particle motion.

  mSini = fpara.Dot(phdir);   // Sine of the incident angle (Used in
                              // Snell's Law)

  // NOTE: Subsequent functioning of this class depends on setting of
  // the Veloc and Density public members.  It is incumbant upon the
  // user of this class to set those members before using the methods
  // of the class to get an R/T outcome.
  //

}


//////
// METHOD:   RTCoef :: GetCoefs()
//
//   A dispatcher function to call the right helper function based on
//   incident raytype.
//
//   The result of this function (actually, of the appropriately
//   dispatched helper funciton) will be the setting of the
//   amplitudes, probabilities, and sines and cosines for all
//   permissible outcome raytypes.  Once these are set, the user can
//   call Choose() to actually select an outcome based on the relative
//   probabilities here determined.
//
//   A default choice is also recorded, for the cases where the
//   parameters result in indeterminant probabilities (such as
//   TotalP==0 or NaN).  In this case, the default should be to
//   reflect as same type as incoming type.  This can happen for
//   exactly-parallel incidence, and when one or more of the seismic
//   velocities is zero (a situation which should be avoided.  Density
//   of zero is OK though.)
//
void RTCoef::GetCoefs(raytype intype) {

  switch (intype) {
  case RAY_P:
    mDefChoice = R_P;
    GetCoefs_PSV(RAY_P);
    break;
  case RAY_SH:
    mDefChoice = R_SH;
    GetCoefs_SH();
    break;
  case RAY_SV:
    mDefChoice = R_SV;
    GetCoefs_PSV(RAY_SV);
    break;
  default:
    // Shouldn't ever get here.
    throw std::exception();
    break;
  }

}


//////
// METHOD:     RTCoef::GetCoefs_PSV()
// CALLED BY:  RTCoef::GetCoefs()
//
//   Gets coefficients AND probabilities for the R/T outcomes
//   allowable by the P-SV interaction.
//
void RTCoef::GetCoefs_PSV(raytype intype) {

                                        // Compute needed parameters
  GetSinesCosinesAndParamsPSV(intype);  // (Outgoing sines, cosines,
                                        // and the Aki and Richards
                                        // parameters from page 149 of
                                        // A&R 1980 kept in member
                                        // variable mAki.)
  // ::::
  // :: Notational Shorthand:
  // :

  const Real rho1 = DensityR;         // Useful aliases to match notation
  //const Real rho2 = DensityT;       // in Aki and Richards 1980 p. 149
  const Real alpha1 = VelocR[RAY_P];
  const Real alpha2 = VelocT[RAY_P];
  const Real beta1  = VelocR[RAY_S];
  const Real beta2  = VelocT[RAY_S];

  // ::::
  // :: Precompute some values:
  // :

  const Real p_sq = mAki.p*mAki.p;

  const Complex cosi1 = mCoso[R_P]  / alpha1;   // "i" for P-angles
  const Complex cosi2 = mCoso[T_P]  / alpha2;
  const Complex cosj1 = mCoso[R_SV] / beta1;    // "j" for SV-angles
  const Complex cosj2 = mCoso[T_SV] / beta2;

  const Real two = 2.0;

  // ::::
  // :: Compute R/T coefficients (relative amplitudes) with phase:
  // :

  Complex Term1;  // Reusable terms
  Complex Term2;  //

  if (intype == RAY_P) {

    // :: P\ -> P/  (Reflect as P)
    Term1 = ((mAki.b * cosi1) - (mAki.c * cosi2));
    Term2 = ((mAki.a) + (mAki.d * cosi1 * cosj2));
    mAmp[R_P] = (Term1 * mAki.F - Term2 * mAki.H * p_sq) / mAki.D;

    // :: P\ -> S/  (Reflect as S)
    Term1 = (mAki.a * mAki.b + mAki.c * mAki.d * cosi2 * cosj2);
    mAmp[R_SV] = -two * cosi1 * Term1 * mAki.p * alpha1 / (beta1 * mAki.D);

    // :: P\ -> P\  (Transmit as P)
    Term1 = two * rho1 * cosi1 * alpha1;   // (Also used in P\->S\)
    mAmp[T_P] =  Term1 * mAki.F / (alpha2 * mAki.D);

    // :: P\ -> S\  (Transmit as S)
    mAmp[T_SV] = Term1 * mAki.H * mAki.p / (beta2 * mAki.D);

  }
  else {  // (else intype is RAY_SV)

    // :: S\ -> P/  (Reflect as P)
    Term1 = (mAki.a*mAki.b + mAki.c*mAki.d*cosi2*cosj2);
    mAmp[R_P] = -two*cosj1*Term1*mAki.p*beta1 / (alpha1*mAki.D);

    // :: S\ -> S/  (Reflect as S)
    Term1 = (mAki.b*cosj1 - mAki.c*cosj2);
    Term2 = (mAki.a + mAki.d*cosi2*cosj1);
    mAmp[R_SV] = -(Term1*mAki.E - Term2*mAki.G*p_sq) / mAki.D;

    // :: S\ -> P\  (Transmit as P)
    Term1 = two*rho1*cosj1*beta1;      // (Also used in S\->S\)
    mAmp[T_P] = -Term1*mAki.G*mAki.p / (alpha2*mAki.D);

    // :: S\ -> S\  (Transmit as S)
    mAmp[T_SV] = Term1*mAki.E / (beta2*mAki.D);

  }

  // ::::
  // :: Compute Probabilities from Amplitudes:
  // :

  GetProbabilitiesFromAmplitudesPSV();

            // Probabilities are now set, and can be used for choosing
            // an outcome.  Phase advance for a given outcome can be
            // determined by looking at the complex phase of the
            // amplitude.  Sines and Cosines are set, and the
            // direction vector of a given outcome can be constructed
            // from the real-parts thereof.  Oooh-rah!
            //
}


//////
// METHOD:     RTCoef::GetCoefs_SH()
// CALLED BY:  RTCoef::GetCoefs()
//
//   Handles the case of an incident SH wave.
//
void RTCoef::GetCoefs_SH() {

  // ::::
  // :: The following outcomes do not couple to incoming SH, and can
  // :: be ignored:
  // ::

  mProb[R_P]  = 0;    // Set these output probabilities
  mProb[R_SV] = 0;    // to zero
  mProb[T_P]  = 0;    //
  mProb[T_SV] = 0;    //

  // ::::
  // :: Notational Shorthand:
  // :

  const Real &  rho1 = DensityR;        // Densities
  const Real &  rho2 = DensityT;
  const Real & beta1 = VelocR[RAY_S];   // Velocities
  const Real & beta2 = VelocT[RAY_S];
  const rt_result_e J1 = R_SH;          // Convention of using J for 
  const rt_result_e J2 = T_SH;          // S-angle and I for P-angles

  const Real two = 2.0;

  // ::::
  // :: Get outgoing sines and cosines:
  // :

  mSino[J1] = mSini;
  mSino[J2] = (beta2/beta1)*mSini;

  mCoso[J1] = sqrt(Complex(1.0 - mSino[J1]*mSino[J1]));
  mCoso[J2] = sqrt(Complex(1.0 - mSino[J2]*mSino[J2]));
                      // Transmitted cosine (J2) will be pure imaginary
                      // if transmitted sine is > 1.0.

  // ::::
  // :: Calculate reflection coefficient and phase:
  // :

  Complex a = rho1*beta1*mCoso[J1];
  Complex b = rho2*beta2*mCoso[J2];

  mAmp[R_SH] = (a - b) / (a + b);
  mAmp[T_SH] = two*a / (a + b);

  // ::::
  // :: Compute Probabilities:
  // :
  // :     P ~~ rho*vel * Re(Cos(angle)) * || Amplitude || ^ 2
  //
  //         I.e., probabilities are proportional to the
  //         magnitude-squared of the amplitude, scaled by the density
  //         and velocity in the medium the ray propagates in, and
  //         scaled by the REAL component of the cosine (this will
  //         zero-out post-critically excluded rays).
  //
  //         In Aki and Richards, it is also scaled by a factor
  //         quantifying the fraction of energy flux from the incoming
  //         wave crossing the interface.  This scaling effectively
  //         normalizes the probabilities and assures that they sum to
  //         1.0. But we don't need normalized probabilities, since
  //         they get normalized anyway in the random chooser, so we
  //         ignore that scaling (and save a few cpu cycles in the
  //         process).
  //

  mProb[R_SH] = rho1*beta1*mCoso[J1].real() * norm(mAmp[R_SH]);
  mProb[T_SH] = rho2*beta2*mCoso[J2].real() * norm(mAmp[T_SH]);

}


//////
// METHOD:     RTCoef::GetSinesCosinesAndParamsPSV()
// CALLED BY:  RTCoef::GetCoefs_P()
//             RTCoef::GetCoefs_SV()
//
//   Computes a bunch of factors needed in P-SV interaction R/T
//   coefficients.  Intype should be one of RAY_P or RAY_SV.
//
void RTCoef::GetSinesCosinesAndParamsPSV(raytype intype) {

  const Real rho1 = DensityR;         // Useful aliases to match notation
  const Real rho2 = DensityT;         // in Aki and Richards 1980 p. 149
  const Real alpha1 = VelocR[RAY_P];
  const Real alpha2 = VelocT[RAY_P];
  const Real beta1  = VelocR[RAY_S];
  const Real beta2  = VelocT[RAY_S];

  const raytype
    rtyp = (intype==RAY_P) ? RAY_P    // Squash SV/SH to
                           : RAY_S;   // just plain S

  // Horizontal Slowness:

  mAki.p = mSini / VelocR[rtyp];

  // Compute Outgoing Sines from Incident Sine:

  mSino[T_P]  = VelocT[RAY_P]  * mAki.p;
  mSino[T_SV] = VelocT[RAY_S] * mAki.p;
  mSino[R_SV] = VelocR[RAY_S] * mAki.p;
  mSino[R_P]  = VelocR[RAY_P]  * mAki.p;

  // Compute Outgoing Cosines:

  mCoso[T_P]  = sqrt(Complex(1.0 - mSino[T_P]*mSino[T_P]));
  mCoso[T_SV] = sqrt(Complex(1.0 - mSino[T_SV]*mSino[T_SV]));
  mCoso[R_SV] = sqrt(Complex(1.0 - mSino[R_SV]*mSino[R_SV]));
  mCoso[R_P]  = sqrt(Complex(1.0 - mSino[R_P]*mSino[R_P]));

                    // For sines that are > 1.0, the cosines
                    // will be imaginary.

  // Compute the remaining parameters:
  //     These ones are all real:

  Real beta1_sq = beta1*beta1;
  Real beta2_sq = beta2*beta2;
  Real p_sq = mAki.p*mAki.p;

  Real tmp1 = rho1 * (1. - 2. * beta1_sq * p_sq);
  Real tmp2 = rho2 * (1. - 2. * beta2_sq * p_sq);
  Real tmp3 = 2. * rho1 * beta1_sq;
  Real tmp4 = 2. * rho2 * beta2_sq;


  mAki.a = tmp2 - tmp1;
  mAki.b = tmp2 + tmp3 * p_sq;
  mAki.c = tmp1 + tmp4 * p_sq;
  mAki.d = tmp4 - tmp3;

  //     And these ones are complex:

  Complex cosi1 = mCoso[R_P]  / alpha1;   // "i" for P-angles
  Complex cosi2 = mCoso[T_P]  / alpha2;
  Complex cosj1 = mCoso[R_SV] / beta1;    // "j" for SV-angles
  Complex cosj2 = mCoso[T_SV] / beta2;

  mAki.E = mAki.b * cosi1 + mAki.c * cosi2;
  mAki.F = mAki.b * cosj1 + mAki.c * cosj2;

  mAki.G = mAki.a - mAki.d * cosi1 * cosj2;
  mAki.H = mAki.a - mAki.d * cosi2 * cosj1;

  mAki.D = mAki.E * mAki.F + mAki.G * mAki.H * p_sq;

}


//////
// METHOD:     RTCoef::GetProbabilitiesFromAmplitudesPSV()
// CALLED BY:  RTCoef::GetCoefs_P()
//             RTCoef::GetCoefs_SV()
//
//   Computes the relative probabilites of each output type allowed by
//   the P-SV interaction.  Follows equation 5.40 in Aki and Richards
//   1980 (p. 151).  Probabilities are not normalized, (sum is not
//   necessarily 1.0), but this doesn't matter, as the random chooser
//   does not expect or require normalized probabilities.
//
void RTCoef::GetProbabilitiesFromAmplitudesPSV() {

  mProb[R_SH] = 0;    // Set these output probabilities to zero,
  mProb[T_SH] = 0;    // as the SH outcomes do not couple with P
                      // or SV.

  const Real rho1 = DensityR;   // Useful aliases to match notation
  const Real rho2 = DensityT;   // in Aki and Richards 1980 p. 149
  const Real alpha1 = VelocR[RAY_P];
  const Real alpha2 = VelocT[RAY_P];
  const Real beta1  = VelocR[RAY_S];
  const Real beta2  = VelocT[RAY_S];

  Real cosi1 = mCoso[R_P].real();   // The real components determine
  Real cosi2 = mCoso[T_P].real();   // energy flux vector.
  Real cosj1 = mCoso[R_SV].real();
  Real cosj2 = mCoso[T_SV].real();

  mProb[R_P]  = rho1*alpha1*cosi1 * mAmp[R_P].norm();
  mProb[R_SV] = rho1*beta1*cosj1  * mAmp[R_SV].norm();
  mProb[T_P]  = rho2*alpha2*cosi2 * mAmp[T_P].norm();
  mProb[T_SV] = rho2*beta2*cosj2  * mAmp[T_SV].norm();

}


//////
// METHOD:   RTCoef :: ChooseSPolType()
//
//   Given a polarization direction, choose SH or SV based on
//   fractional representation.
//
//   Works by comparing pdom (particle direction-of-motion, aka
//   polarization) against fparash, which is assumed to have been
//   previously initialized prior to this function being called.
//
raytype RTCoef::ChooseSPolType(const R3::XYZ & pdom) const {

  raytype  stype;   // Will be either RAY_SH or RAY_SV
  Real    shfrac;   // Fraction of ray energy that is SH wrt interface

  shfrac = pdom.Dot(fparash);
  shfrac *= shfrac;

  Real ran = ((Real) rand() / RAND_MAX);  // 0.0 < ran <= 1.0
  if (ran <= shfrac)
     { stype = RAY_SH; }
  else
     { stype = RAY_SV; }

  return stype;

}


//////
// METHOD:   RTCoef :: Choose()
//
//   Makes a random pick of the output ray based on the relative
//   probabilities coded in mProb[].  Save that pick in mRTChoice,
//   which is then used by reporting methods to report back various
//   aspects of the chosen outcome.
//
//   It is assumed that the mProb[] member has already been populated
//   by a call to the public GetCoefs() method.
//
void RTCoef::Choose() {

  Real PI[RT_NUM];        // Integrated Probability
  Real TotalP;

  PI[0] = mProb[0];
  for (Index i = 1; i < RT_NUM; i++) {
    PI[i] = PI[i-1] + mProb[i];
  }
  TotalP = PI[RT_NUM-1];

  int iran = rand();      // Generate an int in range: [0,RAND_MAX]
  if (iran==0) {iran=1;}  // Now we are in the range:  (0,RAND_MAX]
  Real ran = ((Real) iran / RAND_MAX) * TotalP; // Range: (0.0, TotalP]

  Index choice = RT_NUM-1;
  for (Index i = 0; i < (RT_NUM-1); i++) {
    if (ran <= PI[i]) {
      choice = i;
      break;
    }
  }

  if ((TotalP==0) ||              // If TotalP is zero
      ((TotalP - TotalP)!=0)) {   // or not-finite (i.e., inf or nan)
    choice = mDefChoice;          // Then take default rather than 
  }                               // what the loop chose.  (Handles
                                  // cases where probabilities are
                                  // indeterminant.)

  if (NoTransmit) { // Check for case where user does not want transmission
    if (choice==T_P) {choice=R_P;}    // (used e.g. for free-surface
    if (choice==T_SV) {choice=R_SV;}  // modelling)
    if (choice==T_SH) {choice=R_SH;}  //
  }

  mChoice = (rt_result_e) choice;   // Save result so the reporting
                                    // methods can report on it.

}


//////
// METHOD:   RTCoef :: GetChosenRaytype()
//
//   Returns the raytype of the chosen outcome. We return either RAY_P
//   or RAY_S, and do not distinguish between SV and SH, as these
//   terms are used internally relative to the interface plane and
//   might not have the same meaning to the calling function (which
//   likely defines them w.r.t. an absolute coordinate system or some
//   other reference).
//
//   It is assumed that the Choose() method (which sets mChoice) has
//   already been called.
//
raytype RTCoef::GetChosenRaytype() const {

  if ((mChoice==R_P)||(mChoice==T_P)) {
    return RAY_P;}
  else {
    return RAY_S;
  }

}


//////
// METHOD:   RTCoef :: DidRayTransmit()
//
//   Returns true if the chosen outcome represents a transmitted ray,
//   otherwise returns false for a relfected ray.  It is assumed that
//   the Choose() method (which sets mChoice) has already been called.
//
bool RTCoef::DidRayTransmit() const {

  if ((mChoice==R_P)||(mChoice==R_SV)||(mChoice==R_SH)) {
    return false;}
  else {
    return true;
  }

}


//////
// METHOD:   RTCoef :: GetChosenRayDirection()
//
//   Returns a unit vector pointed in the direction of the outgoing
//   ray.  This method assumes that the Choose() method (which sets
//   mChoice) has already been called.  The result of this method is
//   saved in mChosenRayDir which can then be used by a subsequent
//   call to GetChosenParticleDOM().
//
R3::XYZ RTCoef::GetChosenRayDirection() {

  Real comp_para;     // Component parallel to interface
  Real comp_norm;     // Component normal to interface

  comp_para = mSino[mChoice];
  comp_norm = mCoso[mChoice].real();  // Real part determines direction

  if (comp_para > 1.0) {comp_para=1.0;}

  if ((mChoice==R_P)||(mChoice==R_SV)||(mChoice==R_SH)) {
    comp_norm *= -1;
  }

  mChosenRayDir = fpara.ScaledBy(comp_para)
                + fnorm.ScaledBy(comp_norm);

  return mChosenRayDir;

}


//////
// METHOD:   RTCoef :: GetChosenParticleDOM()
//
//   Returns direction-of-particle-motion for the chosen ray
//   outcome. (Essentially the polarization direction.)  This method
//   assumes that the GetChosenRayDirection() method (which sets
//   mChosenRayDir) has already been called.
//
R3::XYZ RTCoef::GetChosenParticleDOM() {

  R3::XYZ dopm;   // Direction of Particle Motion

  if ((mChoice==T_P)||(mChoice==R_P)) {
    dopm = mChosenRayDir;   // (longitudinal polarization aka P-wave)
  }
  else if ((mChoice==T_SH)||(mChoice==R_SH)) {
    dopm = fparash;         // SH-polarization
  } 
  else if (mChoice==R_SV) {
    dopm = mChosenRayDir.Cross(fparash);  // SV-polarization,
  }                                       // for reflected ray
  else { /* if (mChoice==T_SV) */
    dopm = fparash.Cross(mChosenRayDir);  // SV-polarization,
  }                                       // transmitted ray

                // NOTE: The sign convention for the SV polarization
                // vector is diagrammed in Aki & Richards 1980 fig
                // 5.9 and explained in fig 5.5.  The distinction
                // between the reflected and transmitted case is
                // irrelevant for computing seismic envelopes, but
                // becomes important if we want to produce full
                // waveform seismograms, as does phase tracking
                // (c.f. GetChosenRayPhase()).
                //

  return dopm;

}


//////
// METHOD:   RTCoef :: GetChosenRayPhase()
//
//   Returns the phase-shift that the ray experiences via the R/T
//   interaction.  Works by returning the complex phase of the
//   amplitude coefficient (mAmp[]) of teh chosen outcome.  This
//   method assumes that the Choose() method (which sets mChoice) has
//   already been called.
//
Real RTCoef::GetChosenRayPhase() const {

  // NOTE: In order for this to work properly, special care must be
  // taken for incoming S waves when measuring SV and SH-fraction,
  // because it may be necessary to add a half-cycle shift based on
  // whether the incoming PDOM projects positively or negatively with
  // the SV/SH basis vectors that we use internally. Until this is
  // properly handled, this method is deliberately coded to break,
  // rather than return a result, in order to prevent a careless
  // application of this feature.  This one needs a bit of careful
  // thinking-through.  /* TODO: FIX */
  //

  std::cerr << "Not implemented yet.\n";
  throw std::exception();

}


//////
// METHOD:   RTCoef :: PrintChosenRaytype()
//
//   Prints probabilities and outcome to stdout. Used for diagnostic
//   purposes.
//
void RTCoef::PrintChosenRaytype() const {

  int width = 14;
  static bool header = false;

  if ((!header) || (mSini==0.)) {
    header = true;
    std::cout << "##"
              << std::setw(width) << "Sine_in"
              << std::setw(width) << "Prob_R_P"
              << std::setw(width) << "Prob_T_P"
              << std::setw(width) << "Prob_R_SV"
              << std::setw(width) << "Prob_T_SV"
              << std::setw(width) << "Prob_R_SH"
              << std::setw(width) << "Prob_T_SH"
              << "    Result" << std::endl << std::endl;
  }

  std::cout << "  "
            << std::setw(width) << mSini
            << std::setw(width) << mProb[R_P]
            << std::setw(width) << mProb[T_P]
            << std::setw(width) << mProb[R_SV]
            << std::setw(width) << mProb[T_SV]
            << std::setw(width) << mProb[R_SH]
            << std::setw(width) << mProb[T_SH]
            << "    ";

  switch (mChoice) {
  case R_P:
    std::cout << "Reflection/P" << std::endl;
    break;
  case R_SV:
    std::cout << "Reflection/SV" << std::endl;
    break;
  case R_SH:
    std::cout << "Reflection/SH" << std::endl;
    break;
  case T_P:
    std::cout << "Transmission/P" << std::endl;
    break;
  case T_SV:
    std::cout << "Transmission/SV" << std::endl;
    break;
  case T_SH:
    std::cout << "Transmission/SH" << std::endl;
    break;
  default:
    // Shouldn't ever get here.
    throw std::exception();
    break;
  }

}


//////
// METHOD:   RTCoef :: RunRTCoefTest()   (STATIC)
//
//   Run through a set of input angles for a given set of material
//   parameters and output the resultant set of outcome probabilities.
//
void RTCoef::RunRTCoefTest(Count nSini,
                           Real rho1, Real alpha1, Real beta1,
                           Real rho2, Real alpha2, Real beta2) {

  // ::::
  // :: Step through the three incident raytypes:
  //:
  for (Index irt = 0; irt < RAY_NUMTYPES; irt++) {

    // ::::
    // :: Output a Header Block:
    // :
    std::cout << std::endl
              << "##RTCoef Probability Test:"
              << std::endl << "##" << std::endl;
    std::cout << "##       Input Raytype:  ";
    switch(irt) {
    case RAY_P:
      std::cout << "RAY_P"; break;
    case RAY_SV:
      std::cout << "RAY_SV"; break;
    case RAY_SH:
      std::cout << "RAY_SH"; break;
    }
    std::cout << std::endl << "##" << std::endl
              << "##     Reflection Side:  (rho,alpha,beta) = ( "
              << rho1 << " " << alpha1 << " " << beta1 << " )"
              << std::endl
              << "##   Transmission Side:  (rho,alpha,beta) = ( "
              << rho2 << " " << alpha2 << " " << beta2 << " )"
              << std::endl << "##" << std::endl;


    // ::::
    // :: Range through nSini incidence angles:
    // :
    for (Index isin = 0; isin < nSini; isin++) {
      Real theta = Geometry::Pi90 * ((Real) isin / (nSini-1));
      S2::ThetaPhi norm(0,0);           // Points UP
      S2::ThetaPhi incidence(theta,0);  // Points at an angle
      RTCoef rt(norm, incidence);
      rt.DensityR = rho1;
      rt.DensityT = rho2;
      rt.VelocR[RAY_P] = alpha1;
      rt.VelocT[RAY_P] = alpha2;
      rt.VelocR[RAY_S] = beta1;
      rt.VelocT[RAY_S] = beta2;
      rt.GetCoefs((raytype) irt);
      rt.Choose();
      rt.PrintChosenRaytype();
    }

  }


}
