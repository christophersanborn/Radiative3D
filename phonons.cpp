// phonons.cpp
//
#include <iostream>
#include <cmath>
#include <cstdlib>      /* rand(), RAND_MAX */
#include <cassert>
#include "phonons.hpp"
#include "media.hpp"
#include "rtcoef.hpp"
#include "scatterers.hpp"
#include "dataout.hpp"

//
// CLASS IMPLEMENTATION:  Phonon
//
//  Methods Defined Here:
//
//   o  Move()
//   o  Transform()
//   o  PseudoReflect()
//   o  DirectionOfMotion()
//   o  TransferWithRefraction()
//   o  InsertInto()
//   o  Propagate()
//

//////
// CLASS-STATIC PARAMETERS:
//

unsigned long Phonon::cm_phonon_counter = 0;
unsigned long Phonon::cm_loop_concern = 1048576; // (2^20, ~1M loops)
Real Phonon::cm_slow_concern = 0.001;
Real Phonon::cm_min_theta = 0.0000001;
Real Phonon::cm_max_theta = Geometry::Pi - 0.0000001;
Real Phonon::cm_ttl = 3600.0;    // One hour.  Note: the Model()
                                 // constructor is reponsible for
                                 // setting this value according to
                                 // user preferences. So to affect a
                                 // different "default" value, make
                                 // the change there, not here.

//////
// METHOD:  Phonon :: Move()    (private, at the moment)
//
// WHAT-IT-DOES:
//
//   It "moves" a phonon, updating it's internal location, direction,
//   total pathlength and traveltime members, etc., according to the
//   Travel Record provided in the argument.
//
// ACCESS:
//
//   Method is private, at the moment, because it is only used
//   internally by the Propagate() method.  At present I can think of
//   no need for outside accessibility to this method, but it is
//   conceivable that that could change.  Also, I mark it 'inline' as
//   it is only used in this translation unit in the Propagate() inner
//   loop.  Again though, this could change if I opt for public
//   accessibility of this method.
//
inline void Phonon::Move(const TravelRec & travel) {
  mPathLength += travel.PathLength;
  mTimeAlive  += travel.TravelTime;
  mRecentTravelTime += travel.TravelTime;
  mLoc         = travel.NewLoc;
  mDir         = travel.NewDir;
  mAmplitude  *= travel.Attenuation;
  mMoveCount +=1;
}


//////
// METHOD:  Phonon :: Transform()
//
// WHAT-IT-DOES:
//
//   Transforms the direction and polarization of the phonon by an
//   amount determined by a "relative phonon" taken as an input.
//
// HOW-IT-WORKS:
//
//   We assume the current direction and polarization of the phonon
//   define a coordinate system {a} wherin a3 is a unit vector in the
//   direction of travel, a1 is in the direction of polarization, and
//   a2 is mutually perpendicular as in an rhs coordinate system.
//
//   We further assume that the relative phonon (typically generated
//   by a scattering algorithm) is expressed in a coordinate system
//   {s} wherin s3 is a unit vector in the direction that represents
//   zero deviation from the original path, s1 represents the original
//   polarization direction, and s2 is mutually perpendictular.
//
//   Internally, the phonon direction is represented in an independent
//   "lab frame".  This method updates the phonon direction and
//   polarization to account for the deflection/rotation coded in the
//   relative phonon.
//
//   This method borrows heavily from the code in the SCATRAYPOL
//   and BENDRAY subroutines in Shearer's PSPhonon.
//
// INPUTS:
//
//   rel   -  A "relative" phonon (direction, type, and polarization
//            code a deflection/rotation and/or P<->S conversion)
// OUTPUTS:
//   (none)
//
// SIDE-AFFECTS:
//
//   Internal object state is changed:
//     * dir updated to new phonon direction.
//     * pol updated to new polarization.
//     * type updated to new raytype
//
void Phonon::Transform(const Phonon & rel) {
  
  // Establish coordinate system {a1, a2, a3} to represent the CURRENT
  // direction and orientation of phonon, using a3 as the unit vector
  // in the current direction of travel, a1 as the unit vector in the
  // current S polarization direction, and a2 mutually perpendicular:

  R3::OrthoAxes AA(mDir.Theta(), mDir.Phi(), mPol);   // Magic!

                // AA.S1()  <-->  Direction of S-Polarization (a1) 
                // AA.S2()  <-->  (a2)
                // AA.S3()  <-->  Direction of travel (a3)
                // AA.E1()  <-->  Theta^Hat direction (for reference)
                // AA.E2()  <-->  Phi^Hat direction (for reference)
                //
                //     Vector components of set AA are expressed in
                //     the "lab frame".
                //

  // Establish set of basis vectors {b1, b2, b3} to represent the
  // deflection direction and orientation coded in the relative
  // phonon.  (Ie, this axes set represents the scattering angles
  // given to us in 'Phonon& rel')

  R3::OrthoAxes BB(rel.mDir.Theta(), rel.mDir.Phi(), rel.mPol);

                // BB.S3()  <-->  Deflection direction relative to AA
                // BB.S1()  <-->  Deflected polarization direction,
                //                relative to AA
                //
                //     Vector components of set BB are expressed, as
                //     it were, in the "AA frame".
                //

  // Now establish a set of basis vectors {s1, s2, s3} to represent
  // the SCATTERED phonon direction and orientation, properly
  // expressed in the lab frame: (This is essentially a basis
  // transform of BB from the AA frame to the Lab frame.)

  R3::OrthoAxes SS = AA.Express(BB);  

  // Now collect the new phonon direction and orientation from the {s}
  // OrthoAxes set:

  mDir.Set(SS.Theta(), SS.Phi()); // New direction
  mPol = SS.Rot();                // New orientation
  mType = rel.mType;              // New polarization type

  // And that's it! Happy propagating-

  // TODO: Do we want to call nudge_if_singular here?

  return;

}


//////
// METHOD:  Phonon :: Velocity()
//
//   Returns seismic velocity at current location of phonon.  Depends
//   on mpCell being set.  Barfs unkindly (probably a segfault) if it
//   isn't (e.g., for "relative" phonons used in scattering
//   computations, but this is not the intended use case for this
//   method).
//
Real Phonon::Velocity() const {
  return mpCell->GetVelocAtPoint(mLoc, mType);
}


//////
// METHOD:  Phonon :: DirectionOfMotion()
//
//   Returns the direction of particle motion of the ray.  In the case
//   of P phonons, this is the same as the direction of propagation.
//   In the case of S phonons, it is perpendicular to propagation, and
//   the exact orientation is determined by the polarization angle.
//
// RETURNS:  
//
//   Returns an R3::XYZ object coding the direction of motion.  Return
//   value is guaranteed to be "normal", ie, to have a vector
//   magnitude of 1.0.
//
R3::XYZ Phonon::DirectionOfMotion() const {

  if (mType == RAY_P) {
    return mDir;
  } 
  else {
    R3::OrthoAxes OA(mDir.Theta(), mDir.Phi(), mPol);
    return OA.S1();
  }

}


//////
// METHOD:     Phonon::Refract()
// CALLED-BY:  Phonon::Propagate()
//
//   This is a dispatcher function. It's purpose is to determine what
//   type of refraction handling is required by the situtaion, and to
//   call the appropriate refraction helper functions to handle it.
//
//   On the macro level, this function handles the transmission of
//   phonons from one cell into a neighbor cell across an interface.
//
void Phonon::Refract(const CellFace * pFace) {

  Real vel_eps = 0.00001;   // Threshold below which we don't bother
                            // with Snell's law. (Basically, this
                            // discriminates between stairstep and
                            // continuous (i.e. gradient) velocity
                            // models.)
                            //
                            //  Number value (1e-5) was chosen to be
                            //  smaller than the fractional change due
                            //  to Earth-flattening transformation
                            //  over 100m (0.1km) at Earth's surface.
                            //  Smaller values could be appropriate
                            //  for some circumstances, but run the
                            //  risk of false comparissons due to
                            //  slight off-interface errors.
                            //

  if (pFace->GridDiscontinuity()) {  // If the grid asked for a 1st-
    Refraction_FullRT(pFace);        // order discontinuity here, then
  }                                  // do full R/T
  else {
    if (pFace->VelocityJump(mLoc) > vel_eps) {  // Else ray-bend only,
      Refraction_Bend(pFace); }
    else {                                      // Or do just a direct
      Refraction_Continuous(pFace);             // handover (for the 
    }                                           // continous model case.)
  } //
  //

}


//////
// METHOD:     Phonon::Refraction_Continuous()
// CALLED BY:  Phonon::Refract()
//
//   Handles the transfer of a phonon from its current cell into the
//   cell that adjoins it through pFace.
//
//   This version is for the case when velocity is continuous across
//   the interface.  For cells that model velocity gradients, this is
//   usually the case (baring special discontinuity interfaces such as
//   a Moho layer).
//
//   The phonon is basically handed over to the adjoinin cell without
//   much additional fanfare.
//
void Phonon::Refraction_Continuous(const CellFace * pFace) {

  InsertInto( &(pFace->OtherCell()) );

}


//////
// METHOD:     Phonon::Refraction_Bend()
// CALLED BY:  Phonon::Refract()
//
//   Handles the transfer of a phonon from its current cell into the
//   cell that adjoins it via an identified exit CellFace.
//
//   This version assumes the siesmic velocities are discontinuous
//   across the interface, and calculates ray-bending accordingly.
//   Useful for models based on uniform-velocity cells, whose velocity
//   models consequently become "stair-step" functions.  Falls back to
//   hard reflect if supercritical incidence (total-internal
//   reflection).  In this latter case, we don't change the
//   bookkeeping of the phonons "current cell", and the result is as
//   if the phonon "transferred" back into the same cell.
//
//   This version does NOT do P/S conversions.  Thus this function
//   should not be used to handle grid-indicated discontinuity faces.
//   Another mechanism is provided for that.
//
//   Keeps S-polarization "unchanged" throughout, which is to say that
//   the SH/SV ratio is unchanged, as measured relative to the
//   plane-of-incidence.  (This does NOT mean that the mPol member
//   will not change.  It probably will change, since it is measured
//   w.r.t a global coordinate system, not the plane of incidence.)
//
//   NOTE: because of the way exit-face intersection is computed, we
//   can count on pFace->Normal().Dot(mDir) to always be
//   positive. (ie., we are always traveling "with", not "against",
//   the exit-face normal.)
//
void Phonon::Refraction_Bend(const CellFace * pFace) {

  R3::XYZ fnorm = pFace->Normal(mLoc);
  // Face-Norm: Unit vec normal to exiting CellFace

  R3::XYZ fpara = fnorm.GetInPlaneUnitPerpendicular(mDir);
  R3::XYZ fparash = fnorm.Cross(fpara);
  // Face-Parallel and Face-Parallel-SH: both these unit vectors are
  // mutually perpendicular and are in the plane of the exit face.
  // fpara is also co-planar with the incidence plane defined by mDir
  // and fnorm.  fparash is a useful base vector for quantifying the
  // SH/SV fraction (relative to the CellFace) of the incoming ray.

  Real veli = pFace->Cell().GetVelocAtPoint(mLoc, mType);
  Real velo = pFace->OtherCell().GetVelocAtPoint(mLoc, mType);
  // Get velocities on the incoming and outgoing sides of the face

  // Get sines and cosines of the incoming and outgoing ray
  // directions.  fnorm and fpara form the basis set for this
  // computation.

  Real sini = fpara.Dot(mDir);    // sine of incoming angle
  Real sino = (velo/veli)*sini;   // The Law of Snell (TM)

  bool transfer;    // True if we transmit, false if we reflect
  Real coso;        // cosine of outgoing angle

  // Determine output direction:

  if (sino >= 1.0) {              // Total internal reflection,
    transfer = false;             // transmission blocked
    sino = sini;
    coso = -1.0 * fnorm.Dot(mDir);    // Negative-definite
  }
  else {                          // Transmission occurs
    transfer = true;
    coso = sqrt(1.0 - (sino*sino));   // Positive-definite
  }

  R3::XYZ outdir = fpara.ScaledBy(sino) + fnorm.ScaledBy(coso);

  // If we are an S-wave, then we still need to resolve 
  // polarization:

  R3::XYZ pdomi;    // Particle-Direction-of-Motion, Incoming
  R3::XYZ svbasei;  // Unit vec for SV on incoming ray
  R3::XYZ svbaseo;  // Unit vec for SV outgoing
  Real svcomi;      // S-Vertical Component, Incoming
  Real shcomi;      // S-Horizontal Component, Incoming
                    //   (NOTE: SV/SH components don't actually change
                    //   durring reflection or refraction, just the
                    //   basis vectors change, and then only the SV
                    //   base changes.  The SH basis vector (fparash)
                    //   stays the same.)
  R3::XYZ pdomo;    // Particle-Direction-of-Motion, Outgoing
  Real polout=0;    // Outgoing polarization, default to zero
                    // (arbitrary, but a good in-bounds value in case
                    // of P-wave, which doesn't need a polarization)

  if (mType != RAY_P) {   // If not a P-wave:

    // Get output polarization direction:

    pdomi = DirectionOfMotion();
    svbasei = fparash.Cross(mDir);
    svbaseo = fparash.Cross(outdir);
    shcomi = pdomi.Dot(fparash);
    svcomi = pdomi.Dot(svbasei);
    pdomo = fparash.ScaledBy(shcomi) + svbaseo.ScaledBy(svcomi);
    
    // Compute polarization angle in global coordinate system:

    Real pol_v = pdomo.Dot(outdir.ThetaHat());
    Real pol_h = pdomo.Dot(outdir.PhiHat());
    polout = atan2(pol_h, pol_v);

  } // Gets us a polarization value for S-waves
  ///  (Or leaves default value unchanges for P-waves)

  // Update Phonon motion attributes, based on above:

  mDir.Set(outdir.Theta(), outdir.Phi());
  mPol = polout;

  // Update Linkage, if "transfer" occured:

  if (transfer==true) {
    InsertInto( &(pFace->OtherCell()) );
  } // Otherwise we just stay put in the cell we're 
  ///  already in.

  // Done. It was my pleasure to direct your phonon 
  // today.

}


//////
// METHOD:     Phonon::Refraction_FullRT()
// CALLED BY:  Phonon::Refract()
//
//   Handles the transfer of a phonon from its current cell into the
//   cell that adjoins it through pFace.
//
//   This version is for the specially marked interfaces that are
//   intended to represent first-order discontinuities (velocity
//   jumps) across an interface, e.g. a Moho discontinuity.  Full
//   reflection/transmission and P/S conversion probabilities are
//   handled.
//
//   Contrast this with Refraction_Bend(), which does handle
//   discontinuous velocities across an interface, but treats only the
//   simple bending of rays by Snells law. This latter case is
//   appropriate for so-called stair-step models in which velocities
//   are stepwise-constant functions.  In order to get FullRT
//   treatment on a particular interface, the user must request this
//   via a double-valued grid layer.
//
void Phonon::Refraction_FullRT(const CellFace * pFace) {

  RTCoef rt = pFace->GetRTBasis(mLoc, mDir);  
                      // Gets basis vectors for scattered
                      // rays and velocities, densities

  raytype intype = RAY_P;
  if (mType == RAY_S) {
    intype = rt.ChooseSPolType(DirectionOfMotion());
  }//
  //    We are now either RAY_P, RAY_SH, or RAY_SV.
  //

  rt.GetCoefs(intype);  // Generate output sines and probabilities for
                        // all cases allowed by the incident raytype.
                        // (Results stored internally to rt)

  rt.Choose();          // Pick an outcome from the allowed set.

  // ::::
  // :: Update internal state for new phonon trajectory:
  // :

  raytype outtype = rt.GetChosenRaytype();      // Returns RAY_P or RAY_S
  bool   transmit = rt.DidRayTransmit();
  R3::XYZ  outdir = rt.GetChosenRayDirection();

  mType = outtype;                              // Set Ray Type
  mDir.Set(outdir.Theta(), outdir.Phi());       // Set Ray Direction 

  if (mType == RAY_S) {                         // Set Polarization
    R3::XYZ pdomo                               // (Direction of Particle
          = rt.GetChosenParticleDOM();          // Motion)
    Real pol_v = pdomo.Dot(mDir.ThetaHat());
    Real pol_h = pdomo.Dot(mDir.PhiHat());
    mPol = atan2(pol_h, pol_v);
  }

  if (transmit==true) {                         // And update linkage if
    InsertInto( &(pFace->OtherCell()) );        // we transmitted
  }

  // :
  // :: ...And, that should do it.
  // :                        Done.
  //       In the beginning, there was the Word,
  //         And the Word was God.
}


//////
// METHOD:  Phonon :: InsertInto( MediumCell * )
//
// WHAT-IT-DOES:
//
//   This methods makes the necessary internal changes to "link" a
//   phonon to the MediumCell in which it is to reside for the next
//   few rounds of simulation.  At a basic level, this means setting
//   the internal mpCell and mpScat members, but, depending on
//   situation, this could also mean making trajectory or other
//   changes.
//
//   UPDATE: Trajectory and other such "motion" parameters are handled
//   elsewhere.  This function is now only about linkage.
//
void Phonon::InsertInto(MediumCell * pCell) {

  mpCell = pCell;           // "Link" us to the MediumCell

  mpScat = pCell->GetActiveScatterer();   // And "link" us to the
                                          // Scatterer object active
                                          // within that cell.

}


//////
// METHOD:  Phonon :: Propagate()
//
// WHAT-IT-DOES:
//
//   Takes a phonon which is assumed to exist at somesuch location in
//   somesuch MediumCell and "propagates" that phonon along its
//   directed path, processing all possible occurances (like
//   scattering and reflections) along the way.  Typically, this
//   function would be called in a loop soon after a "new" phonon is
//   generated by a seismic SourceEvent object.
//
// HOW-IT-WORKS:
//
//   (Mostly documented in comments inside the function)
//
// INPUTS:
//   (none)
//
// OUTPUTS:
//   (none)
//
// SIDE-AFFECTS:
//
//   Upon comletion, phonon is in a "dead" state, meaning it has
//   either left the model (and no longer exists inside any valid
//   MediaCell), or because its mTimeAlive has exceeded cm_ttl
//   (meaning phonon has propagated longer than the time-to-live limit
//   on propagation sim-time and has therefore been discarded).
//
//   During propagation, various phonon events, such as scattering
//   incidents, reflections, intersections with collection faces,
//   etc., are reported to the outside world through method calls to
//   the global DataOut object.
//
void Phonon::Propagate() {

  while (true) {    // ===========================
                    //   PROPAGATION INNER LOOP:

    // ::::::
    // :: Test for TIME-OUT and VALIDITY:
    // :

    if (mTimeAlive > cm_ttl) {
      dataout.ReportPhononTimeout(*this);
      break;
    }

    if ((mMoveCount % 128)==127) { // Various VALIDITY checks:
      if (std::isnan(mPathLength)) {
        dataout.ReportInvalidPhonon(*this, DataReporter::INV_PATH_NAN);
        break;
      }
      if (std::isnan(mTimeAlive)) {
        dataout.ReportInvalidPhonon(*this, DataReporter::INV_TIME_NAN);
        break;
      }
      if (mPathLength < 0) {
        dataout.ReportInvalidPhonon(*this, DataReporter::INV_PATH_NEGATIVE);
        break;
      }
      if ((mTimeAlive < 0) || (mRecentTravelTime < 0)) {
        dataout.ReportInvalidPhonon(*this, DataReporter::INV_TIME_NEGATIVE);
        break;
      }
      if (mRecentTravelTime == 0) {
        dataout.ReportInvalidPhonon(*this, DataReporter::INV_STUCK);
        break;
      }
      if (mRecentTravelTime < cm_slow_concern) {
        dataout.ReportInvalidPhonon(*this, DataReporter::INV_SLOW);
        break;
      }
      if (mMoveCount > cm_loop_concern) {
        dataout.ReportInvalidPhonon(*this, DataReporter::INV_LOOP_EXCEED);
        break;
      }
      mRecentTravelTime = 0;
    }

    // ::::::
    // :: Test for SCATTERRING event:
    // :

    TravelRec travel =                  // Travel record to cell
      mpCell->GetPathToBoundary(mType,  // boundary. (Gets overridden
                                mLoc,   // below if we scatter)
                                mDir); 

    if (travel.PathLength == 1.0/0.0){
      dataout.ReportPhononTimeout(*this);
      break;
    }


    Real scatlen = mpScat->GetRandomPathLength(mType);
    Real edgelen = travel.PathLength;


    if (scatlen < edgelen) {  // true if we scattered,
                              // false if we hit cell boundary

      travel = mpCell->AdvanceLength(mType, scatlen, mLoc, mDir);
      this->Move(travel);     // Move us to the scatterer

      Phonon rph = mpScat->GetRandomScatteredRelativePhonon(mType);
      this->Transform(rph);   // Re-orient and re-raytype us according
                              // the results of the scattering event
                              // (coded in 'rph').

      dataout.ReportScatterEvent(*this);
     
      continue;

    } // If we get here, then we did NOT scatter. (We are
      // at a cell boundary)

    this->Move(travel);       // Move us to the cell boundary
   
    // ::::::
    // :: Test for COLLECTION:
    // :

    if (travel.pFace->IsCollectionFace()) {
      dataout.ReportPhononCollected(*this);
    } // Our interaction with a collection-surface has been
      // reported. Phonon does not (necessarily) stop here. Continue
      // propagation handling below:


    // ::::::
    // :: Test for REFLECTION:
    // :

    if (travel.pFace->IsReflectionFace()) { // Use R/T coefficient treatment
                                            // to handle reflection with P/S
      Refraction_FullRT(travel.pFace);      // conversion. Note: Unexpected
      dataout.ReportReflection(*this);      // behavior may result if
      continue;                             // reflection face is not a
                                            // free-surface face.
    } 
      // If we get here, then we did NOT reflect.  Next check for
      // transmission into adjoining cell:


    // ::::::
    // :: Test for NEIGHBOR:
    // :

    if (travel.pFace->HasNeighbor()) {
      MediumCell * oldcell = mpCell;  // Remember where we were
      Refract(travel.pFace);          // Reflect or Refract into next cell 
      if (mpCell == oldcell) {            // If mpCell hasn't changed
        dataout.ReportReflection(*this);  // then we reflected
      } else {
        dataout.ReportCellToCell(*this);  // else we transmitted
      }

      continue;

    } // If we get here, then phonon left cell, but not into a
      // neighbor cell (because there isn't one). We have transmitted
      // into "unmodelled space."  Phonon is, essentially, lost.


    // ::::::
    // :: Phonon LOST:
    // :

    dataout.ReportLostPhonon(*this);
    break;
    

  }   // END PROPAGATION INNER LOOP

////
}// END: Phonon::Propagate()
//
