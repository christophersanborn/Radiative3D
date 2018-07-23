/// @file raytype.hpp
///
/// Defines a simple enum to give symbolic representation to ray
/// polarization types (P, SH, SV).  This used to be defined in
/// phonon.hpp, but including the whole phonon class in headers that
/// only needed the raytype enum became cumbersome.
///
#ifndef RAYTYPE_H_
#define RAYTYPE_H_
//

//////
// ENUM: raytype
///
///  Enumerates polarization modes for Phonon objects and other class
///  methods where wave mode or "ray type" matters.
///
///  On the most basic level, we distinguish between P (longitudinal)
///  and S (transverse) polarization. Most methods only make this
///  distinction. Some methods distinguish further between horizontal
///  and vertical S polarization, yielding the RAY_SH and RAY_SV
///  indices.  This is mostly the case in Phonon generation from a
///  source event or scatter event.  The Phonon object itself tracks S
///  polarization as a continuous value and thus only retains RAY_P or
///  RAY_S status.
///
///  At present time, more complex polarization modes (quasi-P,
///  quasi-S, etc.) are not tracked or represented.
///
enum raytype {RAY_P = 0,   ///< Longitudinal or 'P' polarization
              RAY_S = 1,   ///< Used when only differentiating btw/ P and S
              RAY_SH = 1,  ///< Used when all three types matter
              RAY_SV = 2,  ///< ''
              RAY_NUMTYPES = 3,         ///< P, SH, SV
              RAY_NUMBASICTYPES = 2,    ///< P, S
              RAY_NT = 3,               ///< Alias for RAY_NUMTYPES
              RAY_NBT = 2,              ///< Alias for RAY_NUMBASICTYPES
              RAY_NA = 0}; ///< Used when raytype "not applicable"


///
#endif //#ifndef RAYTYPE_H_
//
