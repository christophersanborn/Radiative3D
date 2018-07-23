// raytype.hpp
//
// Defines a simple enum to give symbolic representation to ray
// polarization types (P, SH, SV).  This used to be defined in
// phonon.hpp, but including the whole phonon class in headers that
// only needed the raytype enum became cumbersome.
//
#ifndef RAYTYPE_H_
#define RAYTYPE_H_
//

enum raytype {RAY_P = 0,
              RAY_S = 1,   // Used when only differentiating btw/ P and S
              RAY_SH = 1,  // Used when all three types matter
              RAY_SV = 2,  // ''
              RAY_NUMTYPES = 3,         // P, SH, SV
              RAY_NUMBASICTYPES = 2,    // P, S
              RAY_NT = 3,               // Alias for RAY_NUMTYPES
              RAY_NBT = 2,              // Alias for RAY_NUMBASICTYPES
              RAY_NA = 0}; // Used when raytype "not applicable"


///
#endif //#ifndef RAYTYPE_H_
//
