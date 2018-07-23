// geom_base.hpp
///
/// The geom library defines namespaces containing classes and other
/// constructs related to the geometry of spaces in various dimensions.
/// This file defines:
///
/// o Namespace: Geometry
///
///    Mostly constants like pi and degrees-to-radians factors, etc.
///    Meant to be of common utility to the other two spaces.
///
///
#ifndef GEOM_BASE_H_
#define GEOM_BASE_H_
//
#include "typedefs.hpp"

namespace Geometry {
//__________________________________________________________________________
//**************************************************************************
// NAMESPACE: Geometry
// PURPOSE:
///
///   To provide some basic constants and other constructs of common
///   utility to all the OTHER geometry based namespaces.
///
//__________________________________________________________________________
//**************************************************************************

//////
// CONSTANTS:
//
const Real Pi    = 3.14159265358979323846;
const Real Pi45  = Pi * 0.25;
const Real Pi90  = Pi * 0.5;
const Real Pi180 = Pi;
const Real Pi270 = Pi * 1.5;
const Real Pi360 = Pi * 2.0; 
const Real RtoD  = 180.0/Pi;
const Real DtoR  = Pi/180.0;
//
const Real Sqrt2    = 1.41421356237309504880;
const Real Sqrt2Inv = 1.0 / Sqrt2;
const Real Sqrt3    = 1.73205080756887729353;
const Real Sqrt3Inv = 1.0 / Sqrt3;
//
/////
};// end namespace: Geometry
///

#endif //#ifndef GEOM_BASE_H_
