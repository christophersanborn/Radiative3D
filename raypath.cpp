// raypath.cpp
//
#include "raypath.hpp"

Real RayArcAttributes::AngleOffsetFromBottom(const R3::XYZ & loc) const {
  const R3::XYZ centerToLoc = Center.VectorTo(loc);
  const Real loctanx = u3.Dot(centerToLoc);
  const Real loctany = u1.Dot(centerToLoc);
  return atan2(loctany,loctanx);
}

std::string RayArcAttributes::str() const {
  std::ostringstream s;
  s << "{RayArc:"<<" {Radius: "<<Radius<<" @ Center: "<<Center.str()<<""
    <<"  u1: "<<u1.str()<<" u3: "<<u3.str()<<"}}";
  return s.str();
}
