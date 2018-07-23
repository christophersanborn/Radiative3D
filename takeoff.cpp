// takeoff.cpp
///@file
///
/// Implement take-off angle classes
///

#include "takeoff.hpp"

TakeoffArray::TakeoffArray(Count size) {
  mvTOA.reserve(size);
}

//////
/// 20 * 4 ^ degree
int NTesselSphere::NumAnglesFromDegree(const int degree) {
  return 20 * pow(4,(degree>0?degree:0));
}

NTesselSphere::NTesselSphere(const int degree) :
  TakeoffArray(NumAnglesFromDegree(degree))
{
  // Twelve nodes make up the vertices of an icosahedron. The nodes
  // can be aligned so that each node extends in only two of the three
  // cartesian basis directions, which forms the basis of the labeling
  // convention we use here. For details see LabNotebook entry at
  // https://rainbow.phys.uconn.edu ...
  // ... /geophysics/wiki/index.php/Radiative3D_//_2013-03_MARCH

  using S2::Node;
  const Real phi = (1.+sqrt(5.0))/2.0;  // Golden ratio when other
                                        // sidelength is 1.0

  Node NF = Node(   1,    0,  phi );
  Node NB = Node(  -1,    0,  phi );
  Node SF = Node(   1,    0, -phi );
  Node SB = Node(  -1,    0, -phi );
  Node FL = Node( phi,   -1,    0 );
  Node FR = Node( phi,    1,    0 );
  Node BL = Node(-phi,   -1,    0 );
  Node BR = Node(-phi,    1,    0 );
  Node RN = Node(   0,  phi,    1 );
  Node RS = Node(   0,  phi,   -1 );
  Node LN = Node(   0, -phi,    1 );
  Node LS = Node(   0, -phi,   -1 );

  // Twenty faces make up the icosahedron.  Here we create a 2-dimensional
  // array to represent the faces as node * tripples.

  //There are altogether 12 faces that touch the "edges" of the three
  //golden rectangles, and then 8 more that span the corners, but not
  //the edges.  We list the edge ones first, in order of N,S,F,B,L,R,
  //then the corner ones, starting with the four in the northern
  //hemisphere, and then the four in the southern.

  const int num_faces = 20;
  Node * fa [num_faces] [3] = {   // Face Array:
    {&NF, &NB, &LN},    // The two that touch the North edge
    {&NF, &NB, &RN},
    {&SF, &SB, &LS},    // The two that touch the South edge
    {&SF, &SB, &RS},
    {&FL, &FR, &NF},    // The two that touch the Front edge
    {&FL, &FR, &SF},
    {&BR, &BL, &NB},    // The two that touch the Back edge
    {&BR, &BL, &SB},
    {&LN, &LS, &FL},    // The two that touch the Left edge
    {&LN, &LS, &BL},
    {&RN, &RS, &FR},    // The two that touch the Right edge
    {&RN, &RS, &BR},
    {&NF, &LN, &FL},    // The four that touch corners in North hemisphere
    {&NF, &RN, &FR},
    {&NB, &LN, &BL},
    {&NB, &RN, &BR},
    {&SF, &LS, &FL},    // The four that touch corners in South hemisphere
    {&SF, &RS, &FR},
    {&SB, &LS, &BL},
    {&SB, &RS, &BR}

  };

  // Now create a temporary fully-subdivided triangle for each face,
  // and harvest the centerpoints arrays:

  for (int ctr=0; ctr<num_faces; ctr++) {
    // For each icosahedron face:
    S2::Triangle tri(*fa[ctr][0], *fa[ctr][1], *fa[ctr][2], degree);
    for (int subctr=0; subctr<tri.NumCenterAngles(); ctr++) {
      // For each minimal subtriangle, get direction and weight: 
      S2::ThetaPhi toadir = tri.CenterAngle(subctr);
      Real weight = solid_to_weight(tri.AngleArea(subctr));
      mvTOA.push_back(TakeoffAngle(toadir,weight));
    }
  }

}

//
