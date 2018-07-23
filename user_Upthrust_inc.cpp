// user_Upthrust_inc.cpp
//
//   Defines a specific model.
//   #include this file in user.cpp
//

//////
// Crustal Upthrust Model
//
// Model a collisional upthrust into the bulge region.
//
// Based on same fan-shaped azimuthally-symmetric design as the NSCP
// CrustPinch model.
//
// Uses RAE (Range, Azimuth, Elevation) coordinate scheme. For version
// built with XYZ coordinates, revert to revision 912.
//
// First, a namespace for helper functions:
//
namespace CrustUpthrustHelper {  

  using Elastic::Q;
  using Elastic::HetSpec;

  EarthCoords::Generic GetFanRAZ(Real, Index, Index);
  GridNode GetPlumbNode(const Grid &, Real, Real, Real, Index, Index);
  void ProcessModelArgs(const std::vector<Real> &,
                        HetSpec &, HetSpec &, HetSpec &, HetSpec &, HetSpec &,
                        Q &, Q &, Q &, Q &, Q &, Real &, Real &, Real &,
                        Real &, Real &, Real &, Real &);

  // The following define the "fan" structure:
  const Count nR = 15;    // How many range points, incl negative ranges
  const Count nRneg = 2;  // How many ranges are behind origin
  const Count nCR = 1;    // How many on eiter side of orig are "close range"
  const Index iAR = 5;    // First index of "active range"
  const Count nAR = 6;    // How many indices enclose "active region"
  const Count nAzis = 6;  // How many azimuths

  Real Ranges[] =   {-120, -60, 60, 120, 220,       // Model func may modify
                     310, 365, 418, 422, 475, 530,  // the "active" range
                     650, 770, 890, 1020};          //
  const Real Azis[] = 
    {45.0, 56.25, 78.75, 101.25, 123.75, 135.0};    // Azis in degrees
  const Real AzisCR[] = 
    { 5.0, 39.00, 73.00, 107.00, 141.00, 175.0};    // Close-range spread

}
//
void CrustUpthrustWCG(Grid & gr, const std::vector<Real> & args) {

  using Elastic::Velocity;
  using Elastic::VpVs;
  using Elastic::Q;
  using Elastic::QmQk;
  using Elastic::HetSpec;
  using Elastic::HSneak;
  using CrustUpthrustHelper::ProcessModelArgs;
  using CrustUpthrustHelper::GetPlumbNode;
  using CrustUpthrustHelper::GetFanRAZ;

  Real sc_nu  = 0.8;    // Default scat params -
  Real sc_eps = 0.01;   // we'll scan args for user values
  Real sc_a   = 4.00;   //
  Real sc_k = 0.8;      //
  Real q = 200;         //

  HetSpec HSSe = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Sedi
  HetSpec HSCr = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Crust
  HetSpec HSAR = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Crust Active Region
  HetSpec HSMo = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Moho
  HetSpec HSMa = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Mantle

  Q QSe = QmQk(q);          //
  Q QCr = QmQk(q);          //
  Q QAR = QmQk(q);          //
  Q QMo = QmQk(q);          //
  Q QMa = QmQk(q);          //

  Real SediThick  = 2.0;    // Thickness of Sediments layer
  Real CrustThick = 30.0;   // Thickness of (unpinched) Crust layer
  Real MohoThick  = 10.0;   // Thickness of Moho layer
                            //
  Real UpthrL = -10;
  Real UpthrR = 0;
  Real IGamma = 1.0;
  Real Shear = 0;

  ProcessModelArgs(args,
                   HSSe, HSCr, HSAR, HSMo, HSMa,
                   QSe,  QCr,  QAR,  QMo,  QMa,
                   SediThick, CrustThick, MohoThick,
                   UpthrL, UpthrR, IGamma, Shear);

  using CrustUpthrustHelper::nR;    // How many range points
  using CrustUpthrustHelper::iAR;   // First index of "active" range
  using CrustUpthrustHelper::nAR;   // How many, incl endpoints
  using CrustUpthrustHelper::nAzis; // How many azimuths

  Count nZ = 8;     // How many z-depths in basic depth model,
  Count nZin = 10;  // Upthrust will require 2 add'nl interpolated depths.

  Real ZBase[]  = {   0.0,              // Basic Z-depth model
                    -SediThick,
                    -(SediThick+CrustThick), 
                    -(SediThick+CrustThick+MohoThick),
                      -80,    // Very subtle discontinuity here.
                      -120,
                      -210,
                      -360
  };

  // First we create a PlumbLine grid, to capture the basic depth
  // structure. This will be used as a template when building the real
  // grid.

  Grid grPL;
  Grid grPLAR = grPL;   // Active region PL, may have enhanced scattering

  grPL.SetSize(1,1,nZ);
  grPL.SetIndexBase(0);
  grPL.SetMapping(Grid::GC_RAE, Grid::GC_CURVED);
  for (Index iz=0; iz<nZ; iz++) {
    grPL.WNode(0,0,iz).SetLocation(0,0,ZBase[iz]);
  }
  grPL.WNode(0,0,0).SetAttributes(VpVs(4.50, 2.60), 2.20, QSe, HSSe);
  // ### Sediments Layer ### //
  grPL.WNode(0,0,1).SetAttributes(VpVs(4.52, 2.61), 2.21, QSe, HSSe);
  grPL.WNode(0,0,1).SetAttributes(VpVs(6.20, 3.58), 2.80, QCr, HSCr);
  // ### Crust Layer ####### //
  // ### Crust Layer ####### //
  // ### Crust Layer ####### //
  grPL.WNode(0,0,2).SetAttributes(VpVs(6.24, 3.60), 2.82, QCr, HSCr);
  grPL.WNode(0,0,2).SetAttributes(VpVs(7.70, 4.44), 3.39, QMo, HSMo);
  // ### Moho Transition ### //
  grPL.WNode(0,0,3).SetAttributes(VpVs(8.00, 4.46), 3.40, QMa, HSMa);
  // # Mantle Layer (Moho to -80 km)
  // # Mantle Layer (Moho to -80 km)
  grPL.WNode(0,0,4).SetAttributes(VpVs(8.040, 4.48), 3.50, QMa, HSMa);
  grPL.WNode(0,0,4).SetAttributes(VpVs(8.045, 4.49), 3.50, QMa, HSMa);
  // # Mantle Layer (-80 to -120 km)
  grPL.WNode(0,0,5).SetAttributes(VpVs(8.051, 4.50), 3.43, QMa, HSMa);
  // # Mantle Layer (-120 to -210 km)
  grPL.WNode(0,0,6).SetAttributes(VpVs(8.301, 4.52), 3.32, QMa, HSMa);
  // # Mantle Layer (-210 to -360 km)
  grPL.WNode(0,0,7).SetAttributes(VpVs(8.848, 4.78), 3.46, QMa, HSMa);

  // Active region, PLAR: Copy and modify PL; (Use QAR and HSAR
  grPLAR = grPL;                           //  instead of QCr and HSCr)
  // ... Sediments Layer ### //
  grPLAR.WNode(0,0,1).ClearAttributes();
  grPLAR.WNode(0,0,1).SetAttributes(VpVs(4.52, 2.61), 2.21, QSe, HSSe);
  grPLAR.WNode(0,0,1).SetAttributes(VpVs(6.20, 3.58), 2.80, QAR, HSAR);
  // ### Crust Layer ####### //
  // ### Crust Layer ####### //
  // ### Crust Layer ####### //
  grPLAR.WNode(0,0,2).ClearAttributes();
  grPLAR.WNode(0,0,2).SetAttributes(VpVs(6.24, 3.60), 2.82, QAR, HSAR);
  grPLAR.WNode(0,0,2).SetAttributes(VpVs(7.70, 4.44), 3.39, QMo, HSMo);
  // ### Moho Transition ... //
  // ... Remaining layers are unmodified //

  //
  // Now we build the real grid, 
  //   using the plumblines as a template:
  //

  gr.SetSize(nR,nAzis,nZin);    // Sets index bounds, incl 2 extra z levels
  gr.SetIndexBase(0);           // When addressing nodes, use base 0
  gr.SetMapping(Grid::GC_RAE, Grid::GC_CURVED);

  for (Index ir = 0; ir < nR; ir++) {
    Real UpthrLS = UpthrL;        // "Selected" values; may be flip-flopped
    Real UpthrRS = UpthrR;        // is thrust direction is reversed. See
    Real ShearS = Shear;          // conditional below.
    Index m; bool bAR;            //
    m = (ir > iAR) ? ir-iAR : 0;  // Index into active region, but
    m = (m < nAR) ? m : nAR-1;    // truncate at iAR_max
    bAR = (ir>=iAR && ir<(iAR+nAR-1)); // Excludes right-most enclosing index;
                                       // HetSpec applies block ahead of index.
    Grid & grSelect = bAR ? grPLAR : grPL; // Reference either the base
                                           // plumbline or the AR one.
    if (UpthrR < UpthrL) {    // If thrust direction is reversed
      m = (nAR-1)-m;          // then flip m value
      UpthrLS = UpthrR;       // ... and Upthr values
      UpthrRS = UpthrL;
      ShearS = -Shear;
    } else if (UpthrR==UpthrL) {  // TODO: breaks for equivalencies not 
      m=0;                        // not equal to zero. Prolly need to make
      ShearS=0;                   // m a /\ triangle pattern.
    }
    for (Index iz = 0; iz < nZin; iz++) {  // Loop over depth index:
      GridNode N = GetPlumbNode(grSelect, UpthrLS, UpthrRS, IGamma, m, iz);
      for (Index iaz = 0; iaz < nAzis; iaz++) { // Loop over azi index:

        Real Z = N.GetRawLoc().x3();
        EarthCoords::Generic RAZ = GetFanRAZ(Z,ir,iaz);

        gr.WNode(ir,iaz,iz).SetLocation(RAZ);
        gr.WNode(ir,iaz,iz).SetAttributes(N.Data(GridNode::GN_ABOVE));
        if (N.IsDiscontinuous()) {
          gr.WNode(ir,iaz,iz).SetAttributes(N.Data(GridNode::GN_BELOW));
        }
        if ((m==2||m==3)&&(iz==2||iz==3)) {
          gr.WNode(ir,iaz,iz).AdjustLocation(-ShearS,0,0);
        } else if ((m==2||m==3)&&(iz==4||iz==5)) {
          gr.WNode(ir,iaz,iz).AdjustLocation(+ShearS,0,0);
        }
      }
    }
  }// End loop over ir, iz, iaz
  //

}
//
// Continue helper namespace: Define functions
//
namespace CrustUpthrustHelper {

//////
// FUNCTION:  GetFanRAZ()
//
//  Get RAZ coords for a fan-shaped model given a range index and azi
//  index and an already determined z coord.
//
//  Uses namescape level globals Azis, AzisCR, nAzis, Ranges, nR,
//  nRneg, and nCR, all defined at the top of this file.
//
EarthCoords::Generic GetFanRAZ(Real z, Index ir, Index iaz) {

  Real azi  = Azis[iaz];  // Azimuth value (default)

  // TODO: Bounds checking on ir, iaz
  // TODO: Sanity check that nRneg <= nR, or better yet determine
  // nRneg by analyzing Ranges[] list.
  // TODO: Sanity check that range value is not zero

  RelIndex m = ir - nRneg; 

  if (m < 0) {  // then R is negative; use reverse sequence for azi's
    if (m+(RelIndex)nCR >= 0) { // then R is negative AND close range
      azi = AzisCR[nAzis-iaz-1];
    } else {          // then R is negative and regular range
      azi = Azis[nAzis-iaz-1];
    }
  } else {      // else R is positive
    if (m < (RelIndex)nCR) {  // and R is close range
      azi = AzisCR[iaz];
    } // else keep default azi
  } //

  return EarthCoords::Generic(Ranges[ir], azi, z);

}//
//

//////
// FUNCTION:  GetPlumbNode()
//
//  Returns appropriate node from a PlumbLine (1-D in depth) model,
//  after inserting two interpolated nodes either in the crust or in
//  the mantle, needed to handle a lateral discontinuity in crust
//  structure that simulates a subducting or upthrusting slab.
//
//  Index m has six allowed values:
//
//   m    Action:
//  ===  =========
//   0    Insert nodes between n=1 and n=2; no Z adjust at moho
//   1     ''     Partial Z adjust at moho, using UpthrL
//   2     ''     Full Z adjust at moho, using UpthrL
//   3    Insert nodes between n=3 and n=4; Full Z adjust at moho (UpthrR)
//   4     ''     Partial Z adjust at moho, using UpthrR
//   5     ''     No Z adjust at moho
//
//  Index n is the z-index.  We use index k to represent the correct
//  corresponding index into grPL, except for two n values for which
//  interpolated GridNodes will need to be constructed.
//
GridNode GetPlumbNode(const Grid & grPL, Real UpthrL, Real UpthrR, 
                      Real IGamma, const Index m, const Index n) {

  // TODO: Sanity check that m not > 5 or < 0

  // Simplest case: Above and below the modification layers
  if ((n<2)||(n>5)) {           // Covers n<2 U n>5; for all m
    Index k = (n>5) ? n-2 : n;  //
    return grPL.Node(0,0,k);
  } // If we get here, then n is in [2,5]
   //

  Real ZAdjFrac = (m<3) ?
                  ((Real)m)/2.0
                : ((Real)(5-m))/2.0;
  ZAdjFrac = std::pow(ZAdjFrac,IGamma); //
  Real ZAdj = (m<3) ?                   // Amount to vertically shift Moho
              UpthrL*ZAdjFrac           // layer (k indices 2 and 3).
            : UpthrR*ZAdjFrac;          //


  // Slightly more complicated:
  if ((n<4)!=(m<3)) {           // XOR: If we make a matrix of m,n values,
                                // this selects the UR and LL quadrants.
                                // For these values, no grid data or Z
                                // value interpolation is needed.
    Index k = (n<4) ? n : n-2;
    GridNode retval = grPL.Node(0,0,k);
    retval.AdjustLocation(0,0,ZAdj);
    return retval;

  } // If we get here, then we're in the UL or LR quadrant
   //


  // And here it gets tricky: (Need to interpolate)

  Real Thick = grPL.Node(0,0,2).GetRawLoc().x3()
             - grPL.Node(0,0,3).GetRawLoc().x3();
  Real ZAbv = (m<3) ? grPL.Node(0,0,1).GetRawLoc().x3()
                    : grPL.Node(0,0,3).GetRawLoc().x3() + ZAdj;
  Real ZBlw = (m<3) ? grPL.Node(0,0,2).GetRawLoc().x3() + ZAdj
                    : grPL.Node(0,0,4).GetRawLoc().x3();
  Real ZMid = (ZAbv+ZBlw)/2;
  Real ZInt;

  if (m==2) {
    ZInt = grPL.Node(0,0,n).GetRawLoc().x3() + UpthrR;
  } else if (m==3) {
    ZInt = grPL.Node(0,0,n-2).GetRawLoc().x3() + UpthrL;
  } else {
    Real nudge = (n==2||n==4) ? 0.5*Thick : -0.5*Thick;
    ZInt = ZMid+nudge;
  } // OK, now I've got the Z value in ZInt.

  Real VpAbv = (m<3) ? grPL.Node(0,0,1).Data(GridNode::GN_BELOW).Vp()
                     : grPL.Node(0,0,3).Data(GridNode::GN_BELOW).Vp();
  Real VpBlw = (m<3) ? grPL.Node(0,0,2).Data(GridNode::GN_ABOVE).Vp()
                     : grPL.Node(0,0,4).Data(GridNode::GN_ABOVE).Vp();
  Real VsAbv = (m<3) ? grPL.Node(0,0,1).Data(GridNode::GN_BELOW).Vs()
                     : grPL.Node(0,0,3).Data(GridNode::GN_BELOW).Vs();
  Real VsBlw = (m<3) ? grPL.Node(0,0,2).Data(GridNode::GN_ABOVE).Vs()
                     : grPL.Node(0,0,4).Data(GridNode::GN_ABOVE).Vs();
  Real RhoAbv = (m<3) ? grPL.Node(0,0,1).Data(GridNode::GN_BELOW).Rho()
                      : grPL.Node(0,0,3).Data(GridNode::GN_BELOW).Rho();
  Real RhoBlw = (m<3) ? grPL.Node(0,0,2).Data(GridNode::GN_ABOVE).Rho()
                      : grPL.Node(0,0,4).Data(GridNode::GN_ABOVE).Rho();
  Elastic::Q Q = (m<3) ? grPL.Node(0,0,1).Data(GridNode::GN_BELOW).getQ()
                       : grPL.Node(0,0,3).Data(GridNode::GN_BELOW).getQ();
  Elastic::HetSpec HS = (m<3) ? grPL.Node(0,0,1).Data(GridNode::GN_BELOW).getHS()
                              : grPL.Node(0,0,3).Data(GridNode::GN_BELOW).getHS();

  // Interpolate elastic values:

  Real WAbv = (ZBlw-ZInt)/(ZBlw-ZAbv);
  Real WBlw = 1-WAbv;
  Real Vp = WAbv*VpAbv+ WBlw*VpBlw;
  Real Vs = WAbv*VsAbv+ WBlw*VsBlw;
  Real Rho = WAbv*RhoAbv+ WBlw*RhoBlw;

  // Now construct and return GridNode

  GridNode retval;
  retval.SetLocation(0,0,ZInt);
  retval.SetAttributes(Elastic::VpVs(Vp,Vs),Rho,Q,HS);
  return retval;

}//
//

//////
// FUNCTION:  ProcessModelArgs()
//
void ProcessModelArgs(const std::vector<Real> & args,
                      HetSpec & HSSe, HetSpec & HSCr, HetSpec & HSPi, 
                      HetSpec & HSMo, HetSpec & HSMa,
                      Q & QSe, Q & QCr, Q & QPi, Q & QMo, Q & QMa,
                      Real & SediThick, Real & CrustThick, Real & MohoThick,
                      Real & PinchFrac, Real & DepthFrac,
                      Real & SediFrac, Real & MohoFrac) {

  using Elastic::QmQk;
  using Elastic::HSneak;

  // Check args array and modify parameters appropriately:

  switch(args.size()) {
  case 0:   // No args given, just use defaults already constructed
            //
    break;

  case 32:  // 5 x neakq + layer thicks + pinch fractions
            //
    PinchFrac = args.at(28);
    DepthFrac = args.at(29);
    SediFrac  = args.at(30);
    MohoFrac  = args.at(31);
    /* FALL THROUGH  - no break; */

  case 28:  // 5 x neakq + layer thicknesses
            //
    SediThick  = args.at(25);
    CrustThick = args.at(26);
    MohoThick  = args.at(27);
    /* FALL THROUGH  - no break; */

  case 25:  // 5 x neakq
            //
    HSSe = HSneak(args.at(0), args.at(1), args.at(2), args.at(3));
    QSe  = QmQk(args.at(4));
    HSCr = HSneak(args.at(5), args.at(6), args.at(7), args.at(8));
    QCr  = QmQk(args.at(9));
    HSPi = HSneak(args.at(10), args.at(11), args.at(12), args.at(13));
    QPi  = QmQk(args.at(14));
    HSMo = HSneak(args.at(15), args.at(16), args.at(17), args.at(18));
    QMo  = QmQk(args.at(19));
    HSMa = HSneak(args.at(20), args.at(21), args.at(22), args.at(23));
    QMa  = QmQk(args.at(24));
    break;

  default:  // Unrecognized pattern of values
            //
    std::cerr << "Error: wrong number of model args passed "
              << "to compiled-in grid-building function.\n";
    exit(1);  // TODO: Raise a meaningful exception instead
    break;

  }

}

}//
// END namespace CrustUpthrustHelper;
//
