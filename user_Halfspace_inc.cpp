// user_Halfspace_inc.cpp
//
//   Defines a specific model.
//   #include this file in user.cpp
//

//////
// HalfspaceCylinder()
//
//   Models a Halfspace model using Cylinder grid metaphor.  Actually models
//   two layers, for *some* simple flexibility, but if both layers set equal
//   it's a proper half-space.  Might extend it with some options for
//   inclining interface between top and bottom layer, and for signally R/T
//   handling at layer interface.
//
//   This grid can be used for, among other things, seismic albedo testing
//   and MLTW method testing.
//
//   Scattering parameters, Q values, velocities, and layer bottoms are
//   specified in the 'args' argument and are provided by the user by
//   supplying the --model-compiled-args switch.
//
//   --model-compiled-args=nu1,eps1,a1,k1,Qs1,
//                         nu2,eps2,a2,k2,Qs2,
//                         Vp1,Vs1,rho1,zbottom1,
//                         Vp2,Vs2,rho2,zbottom2
//
void HalfspaceCylinder(Grid & gr, const std::vector<Real> & args) {

  using Elastic::Velocity;
  using Elastic::VpVs;
  using Elastic::Q;
  using Elastic::QmQk;
  using Elastic::HetSpec;
  using Elastic::HSneak;

  // Grid Dimension:

  gr.SetSize(3,1,3);    // 3 Grid sheets for TWO model layers
  gr.SetIndexBase(0);   //


  // Site locations:

  Real EpiX =  0.0;
  Real EpiY =  0.0;
  Real EASX = 100.0;
  Real EASY =   0.0;
  Real NORX =   0.0;
  Real NORY =  100.0;

  // Scattering Params, defaults and argument processing:

  Real CrustNu  = 0.8;          // Top Layer Defaults
  Real CrustEps = 0.05;
  Real CrustA   = 0.50;
  Real CrustK   = 0.5;
  Real CrustQ   = 1.0/0.0;

  Real CrustVp  = 6.40;
  Real CrustVs  = 3.63;
  Real CrustRho = 2.83;
  Real CrustZb  = -21.5;

  Real MantNu  = CrustNu;       // Mantle region defaults
  Real MantEps = 0.04;
  Real MantA   = 1.0;
  Real MantK   = CrustK;
  Real MantQ   = CrustQ;

  Real MantVp  = 6.40;
  Real MantVs  = 3.63;
  Real MantRho = 2.83;
  Real MantZb  = -100.0;

  // Now check args and assign whatever regional params were
  // specified, depending on how many args passed:

  switch(args.size()) {
  case 0:   // No args given, just use defaults computed above
            //
    break;

  case 18:  // User specified nu,eps,a,k,Q AND elastic velocities
            //
    CrustVp  = args.at(10);
    CrustVs  = args.at(11);
    CrustRho = args.at(12);
    CrustZb  = args.at(13);
    MantVp  = args.at(14);
    MantVs  = args.at(15);
    MantRho = args.at(16);
    MantZb  = args.at(17);
    // fallthrough:

  case 10:  // Specified nu,eps,a,k,Q for three regions
            //
    CrustNu  = args.at(0);
    CrustEps = args.at(1);
    CrustA   = args.at(2);
    CrustK   = args.at(3);
    CrustQ   = args.at(4);
    MantNu  = args.at(5);
    MantEps = args.at(6);
    MantA   = args.at(7);
    MantK   = args.at(8);
    MantQ   = args.at(9);
    break;

  default:  // Unrecognized pattern of values
            //
    std::cerr << "Error: wrong number of model args passed "
              << "to compiled-in grid-building function.\n";
    exit(1);  // TODO: Raise a meaningful exception instead
    break;

  }

  Real EpiZ[] = { 0.0, CrustZb, MantZb };
  Real EASZ[] = { 0.0, CrustZb, MantZb };
  Real NORZ[] = { 0.0, CrustZb, MantZb };
  bool RT = !((CrustVp==MantVp)&&(CrustVs==MantVs)&&(CrustRho==MantRho));

  HetSpec HS;
  Q Q;

  // :::::::::::::::::::::::::
  // ::: *** TOP LAYER *** :::    (Very top layer of Earth model)
  // :::::::::::::::::::::::::

  //     Grid is specified as a series of "sheets", defining the
  //     interface between two volumetric layers.  If only one
  //     attribute set is specified, it applies to layers above AND
  //     below.  If TWO attribute sets are specified, the first
  //     applies to the layer above, and the second to layer
  //     below. (This defines a discontinuous "jump" in properties, as
  //     opposed to a waypoint in a piecewise continuous function of
  //     properties.  Such a discontinuity triggers R/T handling at
  //     the interface.)
  //

  // Sheet 0: between: Surface/Air       <-- Name of "layer" above
  //              and: Top-Layer         <-- Name of "layer" below

  HS = HSneak(CrustNu,CrustEps,CrustA,CrustK);
  Q  = QmQk(CrustQ);    // Parameter values in the top layer

  gr.WNode(0,0,0).SetLocation(EpiX,EpiY,EpiZ[0]);
  gr.WNode(1,0,0).SetLocation(EASX,EASY,EASZ[0]);
  gr.WNode(2,0,0).SetLocation(NORX,NORY,NORZ[0]);

  gr.WNode(0,0,0).SetAttributes(VpVs(CrustVp, CrustVs), CrustRho, Q, HS);

  // Sheet 1: Between: Top-Layer
  //              and: Bottom-Layer
  //

  gr.WNode(0,0,1).SetLocation(EpiX,EpiY,EpiZ[1]);
  gr.WNode(1,0,1).SetLocation(EASX,EASY,EASZ[1]);
  gr.WNode(2,0,1).SetLocation(NORX,NORY,NORZ[1]);

  // ::::::::::::::::::::::::::::
  // ::: *** BOTTOM LAYER *** :::
  // ::::::::::::::::::::::::::::

  //  |  Specify TWO attribute sets for this sheet.  First set defines
  // \|/ bottom of layer above; second set defines top of layer below.
  //  V  Double-valued attrs mark a discontinuity and make an R/T interface.
  if (RT) gr.WNode(0,0,1).SetAttributes(VpVs(CrustVp,CrustVs),CrustRho, Q, HS);
  HS = HSneak(MantNu, MantEps, MantA, MantK);
  Q  = QmQk(MantQ);
  gr.WNode(0,0,1).SetAttributes(VpVs(MantVp, MantVs), MantRho, Q, HS);

  // Sheet 2: Bottom of Model
  //

  gr.WNode(0,0,2).SetLocation(EpiX,EpiY,EpiZ[2]);
  gr.WNode(1,0,2).SetLocation(EASX,EASY,EASZ[2]);
  gr.WNode(2,0,2).SetLocation(NORX,NORY,NORZ[2]);

  gr.WNode(0,0,2).SetAttributes(VpVs(MantVp, MantVs), MantRho, Q, HS);

}
