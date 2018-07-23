// user_LopNorCyl_inc.cpp
//
//   Defines a specific model.
//   #include this file in user.cpp
//

//////
// LopNorCylinder()
//
//   Models the region surrounding Lop Nor test site and stations MAK and
//   WUS as a cylinder grid/model.
//
//   This grid is for generating figures/data for ABQ2014 (June)
//
//   Scattering parameters and potentially Q values are specified in
//   the 'args' argument and are provided by the user by supplying the
//   --model-compiled-args switch.  The switch takes a list of Real
//   values whichdefine the parameters in three regions: the sediments
//   region, the crust region, and the mantle region.  The switch
//   takes 9, 12, or fifteen values. If 9, then they specify nu,eps,a
//   for sediments, crust, mantle, respectively.  If 12, they specify
//   nu,eps,a,kappa, and if 15 they specify nu,eps,a,kappa,Q.  E.g.:
//
//   --model-compiled-args=nu1,eps1,a1,k1,nu2,eps2,a2,ka,nu3,eps3,a3,k3
//
//   which would give the nu, eps, a, and kappa values, but would omit
//   the Q values.  Q is taken to be infinite by default, unless
//   specified otherwise by the user.
//
void LopNorCylinder(Grid & gr, const std::vector<Real> & args) {

  using Elastic::Velocity;
  using Elastic::VpVs;
  using Elastic::Q;
  using Elastic::QmQk;
  using Elastic::HetSpec;
  using Elastic::HSneak;

  // Grid Dimension:

  gr.SetSize(3,1,22);   // 22 Grid Layers -  5 define the crust, the
  gr.SetIndexBase(0);   //                   rest define the upper mantle


  // Site locations:

  Real LopX =  492.31;
  Real LopY = -263.65;
  Real LopZ[] = { 1.050, 0.563, -18.901, -38.365, -47.610 };
  Real MAKX = -102.27;
  Real MAKY =  430.84;
  Real MAKZ[] = { 0.600, 0.118, -16.743, -33.122, -43.720 };
  Real WUSX = -390.04;
  Real WUSY = -167.18;
  Real WUSZ[] = { 1.457, 0.963, -18.812, -38.587, -47.980 };

  Real MantleZ;         // Used below; zero slope on Mantle layers

  // Scattering Params, defaults and argument processing:

  Real CrustNu  = 0.8;          // Crustal Region Defaults
  Real CrustEps = 0.05;         // 
  Real CrustA   = 0.50;         // 
  Real CrustK   = 0.5;          // 
  Real CrustQ   = 1.0/0.0;      // 

  Real SediNu  = CrustNu;       // Sediments layer defaults
  Real SediEps = 0.06;
  Real SediA   = 0.25;
  Real SediK   = CrustK;
  Real SediQ   = CrustQ;

  Real MantNu  = CrustNu;       // Mantle region defaults
  Real MantEps = 0.04;
  Real MantA   = 1.0;
  Real MantK   = CrustK;
  Real MantQ   = CrustQ;

  // Now check args and assign whatever regional params were
  // specified, depending on how many args passed:

  switch(args.size()) {
  case 0:   // No args given, just use defaults computed above
            //
    break;

  case 9:   // User specified nu,eps,a for the three regions.
            // As for k and Q, use defaults computed above.
    SediNu  = args.at(0);
    SediEps = args.at(1);
    SediA   = args.at(2);
    CrustNu  = args.at(3);
    CrustEps = args.at(4);
    CrustA   = args.at(5);
    MantNu  = args.at(6);
    MantEps = args.at(7);
    MantA   = args.at(8);
    break;

  case 12:  // User specified nu,eps,a,k for the three regions.
            // For Q, use defaults computed above.
    SediNu  = args.at(0);
    SediEps = args.at(1);
    SediA   = args.at(2);
    SediK   = args.at(3);
    CrustNu  = args.at(4);
    CrustEps = args.at(5);
    CrustA   = args.at(6);
    CrustK   = args.at(7);
    MantNu  = args.at(8);
    MantEps = args.at(9);
    MantA   = args.at(10);
    MantK   = args.at(11);
    break;

  case 15:  // Specified nu,eps,a,k,Q for three regions
            //
    SediNu  = args.at(0);
    SediEps = args.at(1);
    SediA   = args.at(2);
    SediK   = args.at(3);
    SediQ   = args.at(4);
    CrustNu  = args.at(5);
    CrustEps = args.at(6);
    CrustA   = args.at(7);
    CrustK   = args.at(8);
    CrustQ   = args.at(9);
    MantNu  = args.at(10);
    MantEps = args.at(11);
    MantA   = args.at(12);
    MantK   = args.at(13);
    MantQ   = args.at(14);
    break;

  default:  // Unrecognized pattern of values
            //
    std::cerr << "Error: wrong number of model args passed "
              << "to compiled-in grid-building function.\n";
    exit(1);  // TODO: Raise a meaningful exception instead
    break;

  }

  HetSpec HS;
  Q Q;
   
  // :::::::::::::::::::::::::
  // ::: *** SEDIMENTS *** :::    (Very top layer of Earth model)
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
  //              and: Soft-Sediments    <-- Name of "layer" below

  HS = HSneak(SediNu,SediEps,SediA,SediK);
  Q  = QmQk(SediQ);     // Parameter values in the "sediments" layer

  gr.WNode(0,0,0).SetLocation(LopX,LopY,LopZ[0]);
  gr.WNode(1,0,0).SetLocation(MAKX,MAKY,MAKZ[0]);
  gr.WNode(2,0,0).SetLocation(WUSX,WUSY,WUSZ[0]);

  //                                 v_p   v_s    rho
  gr.WNode(0,0,0).SetAttributes(VpVs(2.50, 1.20), 2.10, Q, HS);


  // :::::::::::::::::::::
  // ::: *** CRUST *** :::
  // :::::::::::::::::::::

  // Sheet 1: Between: Soft-Sediments
  //              and: Upper-Crust
  //                                        

  HS = HSneak(CrustNu, CrustEps, CrustA, CrustK);
  Q  = QmQk(CrustQ);    // Parameter values in the "crust" layers

  gr.WNode(0,0,1).SetLocation(LopX,LopY,LopZ[1]);
  gr.WNode(1,0,1).SetLocation(MAKX,MAKY,MAKZ[1]);
  gr.WNode(2,0,1).SetLocation(WUSX,WUSY,WUSZ[1]);

  //  |  Specify TWO attribute sets for this sheet.  First set defines
  // \|/ bottom of layer above (soft-sediments); second set defines top of
  //  V  layer below (upper-crust).  Double-valued attrs make R/T iface.
  gr.WNode(0,0,1).SetAttributes(VpVs(2.50, 1.20), 2.10, Q, HS);
  gr.WNode(0,0,1).SetAttributes(VpVs(6.13, 3.53), 2.75, Q, HS);

  // Sheet 2: Between Upper Crust, and
  //                  Middle Crust          <-- Attributes apply here

  gr.WNode(0,0,2).SetLocation(LopX,LopY,LopZ[2]);
  gr.WNode(1,0,2).SetLocation(MAKX,MAKY,MAKZ[2]);
  gr.WNode(2,0,2).SetLocation(WUSX,WUSY,WUSZ[2]);

  gr.WNode(0,0,2).SetAttributes(VpVs(6.40, 3.63), 2.83, Q, HS);

  // Sheet 3: Between Middle Crust, and
  //                  Lower Crust           <-- Attributes apply here

  gr.WNode(0,0,3).SetLocation(LopX,LopY,LopZ[3]);
  gr.WNode(1,0,3).SetLocation(MAKX,MAKY,MAKZ[3]);
  gr.WNode(2,0,3).SetLocation(WUSX,WUSY,WUSZ[3]);

  gr.WNode(0,0,3).SetAttributes(VpVs(7.23, 4.00), 3.10, Q, HS);


  // ::::::::::::::::::::::
  // ::: *** MANTLE *** :::
  // ::::::::::::::::::::::

  // Sheet 4: Between Lower Crust, and      <-- First Attribute Set (for R/T)
  //                  Mantle                <-- Second Attribute Set
  //                                        Double-Valued Attrs make R/T iface
  //   Transition from Crust to Mantle
  //

  gr.WNode(0,0,4).SetLocation(LopX,LopY,LopZ[4]);
  gr.WNode(1,0,4).SetLocation(MAKX,MAKY,MAKZ[4]);
  gr.WNode(2,0,4).SetLocation(WUSX,WUSY,WUSZ[4]);

  gr.WNode(0,0,4).SetAttributes(VpVs(7.23, 4.00), 3.10, Q, HS);

  HS = HSneak(MantNu, MantEps, MantA, MantK);
  Q  = QmQk(MantQ);     // Parameter values in the "mantle" layers
                        //
                        //    |
                        //    |
                        //   \|/  values used for UNDERSIDE of 
                        //    V   gridnode sheet (top of first mantle layer) 

  gr.WNode(0,0,4).SetAttributes(VpVs(8.07, 4.63), 3.35, Q, HS);

  // Sheet 5: Between First Mantle Layer, and
  //                  Second Mantle Layer         (Double Valued for R/T)

  MantleZ = -80.0;
  gr.WNode(0,0,5).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,5).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,5).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,5).SetAttributes(VpVs(8.040, 4.480), 3.502, Q, HS);
  gr.WNode(0,0,5).SetAttributes(VpVs(8.045, 4.490), 3.502, Q, HS);

  // Sheet 6:

  MantleZ = -120.0;
  gr.WNode(0,0,6).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,6).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,6).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,6).SetAttributes(VpVs(8.0505, 4.5000), 3.4268, Q, HS);

  // Sheet 7:

  MantleZ = -165.0;
  gr.WNode(0,0,7).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,7).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,7).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,7).SetAttributes(VpVs(8.1750, 4.5090), 3.3711, Q, HS);

  // Sheet 8:

  MantleZ = -210.0;
  gr.WNode(0,0,8).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,8).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,8).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,8).SetAttributes(VpVs(8.3007, 4.5184), 3.3243, Q, HS);

  // Sheet 9:

  MantleZ = -260.0;
  gr.WNode(0,0,9).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,9).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,9).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,9).SetAttributes(VpVs(8.4822, 4.6094), 3.3663, Q, HS);

  // Sheet 10:

  MantleZ = -310.0;
  gr.WNode(0,0,10).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,10).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,10).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,10).SetAttributes(VpVs(8.6650, 4.6964), 3.4110, Q, HS);

  // Sheet 11:

  MantleZ = -360.0;
  gr.WNode(0,0,11).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,11).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,11).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,11).SetAttributes(VpVs(8.8476, 4.7832), 3.4577, Q, HS);

  // Sheet 12:    -- An R/T (double-valued) Layer --

  MantleZ = -410.0;
  gr.WNode(0,0,12).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,12).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,12).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,12).SetAttributes(VpVs(9.0302, 4.8702), 3.5068, Q, HS);
  gr.WNode(0,0,12).SetAttributes(VpVs(9.3601, 5.0806), 3.9317, Q, HS);

  // Sheet 13:

  MantleZ = -460.0;
  gr.WNode(0,0,13).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,13).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,13).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,13).SetAttributes(VpVs(9.5280, 5.1864), 3.9273, Q, HS);

  // Sheet 14:

  MantleZ = -510.0;
  gr.WNode(0,0,14).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,14).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,14).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,14).SetAttributes(VpVs(9.6962, 5.2922), 3.9233, Q, HS);

  // Sheet 15:

  MantleZ = -560.0;
  gr.WNode(0,0,15).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,15).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,15).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,15).SetAttributes(VpVs(9.8640, 5.3989), 3.9218, Q, HS);

  // Sheet 16:

  MantleZ = -610.0;
  gr.WNode(0,0,16).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,16).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,16).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,16).SetAttributes(VpVs(10.0320, 5.5047), 3.9206, Q, HS);

  // Sheet 17:    -- An R/T (double-valued) Layer --

  MantleZ = -660.0;
  gr.WNode(0,0,17).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,17).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,17).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,17).SetAttributes(VpVs(10.2000, 5.6104), 3.9201, Q, HS);
  gr.WNode(0,0,17).SetAttributes(VpVs(10.7909, 5.9607), 4.2387, Q, HS);

  // Sheet 18:

  MantleZ = -710.0;
  gr.WNode(0,0,18).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,18).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,18).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,18).SetAttributes(VpVs(10.9222, 6.0898), 4.2986, Q, HS);

  // Sheet 19:

  MantleZ = -760.0;
  gr.WNode(0,0,19).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,19).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,19).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,19).SetAttributes(VpVs(11.0553, 6.2100), 4.3565, Q, HS);

  // Sheet 20:

  MantleZ = -809.5;
  gr.WNode(0,0,20).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,20).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,20).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,20).SetAttributes(VpVs(11.1355, 6.2424), 4.4118, Q, HS);

  // Sheet 21:

  MantleZ = -859.0;
  gr.WNode(0,0,21).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,21).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,21).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,21).SetAttributes(VpVs(11.2228, 6.2799), 4.4650, Q, HS);


}
