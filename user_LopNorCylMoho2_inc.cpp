// user_LopNorCylMoho_inc.cpp
//
//   Defines a specific model.
//   #include this file in user.cpp
//

//////
// LopNorCylinderMoho2()
//
//   Adds sub-moho structure, alternate velocity profile, a bit more gradual
//
//
void LopNorCylinderMoho2(Grid & gr, const std::vector<Real> & args) {

  using Elastic::Velocity;
  using Elastic::VpVs;
  using Elastic::Q;
  using Elastic::QmQk;
  using Elastic::HetSpec;
  using Elastic::HSneak;

  // Grid Dimension:

  gr.SetSize(3,1,29);   // 29 Grid Sheets -  First 12 define laterally-varying
  gr.SetIndexBase(0);   //                   topography, rest define mantle 
                        //                   from AK135

  // Site locations:

  Real LopX =  492.31;
  Real LopY = -263.65;
  Real LopZ[] = { 1.050, /* . . . . . . . . . . . . */  // Sediments
                  0.563, -18.901,  /* . . . . . . . */  // Crust Proper
                  -38.365, -40.365, -42.365, -44.365,   // Crust Moho-Complex
                  -46.365, -48.365, -50.365, -52.365,   // Mantle Ramp
                  -54.365 }; /* . . . . . . . . . . */  // Mantle Proper
  Real MAKX = -102.27;
  Real MAKY =  430.84;
  Real MAKZ[] = { 0.600,
                  0.118, -16.743,
                  -33.122, -35.122, -37.122, -39.122,
                  -41.122, -43.122, -45.122, -47.122,
                  -49.122 };
  Real WUSX = -390.04;
  Real WUSY = -167.18;
  Real WUSZ[] = { 1.457, 
                  0.963, -18.812,
                  -38.587, -40.587, -42.587, -44.587,
                  -46.587, -48.587, -50.587, -52.587,
                  -54.587 };

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

  Real MohoNu  = CrustNu;       // Moho/transition region defaults
  Real MohoEps = CrustEps;
  Real MohoA   = CrustA;
  Real MohoK   = CrustK;
  Real MohoQ   = CrustQ;

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

  case 20:  // Specified nu,eps,a,k,Q for four regions
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
    MohoNu  = args.at(15);
    MohoEps = args.at(16);
    MohoA   = args.at(17);
    MohoK   = args.at(18);
    MohoQ   = args.at(19);
    break;

  default:  // Unrecognized pattern of values
            //
    std::cerr << "Error: wrong number of model args passed "
              << "to compiled-in grid-building function.\n";
    exit(1);  // TODO: Raise a meaningful exception instead
    break;

  }

  HetSpec HS;
  HetSpec HS_abv;   // Used at discontinuities, refers to volume above
  Elastic::Q Q;
  Elastic::Q Q_abv;

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

  HS_abv = HS; Q_abv = Q;
  HS = HSneak(CrustNu, CrustEps, CrustA, CrustK);
  Q  = QmQk(CrustQ);    // Parameter values in the "crust" layers

  gr.WNode(0,0,1).SetLocation(LopX,LopY,LopZ[1]);
  gr.WNode(1,0,1).SetLocation(MAKX,MAKY,MAKZ[1]);
  gr.WNode(2,0,1).SetLocation(WUSX,WUSY,WUSZ[1]);

  //  |  Specify TWO attribute sets for this sheet.  First set defines
  // \|/ bottom of layer above (soft-sediments); second set defines top of
  //  V  layer below (upper-crust).  Double-valued attrs make R/T iface.
  gr.WNode(0,0,1).SetAttributes(VpVs(2.50, 1.20), 2.10, Q_abv, HS_abv);
  gr.WNode(0,0,1).SetAttributes(VpVs(6.13, 3.53), 2.75, Q, HS);

  // Sheet 2: Between Upper Crust, and
  //                  Lower Crust          <-- Attributes apply here

  gr.WNode(0,0,2).SetLocation(LopX,LopY,LopZ[2]);
  gr.WNode(1,0,2).SetLocation(MAKX,MAKY,MAKZ[2]);
  gr.WNode(2,0,2).SetLocation(WUSX,WUSY,WUSZ[2]);

  gr.WNode(0,0,2).SetAttributes(VpVs(6.40, 3.63), 2.83, Q, HS);

  // ::::::::::::::::::::::::::
  // ::: *** TRANSITION *** :::
  // ::::::::::::::::::::::::::

  // Sheet 3: Between Lower Crust, and
  //                  Moho-complexity 1       (Double-valued for R/T)

  HS_abv = HS; Q_abv = Q; // (Remembered for region-above attributes)
  HS = HSneak(MohoNu, MohoEps, MohoA, MohoK); // Properties apply to
  Q  = QmQk(MohoQ);                           // Moho-complexity region

  gr.WNode(0,0,3).SetLocation(LopX,LopY,LopZ[3]);
  gr.WNode(1,0,3).SetLocation(MAKX,MAKY,MAKZ[3]);
  gr.WNode(2,0,3).SetLocation(WUSX,WUSY,WUSZ[3]);

  gr.WNode(0,0,3).SetAttributes(VpVs(6.40, 3.63), 3.10, Q_abv, HS_abv);
  gr.WNode(0,0,3).SetAttributes(VpVs(7.081, 3.885), 3.15, Q, HS);

  // Sheet 4: Between Moho-plex 1, and
  //                  Moho-plex 2              (Double Valued for R/T)

  gr.WNode(0,0,4).SetLocation(LopX,LopY,LopZ[4]);
  gr.WNode(1,0,4).SetLocation(MAKX,MAKY,MAKZ[4]);
  gr.WNode(2,0,4).SetLocation(WUSX,WUSY,WUSZ[4]);

  gr.WNode(0,0,4).SetAttributes(VpVs(7.081, 3.885), 3.15, Q, HS);
  gr.WNode(0,0,4).SetAttributes(VpVs(6.991, 3.835), 3.13, Q, HS);

  // Sheet 5: Between Moho-plex 2, and
  //                  Moho-plex 3              (Double Valued for R/T)

  gr.WNode(0,0,5).SetLocation(LopX,LopY,LopZ[5]);
  gr.WNode(1,0,5).SetLocation(MAKX,MAKY,MAKZ[5]);
  gr.WNode(2,0,5).SetLocation(WUSX,WUSY,WUSZ[5]);

  gr.WNode(0,0,5).SetAttributes(VpVs(6.991, 3.835), 3.13, Q, HS);
  gr.WNode(0,0,5).SetAttributes(VpVs(7.291, 3.985), 3.22, Q, HS);

  // Sheet 6: Between Moho-plex 3, and
  //                  Moho-plex 4              (Double Valued for R/T)

  gr.WNode(0,0,6).SetLocation(LopX,LopY,LopZ[6]);
  gr.WNode(1,0,6).SetLocation(MAKX,MAKY,MAKZ[6]);
  gr.WNode(2,0,6).SetLocation(WUSX,WUSY,WUSZ[6]);

  gr.WNode(0,0,6).SetAttributes(VpVs(7.291, 3.985), 3.22, Q, HS);
  gr.WNode(0,0,6).SetAttributes(VpVs(7.141, 3.875), 3.15, Q, HS);


  // :::::::::::::::::::::::::::
  // ::: *** MANTLE RAMP *** :::
  // :::::::::::::::::::::::::::

  // Sheet 7: Between Moho-plex 4, and
  //                  Mantle Ramp 1            (Double Valued for R/T)

  HS_abv = HS; Q_abv = Q; // (Remembered for region-above attributes)
  HS = HSneak(MohoNu, MohoEps, MohoA, MohoK); // Properties apply to
  Q  = QmQk(MohoQ);                           // Mantle-ramp region

  gr.WNode(0,0,7).SetLocation(LopX,LopY,LopZ[7]);
  gr.WNode(1,0,7).SetLocation(MAKX,MAKY,MAKZ[7]);
  gr.WNode(2,0,7).SetLocation(WUSX,WUSY,WUSZ[7]);

  gr.WNode(0,0,7).SetAttributes(VpVs(7.141, 3.875), 3.22, Q_abv, HS_abv);
  gr.WNode(0,0,7).SetAttributes(VpVs(7.624, 4.183), 3.502, Q, HS);

  // Sheet 8: Between Mantle Ramp 1, and
  //                  Mantle Ramp 2            (Double Valued for R/T)

  gr.WNode(0,0,8).SetLocation(LopX,LopY,LopZ[8]);
  gr.WNode(1,0,8).SetLocation(MAKX,MAKY,MAKZ[8]);
  gr.WNode(2,0,8).SetLocation(WUSX,WUSY,WUSZ[8]);

  gr.WNode(0,0,8).SetAttributes(VpVs(7.624, 4.183), 3.502, Q, HS);
  gr.WNode(0,0,8).SetAttributes(VpVs(7.728, 4.257), 3.502, Q, HS);

  // Sheet 9: Between Mantle Ramp 2, and
  //                  Mantle Ramp 3            (Double Valued for R/T)

  gr.WNode(0,0,9).SetLocation(LopX,LopY,LopZ[9]);
  gr.WNode(1,0,9).SetLocation(MAKX,MAKY,MAKZ[9]);
  gr.WNode(2,0,9).SetLocation(WUSX,WUSY,WUSZ[9]);

  gr.WNode(0,0,9).SetAttributes(VpVs(7.728, 4.257), 3.502, Q, HS);
  gr.WNode(0,0,9).SetAttributes(VpVs(7.832, 4.331), 3.502, Q, HS);

  // Sheet 10: Between Mantle Ramp 3, and
  //                   Mantle Ramp 4           (Double Valued for R/T)

  gr.WNode(0,0,10).SetLocation(LopX,LopY,LopZ[10]);
  gr.WNode(1,0,10).SetLocation(MAKX,MAKY,MAKZ[10]);
  gr.WNode(2,0,10).SetLocation(WUSX,WUSY,WUSZ[10]);

  gr.WNode(0,0,10).SetAttributes(VpVs(7.832, 4.331), 3.502, Q, HS);
  gr.WNode(0,0,10).SetAttributes(VpVs(7.936, 4.406), 3.502, Q, HS);


  // ::::::::::::::::::::::
  // ::: *** MANTLE *** :::
  // ::::::::::::::::::::::

  // Sheet 11: Between Mantle-Ramp 4, and
  //                   Mantle Proper           (Double Valued for R/T)

  HS_abv = HS; Q_abv = Q; // (Remembered for region-above attributes)
  HS = HSneak(MantNu, MantEps, MantA, MantK); // Properties apply to
  Q  = QmQk(MantQ);                           // Mantle region

  gr.WNode(0,0,11).SetLocation(LopX,LopY,LopZ[11]);
  gr.WNode(1,0,11).SetLocation(MAKX,MAKY,MAKZ[11]);
  gr.WNode(2,0,11).SetLocation(WUSX,WUSY,WUSZ[11]);

  gr.WNode(0,0,11).SetAttributes(VpVs(7.936, 4.406), 3.502, Q_abv, HS_abv);
  gr.WNode(0,0,11).SetAttributes(VpVs(8.040, 4.480), 3.502, Q, HS);

  // Sheet 12: Between First Mantle Layer, and
  //                   Second Mantle Layer     (Double Valued for R/T)

  MantleZ = -80.0;
  gr.WNode(0,0,12).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,12).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,12).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,12).SetAttributes(VpVs(8.040, 4.480), 3.502, Q, HS);
  gr.WNode(0,0,12).SetAttributes(VpVs(8.045, 4.490), 3.502, Q, HS);

  // Sheet 13:

  MantleZ = -120.0;
  gr.WNode(0,0,13).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,13).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,13).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,13).SetAttributes(VpVs(8.0505, 4.5000), 3.4268, Q, HS);

  // Sheet 14:

  MantleZ = -165.0;
  gr.WNode(0,0,14).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,14).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,14).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,14).SetAttributes(VpVs(8.1750, 4.5090), 3.3711, Q, HS);

  // Sheet 15:

  MantleZ = -210.0;
  gr.WNode(0,0,15).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,15).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,15).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,15).SetAttributes(VpVs(8.3007, 4.5184), 3.3243, Q, HS);

  // Sheet 16:

  MantleZ = -260.0;
  gr.WNode(0,0,16).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,16).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,16).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,16).SetAttributes(VpVs(8.4822, 4.6094), 3.3663, Q, HS);

  // Sheet 17:

  MantleZ = -310.0;
  gr.WNode(0,0,17).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,17).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,17).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,17).SetAttributes(VpVs(8.6650, 4.6964), 3.4110, Q, HS);

  // Sheet 18:

  MantleZ = -360.0;
  gr.WNode(0,0,18).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,18).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,18).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,18).SetAttributes(VpVs(8.8476, 4.7832), 3.4577, Q, HS);

  // Sheet 19:    -- An R/T (double-valued) Layer --

  MantleZ = -410.0;
  gr.WNode(0,0,19).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,19).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,19).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,19).SetAttributes(VpVs(9.0302, 4.8702), 3.5068, Q, HS);
  gr.WNode(0,0,19).SetAttributes(VpVs(9.3601, 5.0806), 3.9317, Q, HS);

  // Sheet 20:

  MantleZ = -460.0;
  gr.WNode(0,0,20).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,20).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,20).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,20).SetAttributes(VpVs(9.5280, 5.1864), 3.9273, Q, HS);

  // Sheet 21:

  MantleZ = -510.0;
  gr.WNode(0,0,21).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,21).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,21).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,21).SetAttributes(VpVs(9.6962, 5.2922), 3.9233, Q, HS);

  // Sheet 22:

  MantleZ = -560.0;
  gr.WNode(0,0,22).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,22).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,22).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,22).SetAttributes(VpVs(9.8640, 5.3989), 3.9218, Q, HS);

  // Sheet 23:

  MantleZ = -610.0;
  gr.WNode(0,0,23).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,23).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,23).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,23).SetAttributes(VpVs(10.0320, 5.5047), 3.9206, Q, HS);

  // Sheet 24:    -- An R/T (double-valued) Layer --

  MantleZ = -660.0;
  gr.WNode(0,0,24).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,24).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,24).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,24).SetAttributes(VpVs(10.2000, 5.6104), 3.9201, Q, HS);
  gr.WNode(0,0,24).SetAttributes(VpVs(10.7909, 5.9607), 4.2387, Q, HS);

  // Sheet 25:

  MantleZ = -710.0;
  gr.WNode(0,0,25).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,25).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,25).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,25).SetAttributes(VpVs(10.9222, 6.0898), 4.2986, Q, HS);

 
  // Sheet 26:

  MantleZ = -760.0;
  gr.WNode(0,0,26).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,26).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,26).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,26).SetAttributes(VpVs(11.0553, 6.2100), 4.3565, Q, HS);

  // Sheet 27:

  MantleZ = -809.5;
  gr.WNode(0,0,27).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,27).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,27).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,27).SetAttributes(VpVs(11.1355, 6.2424), 4.4118, Q, HS);

  // Sheet 28:

  MantleZ = -859.0;
  gr.WNode(0,0,28).SetLocation(LopX,LopY,MantleZ);
  gr.WNode(1,0,28).SetLocation(MAKX,MAKY,MantleZ);
  gr.WNode(2,0,28).SetLocation(WUSX,WUSY,MantleZ);

  gr.WNode(0,0,28).SetAttributes(VpVs(11.2228, 6.2799), 4.4650, Q, HS);


}
