// user_NSCP_inc.cpp
//
//   Defines a specific model.
//   #include this file in user.cpp
//

//////
// Crustal Pinch Model
//
// Uses RAE (Range, Azimuth, Elevation) coordinate scheme. For version
// built with XYZ coordinates, revert to revision 912.
//
void CrustPinchWCG(Grid & gr, const std::vector<Real> & args) {

  using Elastic::Velocity;
  using Elastic::VpVs;
  using Elastic::Q;
  using Elastic::QmQk;
  using Elastic::HetSpec;
  using Elastic::HSneak;

  Real sc_nu  = 0.8;    // Default scat params -
  Real sc_eps = 0.01;   // we'll scan args for user values
  Real sc_a   = 4.00;   //
  Real sc_k = 0.8;      //
  Real q = 200;         //

  HetSpec HSSe = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Sedi
  HetSpec HSCr = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Crust
  HetSpec HSPi = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Crust Pinched
  HetSpec HSMo = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Moho
  HetSpec HSMa = HSneak(sc_nu, sc_eps, sc_a, sc_k); // Mantle

  Q QSe = QmQk(q);          //
  Q QCr = QmQk(q);          //
  Q QPi = QmQk(q);          //
  Q QMo = QmQk(q);          //
  Q QMa = QmQk(q);          //

  Real SediThick  = 2.0;    // Thickness of Sediments layer
  Real CrustThick = 30.0;   // Thickness of (unpinched) Crust layer
  Real MohoThick  = 10.0;   // Thickness of Moho layer
                            //
  Real PinchFrac  = 0.70;   // Thickness of "pinched" crust as
                            // fraction of un-pinched thickness.
  Real DepthFrac  = 0.50;   // Depth of pinched slab.
                            //   0.0: Aligned top
                            //   0.5: Common midline
                            //   1.0: Aligned bottom
  Real SediFrac   = 1.0;    // Controls sediments thickness in pinched
                            // region
  Real MohoFrac   = 1.0;    // Controls moho thickness in pinched
                            // region.

  // Now check args array:

  switch(args.size()) {
  case 0:   // No args given, just use defaults already constructed
            //
    break;

  case 32:  // 4 x neakq + layer thicks + pinch fractions
            //
    PinchFrac = args.at(28);
    DepthFrac = args.at(29);
    SediFrac  = args.at(30);
    MohoFrac  = args.at(31);
    /* FALL THROUGH  - no break; */

  case 28:  // 4 x neakq + layer thicknesses
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

  Count nZ = 8;     // How many z-depths
  Count nR = 14;    // How many range points
  Count nAzis = 6;  // How many azimuths

  Real Azis[] =   {45.0, 56.25, 78.75, 101.25, 123.75, 135.0};// Azis in degrees
  Real AzisCR[] = { 5.0, 39.00, 73.00, 107.00, 141.00, 175.0};// Close-range spread

  Real ZBase[]  = {   0.0,              // Un-pinched Z-depths
                    -SediThick,
                    -(SediThick+CrustThick), 
                    -(SediThick+CrustThick+MohoThick),
                      -80,    // Very subtle discontinuity here.
                      -120,
                      -210,
                      -360
  };


  Real PinchThick = PinchFrac*CrustThick;     // Compute pinched z-depths
  Real PinchTop = - SediThick                 // ******
                  - ((CrustThick-PinchThick)  // ***
                     * DepthFrac);
  Real SediFill = std::max(SediFrac*(-SediThick - PinchTop),Real(0));
  Real SediPinchTop = (PinchTop + (SediThick+SediFill));
  Real MohoPinchTop = (PinchTop - PinchThick);
  Real MantlePinchTop = (MohoPinchTop - MohoFrac*MohoThick);

  Real ZPinch[] = { SediPinchTop,       // Pinched Z-depths
                    PinchTop, 
                    MohoPinchTop,
                    MantlePinchTop,
                    ZBase[4],
                    ZBase[5],
                    ZBase[6],
                    ZBase[7]
 };


  gr.SetSize(nR,nAzis,nZ);      // Sets index bounds
  gr.SetIndexBase(0);           // When addressing nodes, use base 0
  gr.SetMapping(Grid::GC_RAE, Grid::GC_CURVED);

  for (Index iaz = 0; iaz < nAzis; iaz++) {

    for (Index iz = 0; iz < nZ; iz++) {

      Real azi  = Azis[iaz];              // Azimuth value
      Real azin = Azis[nAzis-iaz-1];      // picked in reverse order
      Real azicr  = AzisCR[iaz];          // Azi's for close-range
      Real azicrn = AzisCR[nAzis-iaz-1];  // Close range reverse order
      //Real aziir1 = (azicr+azi)/2;        // Intermediate range
      //Real aziir2 = (azicr+2*azi)/3;      //  ''

      gr.WNode( 0,iaz,iz).SetLocation ( -120,  azin,  ZBase[iz]  );
      gr.WNode( 1,iaz,iz).SetLocation (  -60, azicrn, ZBase[iz]  );
      gr.WNode( 2,iaz,iz).SetLocation (   60,  azicr, ZBase[iz]  );
      gr.WNode( 3,iaz,iz).SetLocation (  120,   azi,  ZBase[iz]  );
      gr.WNode( 4,iaz,iz).SetLocation (  220,   azi,  ZBase[iz]  );
      gr.WNode( 5,iaz,iz).SetLocation (  310,   azi,  ZBase[iz]  );

      gr.WNode( 6,iaz,iz).SetLocation (  370,   azi,  ZPinch[iz]  );
      gr.WNode( 7,iaz,iz).SetLocation (  420,   azi,  ZPinch[iz]  );
      gr.WNode( 8,iaz,iz).SetLocation (  470,   azi,  ZPinch[iz]  );

      gr.WNode( 9,iaz,iz).SetLocation (  530,   azi,  ZBase[iz]  );
      gr.WNode(10,iaz,iz).SetLocation (  650,   azi,  ZBase[iz]  );
      gr.WNode(11,iaz,iz).SetLocation (  770,   azi,  ZBase[iz]  );
      gr.WNode(12,iaz,iz).SetLocation (  890,   azi,  ZBase[iz]  );
      gr.WNode(13,iaz,iz).SetLocation ( 1020,   azi,  ZBase[iz]  );

    }

    for (Index ir = 0; ir < nR; ir++) { // Attributes vary only in depth index

      Q QCrtmp; HetSpec HSCrtmp;
      if ((ir>=6) && (ir<8)) {  // If defining a "pinched" block:
        QCrtmp = QPi;           // Then use alternate HS, Q for crust layer.
        HSCrtmp = HSPi;         //
      } else {
        QCrtmp = QCr;           // Otherwise regular HS and Q.
        HSCrtmp = HSCr;         //
      }

      gr.WNode(ir,iaz,0) .SetAttributes( VpVs( 4.50, 2.60), 2.20, QSe, HSSe );
        // ### Sediments Layer ### //
      gr.WNode(ir,iaz,1) .SetAttributes( VpVs( 4.52, 2.61), 2.21, QSe, HSSe );
      gr.WNode(ir,iaz,1) .SetAttributes( VpVs( 6.20, 3.58), 2.80, QCrtmp, HSCrtmp );
        // ### Crust Layer ####### //
      gr.WNode(ir,iaz,2) .SetAttributes( VpVs( 6.24, 3.60), 2.82, QCrtmp, HSCrtmp );
      gr.WNode(ir,iaz,2) .SetAttributes( VpVs( 7.70, 4.44), 3.39, QMo, HSMo );
        // ### Moho Transition ### //
      gr.WNode(ir,iaz,3) .SetAttributes( VpVs( 8.00, 4.46), 3.40, QMa, HSMa );
        // # Mantle Layer (Moho to -80 km)
      gr.WNode(ir,iaz,4) .SetAttributes( VpVs( 8.040, 4.48), 3.50, QMa, HSMa );
      gr.WNode(ir,iaz,4) .SetAttributes( VpVs( 8.045, 4.49), 3.50, QMa, HSMa );
        // # Mantle Layer (-80 to -120 km)
      gr.WNode(ir,iaz,5) .SetAttributes( VpVs( 8.051, 4.50), 3.43, QMa, HSMa );
        // # Mantle Layer (-120 to -210 km)
      gr.WNode(ir,iaz,6) .SetAttributes( VpVs( 8.301, 4.52), 3.32, QMa, HSMa );
        // # Mantle Layer (-210 to -360 km)
      gr.WNode(ir,iaz,7) .SetAttributes( VpVs( 8.848, 4.78), 3.46, QMa, HSMa );

    }

  }

}
