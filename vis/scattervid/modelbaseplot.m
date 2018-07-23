## modelbaseplot.m
##
## This one forked from vis/seisplot/modelplot.m and modified to work
## on an existing axes, used as a backdrop for elevation scatter
## plotting.
##
## Given grid input, plot a model cross-section indicating velocity
## and/or other attribute profiles.
##
## INPUT:
##
## GRIDFILE FORMAT:
##
##   Well, we will want to read the "real" grid file format at some
##   point, but for now, assume this file is the grid dump written to
##   std-out by Radiative3D - just the lines describing actual nodes
##   (i.e. the ones beginning with with an index tripple),
##   pre-processed to replaces the stars (indicating no attrbute data)
##   with 0's so that the load command won't choke.  This
##   pre-processing can be achieved with \\ cat griddump-raw.asc | sed
##   's/\*\*\*/000/g' > griddump.asc //.
##
##   Note that for now we are only dealing with "cylinder" grid
##   models.
##
##
function modelbaseplot (
             gridfile,           # name of file with grid data
             pair = [1 2],       # which depth plumbs to use (ix values)
             mparams = "out_mparams.octv", # contains source params
             seispat = "",   # seis-file pattern (e.g. "seis_%03d.octv")
             seisrange = [], # which ones to plot (eg. [160,169:10:319])
             txtkey = {"LOP","MAK","WUS"}  # Station names (corresponds
                                           # to elements of pair).
           )

  GG_Raw = load(gridfile);
  if (exist(mparams,"file")==2)
    MPAR = load(mparams);
  end


  ibase = 0;                            # (Assume zero, but might not be.)
  NX = (1-ibase)+max(GG_Raw,[],1)(1);   # Grid dimensions
  NY = (1-ibase)+max(GG_Raw,[],1)(2);   #
  NZ = (1-ibase)+max(GG_Raw,[],1)(3);   #

  NAttrDefs = size(GG_Raw,1);

  GG = zeros(NX, NY, NZ, 12);           # 12 Attributes (XYZ, VpVs,
                                        # Dens, QsQp, NEAK) in 3 indices.

  for (j=1:NAttrDefs)
      ix = (1-ibase)+GG_Raw(j,1);
      iy = (1-ibase)+GG_Raw(j,2);
      iz = (1-ibase)+GG_Raw(j,3);
      for (k=(1:12))
        GG(ix,iy,iz,k) = GG_Raw(j,3+k);
      end
  end


  ##
  ## Get depths and lateral X's in local coords
  ##

  XY1 = [GG(pair(1),1,1,1) GG(pair(1),1,1,2)];
  XY2 = [GG(pair(2),1,1,1) GG(pair(2),1,1,2)];
  XYdiff = XY2-XY1;
  lateraldist = sqrt(XYdiff*XYdiff');
  XYdiffnorm = XYdiff/lateraldist;
  X1 = 0;
  X2 = lateraldist;
  Z1 = zeros(NZ,1);
  Z2 = zeros(NZ,1);
  VP = zeros(NZ,1);
  VS = zeros(NZ,1);
  ScEps = zeros(NZ,1);

  for (j=1:NZ)
      Z1(j) = GG(pair(1),1,j,3);
      Z2(j) = GG(pair(2),1,j,3);
      VP(j) = GG(1,1,j,4);
      VS(j) = GG(1,1,j,5);
      ScEps(j) = GG(1,1,j,10);
  end
  Z0 = Z1;  # Assumes Z1 @ x=0; more complex if this not the case
  Zslope = (Z2-Z1)/lateraldist;

  ##
  ## Viewport:
  ##

  viewwidth = 1600;                     # Geographical units (km)
  viewdepth = 90;                       # Will add 15% of "above sea level"
  viewoverrun = viewwidth-lateraldist;
  XL = -viewoverrun/2;                  # X-coord of left edge of viewport
  XR = lateraldist+viewoverrun/2;       #    ''   '' right ''
  viewport = [XL, XR, -viewdepth, 0.15*viewdepth];
  # override left right: (we're not setting up axis anymore)
  XLIM=xlim();
  XL=XLIM(1);
  XR=XLIM(2);

  # Extrapolate Z's:

  ZL = Z0 + XL*Zslope;
  ZR = Z0 + XR*Zslope;

  ##
  ## Setup Plot Space:
  ##

  figwidth  = 5.0;
  figheight = 3.75;
  linethin   = 0.5;
  linethick1 = 3.75;
  fonttitle  = 9.5;
  fontxylabel = 9.0;
  fontaxes = 9.0;
  fontinset  = 9.0;


  ##
  ## Start Plotting:
  ##

  ## Plot regions (in reverse order so sedi layer shows)
  for (j=(NZ-1):-1:1)
      xll = XL;
      xul = XL;
      xur = XR;
      xlr = XR;
      zll = ZL(j+1);
      zul = ZL(j);
      zur = ZR(j);
      zlr = ZR(j+1);
      c = vpcolor(VP(j));
      h = fill([xll xul xur xlr xll],[zll zul zur zlr zll],c);
      set(h, "linewidth",0, "edgecolor",c);
  end

  ## Mark key locations:
  line([X1,X1],[Z1(1),Z1(NZ)], "linestyle", ":");
  line([X2,X2],[Z2(1),Z2(NZ)], "linestyle", ":");
  text(X1,Z1(1)+2, txtkey{pair(1)}, "horizontalalignment", "center",
                            "verticalalignment", "bottom");
  text(X2,Z2(1)+2, txtkey{pair(2)}, "horizontalalignment", "center",
                            "verticalalignment", "bottom");


  ## Plot layer lines:
  for (j=1:NZ)
      line([X1,X2], [Z1(j),Z2(j)]);
  end


  ## Plot V curves:
#  plotvcurve(XL,0.10*viewwidth,ZL,Zslope,ScEps,0.08,[1.0,0.5,0.0],3.0,"eps");
#  plotvcurve(XL,0.15*viewwidth,ZL,Zslope,VS,9.0,[0.0,0.0,0.8],3.0,"vs");
#  plotvcurve(XL,0.15*viewwidth,ZL,Zslope,VP,9.0,[0.8,0.0,0.0],3.0,"vp");

  ## Plot seismometers:
  for j=seisrange
    seisname=sprintf(seispat,j);
    if (exist(seisname,"file")==2)
      SEIS=load(seisname);
      if (isfield(SEIS,"SEIS")) SEIS=SEIS.SEIS; end # Catch old format
      SEISXY = SEIS.Location(1:2);
      SEISz  = SEIS.Location(3);
      SEISx  = (SEISXY-XY1)*XYdiffnorm';
      drawcirc(SEISx,SEISz,0.2,6,[0,0.5,0]);
    end
  end

  ## Plot event source
  if (isfield(MPAR,"EventSourceLoc"))
     ELXY = MPAR.EventSourceLoc(1:2);
     ELx  = (ELXY-XY1)*XYdiffnorm';
     ELz  = MPAR.EventSourceLoc(3);
     drawcirc(ELx,ELz,0.75,8,[1 0 0]);
  end
  
  
  ## Titles and labels:
#  titletxt=sprintf("Cross Section from %s to %s", txtkey{1}, txtkey{2});
#  title(titletxt);
#  xlabel(sprintf("Range from %s (km)", txtkey{1}));
#  ylabel("Depth (km)");


end


###
# function: plotvcurve()
# Plot a staistep "velocity" curve into current axes.
# Note: pure stairstep, no gradients.
# Tweaks Z-depths so velocity jumps appear at actual interface.
function H = plotvcurve(X0, Xfull, ZZin, Zslope, VVin, Vfull, color, thick=4.0, desc="")

  NZin = length(ZZin);
  range_z = floor(1.5:0.5:(0.5+NZin));
  range_v = floor(1:0.5:NZin);
  ZZ = ZZin(range_z);
  ZS = Zslope(range_z);
  VV = Xfull*VVin(range_v)/Vfull;

  XX = VV+X0;
  NX = length(XX);
  XXflip=XX([1,(2:NX)-(-1).^(1:(NX-1))]);
  XXmid=(XX+XXflip)*0.5;
  Zdelta=ZS.*(XXmid-X0);

  H = line(XX,ZZ+Zdelta,"color",color,"linewidth",thick);

  Zbot=axis()(3);
  Xbot=XX(find(ZZ<=Zbot,1));
  text(Xbot, Zbot, desc, 
       "fontsize", 5, "fontweight", "bold",
       "horizontalalignment","left", 
       "verticalalignment", "bottom");

end


###
# function: drawcirc()
# Draw a circle of radius R centered at X,Y.
# Assumptions: Plot figure and axes are already set up and "hold" is
# on.
function H = drawcirc(X, Y, R, linewidth=1, linecolor=[0,0,0])

  persistent UnitCircX; 
  persistent UnitCircY;
  persistent binit = false;

  if (!binit) 
    UnitCircPh=linspace(0,2*pi,33);
    UnitCircX=cos(UnitCircPh);
    UnitCircY=sin(UnitCircPh);
    binit=true;
  end

  viewport=axis();
  asp_x = (viewport(2)-viewport(1))*3;
  asp_y = viewport(4)-viewport(3)*4;
  aspect = asp_x/asp_y;  # quick-and-dirty, could be better
                         # (paperposition and position properties also 
                         # have bearing on this.)

  H = plot (aspect*R*UnitCircX+X, R*UnitCircY+Y, "linewidth", linewidth, 
                                          "color", linecolor);

end

###
#
function color = vpcolor(vp)

  cfun = [2.00   0.5, 0.3, 0.0;
          4.99   0.8, 0.5, 0.0;
          5.00   0.8, 0.6, 0.4;
          7.99   0.9, 0.8, 0.7;
          8.00   0.7, 0.4, 0.2;
          8.30   0.9, 0.6, 0.0;
          9.20   1.0, 0.8, 0.0;
          9.60   1.0, 1.0, 0.0
         ];

  color = [0,0,0]; # error color (shouldn't happen)
  if (vp <= cfun(1,1))
    color = cfun(1,2:4);
  elseif (vp <= cfun(end,1))
    for (j=1:(size(cfun,1)-1))
        if (vp <= cfun(j+1,1))
           frac=(vp-cfun(j,1))/(cfun(j+1)-cfun(j));
           color = frac * cfun(j+1,2:4) + ...
                   (1-frac) * cfun(j,2:4);
           break;
        end
    end
  else
    color = cfun(end,2:4);
  end

end
