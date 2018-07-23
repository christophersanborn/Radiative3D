# seismap.m
#
# Plot a circle representing the gather radii for each seismometer in
# a list of seis files.
#
function seismap(seisglob="seis_*.octv",
                 paramfile="out_mparams.octv",
                 window=[-600,600,-350,600],
                 labelthese = [16 32],
                 labels = {"WUS", "MAK", "LOP"})

  # Define a unit circle, for plotting purposes (see drawcirc
  # function)
  global UnitCircX;
  global UnitCircY;
  UnitCircPh=linspace(0,2*pi,121);#181);
  UnitCircX=cos(UnitCircPh);
  UnitCircY=sin(UnitCircPh);

  # Get event location:
  MParams = load(paramfile);
  ELocX = MParams.EventSourceLocUCS(1); # User coords, ...Hoping they're XYZ.
  ELocY = MParams.EventSourceLocUCS(2); # Will need to re-write if not.

  # Initialize plot
  figinit(5.0, 4.00,                # Clear and initialize 
       "paperunits", "inches",      # a 5.0" x 3.75" figure
       "visible", "off");
  axes("fontname", "Courier",       # Create and setup an axes object
       #"fontsize", fontaxes,       #
       "nextplot", "add");  # equiv to "hold on"

  axis(window);
  title("Seismometer placement and gather radii", "fontname", "Helvetica");
  xlabel("Distance (km)", "fontname", "Helvetica");

  # Plot range circles:
  rangecirccolor=[0.8,0.8,0.8];
  rangecircwidth=1.5;
  for R = (1:9)*100
    drawcirc(ELocX,ELocY,R,rangecircwidth,rangecirccolor);
    h=text(ELocX,ELocY+R,sprintf("%g km",R),"horizontalalignment","left",
           "verticalalignment","bottom", "color", rangecirccolor,
           "clipping", "on");  # clipping on apparently broken...
    if ((ELocY+R)>window(4))
       delete(h);              # mimic what clipping is *supposed* to do
    end
  end

  # Cardinal Directions:
  directionwidth=1.5;
  directionline=":";
  directioncolor=rangecirccolor;
  line([ELocX-920,ELocX+920],[ELocY,ELocY],"linewidth",directionwidth,
                                           "linestyle",directionline,
                                           "color", directioncolor);
  line([ELocX,ELocX],[ELocY-920,ELocY+920],"linewidth",directionwidth,
                                           "linestyle",directionline,
                                           "color", directioncolor);
  line(sqrt(0.5)*[-920,+920]+ELocX,sqrt(0.5)*[+920,-920]+ELocY,
                                           "linewidth",directionwidth,
                                           "linestyle",directionline,
                                           "color", directioncolor);
  line(sqrt(0.5)*[-920,+920]+ELocX,sqrt(0.5)*[-920,+920]+ELocY,
                                           "linewidth",directionwidth,
                                           "linestyle",directionline,
                                           "color", directioncolor);
  i=0;
  for seisfile = {glob(seisglob){:}}
    i++;
    SEIS = load(seisfile{1});
    if (isfield(SEIS,"SEIS")) SEIS=SEIS.SEIS; end # Catch old format
    SLocX = SEIS.Location(1);
    SLocY = SEIS.Location(2);
    SRadS = SEIS.GatherRadius(2,2);
    SRadP = SEIS.GatherRadius(1,2);
    linewidth = 2.5;
    colorS = [0.0,0.0,0.6];
    colorP = [0.6,0.0,0.0];
    HSS = drawcirc(SLocX,SLocY,SRadS,linewidth,colorS,4);
    HSP = drawcirc(SLocX,SLocY,SRadP,linewidth,colorP,4);
    j = find(labelthese==i,1);
    if (length(j)==1)
       Ltxt = labels(j);
       text(SLocX,SLocY+max(SRadP,SRadS),Ltxt,"horizontalalignment","center",
            "verticalalignment","bottom", "fontweight", "bold");
    end
  end
  if (length(labels)>length(labelthese))
    text(ELocX,ELocY,labels(end),"horizontalalignment","left",
         "verticalalignment","top", "fontweight", "bold");
  end

  # Legend
  h=legend([HSS,HSP],"S Gather", "P Gather", "location", "northwest");
  set(h,"fontname", "Courier");

end

###
# function: drawcirc()
# Draw a circle of radius R centered at X,Y.
# Assumptions: Plot figure and axes are already set up and "hold" is
# on, and that UnitCircX,Y are already initialized.
function H = drawcirc(X, Y, R, linewidth=1, linecolor=[0,0,0], skip=1)

  global UnitCircX;  # If desire to use skip=2, ensure 
  global UnitCircY;  # these initialized with odd number of elements
  H = plot (R*UnitCircX(1:skip:end)+X,
            R*UnitCircY(1:skip:end)+Y, 
            "linewidth", linewidth,  "color", linecolor);

end
