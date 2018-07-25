{{Radiative3D}}
{{TOC limit|3}}
=== NAME ===
<blockquote>
[[Radiative3D]] - A program that uses [[radiative transport]] theory to model seismic energy propagation in 3D Earth models.
</blockquote>

=== SYNOPSIS ===
<blockquote>
When run directly, Radiative3D is invoked as follows:

: <tt>./main [OPTIONS]...</tt>

Often, however, Radiative3D runs are scripted, and a few example "do-scripts" are included as samples.  The scripts: (1) set up model and simulation parameters and assemble them as command line options, (2) start the simulation in a timestamped output directory, (3) record program output and other information into log files, and (4) run visualization scripts after run completion to produces figures or videos. One such example do-script is <tt>do-crustpinch.sh</tt>.  It is recommended to copy the script and modify to taste.  It is run in the following way:

: <tt>./do-crustpinch.sh [RUN_ID]</tt>

where ''<tt>RUN_ID</tt>'' is an optional alphanumeric identifier that will become part of the name of the output directory. 
</blockquote>

=== DESCRIPTION ===
<blockquote>
[[Radiative3D]] models energy propagation through both deterministic and statistical structure which are specified and characterized on a gridded dataset referred to as the "Earth model".  Arguments to '''Radiative3D''' control various aspects of the operation, including specification of the Earth model, specification of the broader "physics model", and the selection and control of various inputs and outputs.

==== Arguments ====
<blockquote>

<tt>-F, --frequency=''FREQ''</tt>
: Sets the frequency of all generated phonons.  Radiative3D operates at a single frequency.  To generate broadband results, do multiple data runs and combine the results.

<tt>-N, --num-phonons=''NUM_PHONONS''</tt>
: Number of phonons to emit from source before exiting.  Can take a multiplier suffix of <tt>K</tt>, <tt>M</tt>, or <tt>B</tt> for thousand, million, and billion, respectively.  E.g., use <tt>--num-phonons=10M</tt> for ten million phonons.

<tt>-T, --timetolive=''time''</tt>
: Sets the time to live for each phonon, or in other words the amount of sim-time to model.

<tt>-A, --toa-degree=''TOA_DEGREE''</tt>
: Sets the take-off angle degree, which determines the number of discrete [[take-off angle]]s that PhononSource objects (e.g. event sources and scattering sources) have to choose from.  Increasing this number increases the angular fineness of these objects.  "Degree," in this context, is an integer that represents the number of times the angular space is subdivided, and in the default (and at present only) method for discretizing take-off angles, the algorithm starts with 20 angles and subdivides by four with each iteration (degree), such that N_TOA = 20 * 4 ^ (degree).

===== ''Grid args:'' =====
<tt>--grid-file=''FILENAME''</tt> ''(not yet implemented)''
: This specifies that the model grid should be read from file ''FILENAME'' and used for the construction of the Earth model.  The grid file format is described in: ''[[Grid formats in Radiative3D]]''. (Note: at the current time, reading of grid files is not yet supported, and users should hard-code models. See the <tt>--grid-compiled</tt> option.)

<tt>--grid-compiled[=''INDEX'']</tt>
: Specifies that the model grid to be used is compiled-in and defined in the function <tt>Grid::ConstructGridManual()</tt> located in the user-modifiable source file <tt>user.cpp</tt>.  Defining a model in source code allows full access to the [[grid building API]] and can facilitate construction of quite complex models, in comparison to models defined via [[Grid formats in Radiative3D|grid files]] (c.f. <tt>--grid-file</tt> option).  It can also enable aspects of the model to be parameterized by a set of model arguments provided by the <tt>--model-args</tt> option in order to facilitate scripting of structural variations in a model.  The optional integer-valued ''INDEX'' parameter, if included, will be passed to <tt>Grid::ConstructGridManual()</tt> and can be used as a selection index if the user wishes to define multiple user-coded models and choose between them at runtime without need of a recompile.

<tt>--dump-grid</tt>
: If specified, a parsable plaintext dump of the grid will be written to <tt>stdout</tt> prior to simulation start.

===== ''Model args:'' =====
<tt>--model-args=''ARG_LIST''</tt> <br>
<tt>--model-compiled-args=''ARG_LIST''</tt> ''(deprecated form)''
: Some model grids allow parameterized values, enabled structural variations of the model to be easily scripted. This option allows a comma-separated list of double-precision values to be passed in on the command line and made available as parameters to these grids.

<tt>--range=''RANGE''</tt> <br>
<tt>--cylinder-range=''RANGE''</tt> ''(deprecated form)''
: When using a LAYERED model defined by a plumb-line grid, the model cells are horizontal layers separated by planar tilted interfaces.  The otherwise infinite lateral extent of these layers is truncated by imposing a cylindrical bounding radius centered on a vertical through the origin.  This option sets the bounding radius, effectively defining the distance range that is modeled.  This option only applies to LAYERED models and is ignored for TETRA models based on tetrahedral model cells.

<tt>--flatten</tt>  ''*** DEPRECATED ***''
: Apply an [[Earth-flattening transformation]] to depth coordinates and velocities when constructing the model.  This is intended to enhance ray-turning in layered Earth models as a proxy for the spherical curvature that is lacking from those models.  See [[Aki and Richards 1980 - Quantitative Seismology|Aki and Richards]] for more info on EFT's. (Note: the velocity transform is performed only for layered model types, but the depth transform occurs for both layered or tetrahedral based models, though it should be noted it is generally unnecessary for tetra models, as the model mesh can be constructed to follow real Earth curvature.)

: Deprecation note:  the <tt>--flatten</tt> option is deprecated, as coordinate mapping choices, including flattening, are now attributes of the Grid object, and as such are specified in the grid definition file, or via Grid class methods if hard-coding a grid.  If a mapping choice is specified in this way, then the <tt>--flatten</tt> option will have no effect and will be silently ignored.  This option may be removed in a future release.

===== ''Sim args:'' =====
<tt>--overridemfp=''P'',''S''</tt>
: If specified, use the provided mean free path values for scattering model-wide, despite what the calculated mean free paths would be from the otherwise specified heterogeneity parameters.  Scattering shapes (deflection profiles) are unaffected by this argument.

<tt>--nodeflect</tt>
: If specified, scattering deflections are "squashed," meaning that while scattering events still occur, they do not result in any deflection or other modifications (raytype or polarization changes) to the phonon's ray trajectory.  This is commonly used, sometimes in concert with <tt>--overridemfp</tt> for producing videos (which plot scattering events), to easily visualize "clean" (i.e. no scattering) wavefront propagation. Or in other words, this option changes scattering events into mere "checkpoint" events used to illustrate evolving wavefronts.

===== ''Event args:'' =====
<tt>-L, --source-loc=''X'',''Y'',''Z''</tt>
: Specifies the location of the event source in cartesian X,Y,Z coordinates.  In the future, other coordinate systems may be supported (such as Lat, Lon, and depth).

<tt>-E, --source=''SOURCE-ARGS''</tt>
: Determines the event source type via the provided source arguments.  The format of the source args is a keyword followed by a comma-sepated (no whitespace) list of arguments that are required by the given keyword.  The keywords are as follows:
:
: <tt>--source=EQ</tt>
:: Specifies a "basic earthquake" source.  The source type will be a simple double-couple with principle axes in the X-Y plane.  (Ie, energy radiates predominantly in the X-Y directions, with little energy radiated virtically from the source.)
: <tt>--source=EXPL</tt>
:: Specifies an "explosion" source, which radiates equal energy in all directions (isotropic energy release).
: <tt>--source=USGS,''rr'',''tt'',''pp'',''rt'',''rp'',''tp''</tt>
:: Specifies the full moment-tensor using the USGS standard of specifying elements in a reference with axes pointing UP, SOUTH, and EAST (NOTE: NEED TO FIGURE OUT WHAT THE DIRECTIONS ACTUALLY ARE...)  This method of specifying depends on a specified Earth-center and polar-axis being defined—these parameters are set by the <tt>--earth-frame</tt> argument.
: <tt>--source=SDR,''strike'',''dip'',''rake'',''[isofrac|isoangle'',''[moment]]''</tt>
:: Specifies a double-couple [[moment tensor]] with an optional isotropic component added in.  The double-couple orientation is given in terms of ''<tt>strike</tt>'', ''<tt>dip</tt>'', and ''<tt>rake</tt>'' angles, given in degrees.  The optional fifth argument, ''<tt>moment</tt>'', sets the magnitude of the moment tensor, which otherwise defaults to 1.0.  This parameter has no bearing on the actual simulation, but does affect the numbers displayed when the moment tensor is output to the console (the parameter has purely aesthetic utility).  The optional fourth argument, ''<tt>isofrac</tt>'', specifies the isotropic component of the moment tensor as a fraction of the total squared-magnitude of the moment tensor. Valid values are in the range [-1.0 to 1.0], with 0.0 being the default.  Positive values indicate an explosive component, whereas negative values indicate an implosive component. Zero means pure double-couple, or no isotropic component.  This argument can alternatively be specified as an angle, rather than a fraction, owing to the fact that [[isotropic]] and [[deviatoric]] moments form orthogonal subspaces of the total moment tensor vector space.  If this manner of specification is desired, the value ''<tt>isoangle</tt>'' should be given in the split range of [-90, -1) or (1, +90].  (The range between -1 and 1 is inaccessible because it will be interpreted as ''<tt>isofrac</tt>'' rather than ''<tt>isoangle</tt>'', but this is unlikely to be problematic in practice.)

===== ''Reporting / data output args:'' =====
<tt>--reports[=''REPORTS_FLAGS'']</tt>
: Enables real-time event reporting.  These are reports of various simulation events detailing the trajectories of individual phonons, such as generation at source, reflection off interfaces, scattering events, collection by seismometers, etc., written to an output stream as they are simulated.  The aggregate of these micro-reports can be used for, among other things, producing videos of energy propagation within the Earth model.  Multiple keywords can be provided as a comma-separated list in order to specify which event types are desired. Event reports are written to <tt>stdout</tt> unless directed to a file by the <tt>--report-file</tt> option.  The keywords are:
: <tt>--reports</tt>, <tt>--reports=ALL_ON</tt>
:: Turns ON all real-time reports. (If called without any keyword, <tt>ALL_ON</tt> is assumed.)
: <tt>--reports=ALL_OFF</tt>
:: Turns OFF all real-time reports. ''(DEFAULT)''
: <tt>--reports=GEN,SCT,REF,CEL,COL,LST,TMO</tt>
:: This shows all event types explicitly turned on.  They are, respectively: Generate, scatter, reflection, cell-to-cell transfer, collection (interaction with a surface where a seismometer may reside), phonon lost (departing boundary of model), and phonon timeout (phonon lifetime exceeds window of interest as determined by the <tt>--timetolive</tt> option).
: <tt>--reports=SCATTERS</tt>
:: This is a shortcut equivalent to <tt>--reports=GEN,SCT,REF</tt> to report generate, scatter, and reflection events.  Useful for visualizing bulk energy transport, e.g. in video form.

<tt>--report-file=''filename.dat''</tt>
: Specifies the name of the file that real-time event reports will be written to.  Without this option, reports go to <tt>stdout</tt>.  If option is set, reports go only to the file and will not be printed to <tt>stdout</tt>.

<tt>--output-dir=''path''</tt>
: Specifies a directory name in which in which all files created by Radiative3D will be created.  If this option is not specified, then all files are relative to the current working directory.

<tt>--mparams-outfile=''filename.octv''</tt>
: If this is specified, various model and mission parameters will be written to ''filename.octv'' in a format readable by GNU Octave for easy processing.

<tt>--summary=''SUMMARY_FLAGS''</tt>
: Enables or disables various post-simulation summary reports.  These include such things as the trace data from the various seismometers that were requested. Multiple keywords can be provided as a comma-separated list.  The keywords are:
: <tt>--summary=ON</tt>
:: Turns ON all summary reports
: <tt>--summary=OFF</tt>
:: Turns OFF all summary reports

===== ''Seismometer args:'' =====
<tt>--binsize=''number''</tt>
: When binning energy for seismogram output, this is the width of each bin in time units (typically seconds).  (Mutually exclusive with <tt>--binspercycle</tt>.)

<tt>--binspercycle=''number''</tt>
: An alternate way to specify bin size.  This lets you specify the number of time bins per cycle at the selected frequency.  Can be fractional.  E.g. 1.5 time bins per cycle means three bins will cover two cycles worth of time.  (Mutually exclusive with <tt>--binsize</tt>.)

<tt>--seismometer=''SEIS_ARGS''</tt>
: Adds a seismometer according to ''SEIS_ARGS''.  This flag can be passed multiple times.

<tt>--seis-array=''AZ,R0,gap,gather,#seismometers''</tt>
: Adds an array of seismometers according to ''SEIS_ARRAY_ARGS''.  This flag can be passed multiple times.

<tt>--seis-p2p=''X1,Y1,Z1,X2,Y2,Z2,R0,gather,num_seismometers''</tt><br>
<tt>--seis-p2p=''X1,Y1,Z1,X2,Y2,Z2,R0,gather1,gather2,num_seismometers''</tt><br>
<tt>--seis-p2pw=''X1,Y1,Z1,X2,Y2,Z2,R0,gather,num_seismometers''</tt><br>
<tt>--seis-p2pw=''X1,Y1,Z1,X2,Y2,Z2,R0,gather1,gather2,num_seismometers''</tt><br>
: Adds a "point-to-point" array along the line (X1,Y1,Z1) to (X2,Y2,Z2).  If two gather radii are specified, then the gather radius increases linearly from ''gather1'' at the beginning point to ''gather2'' at the end point.  If the option is specified as <tt>--seis-p2pw</tt>, then the gather radii are specified in elastic wavelengths, which are a function of elastic velocity at the seismometer's emplacement location at the specified frequency (and will result in differing gather radii for P and S wave modes).  Otherwise, the gather radii are specified in the model's chosen length unit, typically kilometers.  (Note however that, at present, the <tt>--seis-p2p</tt> option, if specified with only a single gather radius, still assumes this radius is in wavelengths.  This is for backwords compatibility with old scripts, but this usage is deprecated and will be changed in a future version.)

: This option can be specified multiple times.

</blockquote>

==== Earth models ====
Earth models are constructed as a mesh of nodal locations at which known or hypothesized material properties are specified. Based on these nodal locations, Radiative3D divides space into a set of model cells inside which material properties are interpolated based on the values specified at the nodes.  The material properties specified include the elastic properties of ''Vp'', ''Vs'', ''density'', and ''Q'', and four additional properties defining the [[heterogeneity spectrum]] of the material: ''nu'', ''epsilon'', ''a'', and ''kappa''.

Radiative3D currently supports two distinct modeling and meshing approaches: stacked layers of spatially uniform elastic properties, or a warped cartesian grid (WCG) of nodes in which material properties vary linearly based on a tetrahedral tessellation of the space using the node locations as the corners of the tetrahedral model cells.  The former approach is more akin to "1-D" models, and the latter approach fully utilizes the 3-D modeling capabilities of Radiative3D.

Currently, Radiative3D lacks file importers for reading Earth models into memory. To construct Earth models, the model structure must be hard coded using our model-building API.  Changes to the model structure thus require a recompile.  These "user-coded" models are coded in the source file called <tt>user.cpp</tt>, and inspection of that file will reveal several example models to choose from, and explain the model building API.

==== Dimensions and units ====
For user interaction (input and output), [[Radiative3D]] is generally agnostic about choice of physical units.  The only requirements are consistency of choice, and an understanding that the units you use for input will determine the units you get for output.  Specific details are as follows:

===== ''Distance and time:'' =====
[[Radiative3D]] is designed for Earth-scale models and temporal resolutions of ~1 Hz or greater.  As such there is a general presumption that length quantities will be in kilometers and time quantities will be in seconds, but there is nothing specifically constraining this.  What you use for input is what you will get for output.  E.g., if you use kilometers for your model dimensions, and want temporal output in seconds, then make sure you use kilometers/second as your velocity units.

===== ''Frequency:'' =====
For input/output, [[Radiative3D]] uses ordinary frequency measurements (i.e., cycles-per-time, not radians-per-time).  Frequency measures will share the time dimension with velocities.  I.e, if you are inputting velocities in length-per-second, then input frequency using cycles/second, or Hertz.

===== ''Density:'' =====
For most use cases, [[Radiative3D]] sets no restriction on choice of units for densities.  Where densities are used in computation, only ratios come into play.  Thus it is user choice what units to use.

The only exception to this is when Earth models are defined in terms of shear and bulk moduli (a feature not yet supported), instead of P and S velocities.  In this case, the choice of density units must be in agreement with the choice of moduli units to produce the desired length and time units.

===== ''Attenuation Q:'' =====
[[Radiative3D]] interprets Q values as the inverse of the per-radian fractional rate-of-energy-loss due to intrinsic (e.g. thermal) attenuation. For instance, for a wavefront in a uniform-Q medium, where the energy decays with time ''t'' as:

<math>E(t) = E_0\exp\left[-\frac{\omega t}{Q}\right]</math>,

inverse Q is equivalent to:

<math>Q^{-1} \equiv \frac{-1}{E} \; \frac{\mathrm{d}}{\mathrm{d}(\omega t)} E(t)</math>.

Or, thinking of Q as an angular time constant, Q represents the number of radians ''(ω*t)'' of wavefront evolution that result in a 1/''e'' energy attenuation of the wave.

Thinking in terms of ordinary cycles, instead of radians, Q is the number of cycles ''(f*t)'' of wavefront evolution resulting in an exp[-2π] energy attenuation (about a 99.8% reduction) of the wave.

</blockquote>

=== EXAMPLES ===
<blockquote>
==== Crust pinch: envelopes ====
An envelope and travel-time curve simulation in a crustal pinch model can be run via the example script called <tt>do-crustpinch.sh</tt>.  From the Radiative3D install directory, run:
: <tt>./do-crustpinch.sh test01</tt>

The script will create an output directory at <tt>./data/YYYYMMDD-HHMMSS-test01-R3D/</tt> and will run Radiative3D inside that directory, collecting all output into various log files and data files.  Upon completion of the simulation, it will process the output and produce a series of figure files depicting envelope traces and travel time curves.  The output directory name is a combination of a time stamp and the run ID provided as the first argument to the script ("test01" in this case).

Note that to get good quality envelopes, tens of millions of phonons will be propagated, and you can expect the simulation to run for a few hours.  For this reason, it may be prudent to run the simulation inside a GNU screen session, especially if you are running it on a remote computer.

==== Crust pinch: video ====
A video simulation in a crustal pinch model can be run via the example script called <tt>do-crustpinch-vids.sh</tt>.  From the Radiative3D install directory, run:
: <tt>./do-crustpinch-vids.sh test02</tt>

The script will create an output directory at <tt>./data/YYYYMMDD-HHMMSS-test02-R3VID/</tt> and will run Radiative3D inside that directory, collecting all output into various log files and data files.  Upon completion of the simulation, it will process the output and produce a series of still frames depicting phonon propagation as a function of time, and will concatenate those frames onto two videos, one depicting a viewpoint from above, and the other depicting a side-elevation viewpoint.

To get good quality images, a few tens of thousands of phonons will be propagated.  You can expect the simulation to run for about ten to twenty minutes, and the video generation to take another twenty minutes after that.  Videos are generally much less computationally intensive to produce than envelopes and travel-time curves are, because for videos every phonon is useful, whereas for envelopes only those phonons that interact with a given seismometer are useful in producing the envelope trace.

</blockquote>

=== INSTALLATION ===
<blockquote>
==== Via Subversion ====
Source code for Radiative3D can be obtained via Subversion.  When a version 1.0 release becomes available, the following command will retrieve the source code from our repositories:

: <tt>svn co https://rainbow.phys.uconn.edu/svn/Radiative3D/tags/Radiative3D-v1.0.0/ [destdir]</tt>

This will download the v1.0.0 release and place the contents in destination directory <tt>destdir</tt>, or else it will create directory <tt>Radiative3D-v1.0.0</tt> if the optional argument <tt>[destdir]</tt> is omitted.

''NOTE: At the present time, we are not ready for a version 1.0 release.  If you want access to the code as it stands right now, please contact the authors, and we will provide you the repository address of the latest development version, along with a long list of "caveat emptor" warnings about using alpha versions of a not-yet-validated codebase.  We do welcome your use of the code, however, as well as any feedback you may wish to provide.''

Radiative3D is written in C++ and compilation is automated with GNU Make.  Compiling is simply a matter of typing "<tt>make</tt>" on the command line.

==== Requirements ====
Radiative3D itself has no requirements other than a reasonably recent C++ compiler.  It is tested to compile and run on the following platforms: Fedora Linux, Ubuntu Linux, and Macintosh OS X.  While we strive for platform-independence in the development of Radiative3D, the accompanying scripts that automate Radiative3D simulation runs and produce visualizations from its output may have cross-platform issues.  Linux is the ''recommended'' platform for these scripts. The following packages are required to run the automation and visualization scripts in unmodified form:

* '''Bash''' interpreter.
**This is standard on Linux and OS X systems, but not on Windows.  The scripts are tested under Linux but may use some non-standard commands under OS X.  Cygwin may provide a usable Bash interpreter for Windows.
* '''GNU Octave''' — for figure generation
**GNU Octave is an open-source implementation of the Matlab programming language, and is available for Linux, OS X, and Windows.  Octave is usually available through the package management utilities on modern Linux systems.  On OS X it can be installed via Homebrew, a third-party open source package management system.
* '''Gnuplot'''
** This is a standard plotting utility on Unix systems and provides the back end to Octave's plotting system. It is available through the system package management tools or via Homebrew on Mac.
* '''LaTeX''' (recommended)
** The visualization scripts use <tt>pdfcrop</tt> to get rid of whitespace in pdf figure output. This tool is typically part of a LaTeX distribution. On Linux, you may need to install package "texlive-extra-utils" to get it. On OS X it should be included in the full installation of MacTex. 
* '''epstool''' and '''fig2dev''' (optional)
** These aren't actually needed by the plotting scripts, but Octave will complain about their absence during plotting if they are not installed.  ''epstool'' is its own package on the Linux systems we have tested, and fig2dev is part of the ''transfig'' package.  Both packages are available through Homebrew on OS X.
* '''ffmpeg''' (optional; for video generation scripts)
** ffmpeg is an open-source video encoding tool that we use to concatenate still frames into videos.  It is available for Linux, OS X, and Windows.

</blockquote>

=== BUGS ===
<blockquote>
==== Limitations of algorithm ====
The code employs [[radiative transport]] for modeling the propagation of energy within three-dimensional solid bodies.  Radiative transport is an extension of [[ray theory]], which can be thought of as the limiting case of full-wave treatment for infinite frequencies.  Since we use it to model finite-frequency physics, ray theory, (and hence radiative transport), may be inappropriate and give erroneous results for certain Earth-model configurations.  In particular, radiative transport may have trouble with: steep velocity gradients or model features with sizes on the order of a wavelength or smaller, free-surface effects, sharp interface effects, and the propagation of free-surface or interface waves (diffracted phases), or other instances in which waves would propagate along a 2D manifold rather than a 3D space.

The code models phonons at a single frequency per run.  To produce broadband synthetics, multiple simulation runs must be conducted, and the results summed.  Dispersion is not modeled by the code, however it can be explicitly input by using custom Earth models for each desired frequency.

==== Validation status ====
Radiative3D is in need of validation.  Until it is validated, it should be assumed that results produced with Radiative3D not only may contain, but ''do'' contain, errors and inaccuracies in multiple areas, and that extreme caution should be used in drawing any definitive conclusions from the results.  Validation is an important aspect of code development, and it would be foolish to assume that logical and/or implementation errors will not be found — and corrected — during this stage.

Validation would consist of an iterative process of code review and code testing, in which models of varying degrees of complexity are simulated and results compared against expectations.  In some cases, test models will be simple enough that expected results can be computed analytically. In other cases, the test models will be complex, but can be simulated by other methods, and expectations arrived at in that manner.

It should also be noted that validation, once complete, will not be an assertion that the results produced by Radiative3D are physically correct, but only that they can be reasonably expected to portray an accurate picture of the radiative transport method applied to an Earth model.  [[Radiative transport]] itself, and any limitations of the supplied Earth model, both impose constraints on the physical realism that can be expected from the resultant simulation output.

As regards the types of issues that may be revealed during validation, our greatest concern is with the apportionment of energy between the three components of motion for seismic envelope output. (E.g., between radial, transverse, and vertical components.)  There are known (but correctable) issues with the current handling of S-wave polarization for curved ray paths, which will affect channel apportionment.  Additionally, the authors would like to do a close review of the energy apportionment when rays reflect off the surface.  We have a greater degree of confidence in the timing of phase arrivals and in the summed (three-component) energy arriving in any given time window, though these must also be subject to the validation process.

==== Portability ====
Radiative3D should compile and run under multiple platforms (Linux, OS X, Windows...), however at present the scripts and visualization routines are not platform-independent. The automation scripts are in Bash, and use some commands that are not standard on OS X, which means Linux is the only 100% supported platform.  The graphics are produced using GNU Octave, which, while nominally platform-independent, it has been our empirical observation that graphical results differ from one platform to another, and even from one Linux distribution to another.  The graphics subsystem in Octave seems to be in a continual state of flux.

For this reason, it is our belief that reimplementing our graphics and automation scripts in a more consistently cross-platform language, such as Python, would be a very good idea.  For now, though, this is on our to-do list, but is not a top priority.  

</blockquote>

=== AUTHORS ===
<blockquote>
[[Radiative3D]] is written by [[User:ChristopherSanborn|Christopher J. Sanborn]] and the [[Solid Earth Geophysics Research Group]] at UConn, including significant contributions by [[User:Swalsh|Steven Walsh]].  The code implements algorithms and takes inspiration from two major forerunner projects: [[Raytrace3d]] by William Menke and [[PSPhonon]] by Shearer and Earle.
</blockquote>

=== COPYRIGHT ===
<blockquote>
The authors have yet to decide on an official copyright to use.  For now, if you have been given a copy of the code by one of the authors, you can assume you have permission to use it and publish results from it.  If you wish to share the code with others, please contact the authors first for permission.  The authors retain ALL RIGHTS to the code and disclaim ALL LIABILITIES for use of the code.
</blockquote>
