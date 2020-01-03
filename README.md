# Radiative3D
Radiative transport in 3D Earth models

Radiative3D is a 3D radiative transport software tool being developed by Christopher Sanborn and the Solid Earth Geophysics Research Group at the University of Connecticut. Radiative3D can be used to produce synthetic waveforms, travel-time curves, or volumetric visualizations of energy propagation through three-dimensional Earth models. Radiative3D uses ray tracing to simulate propagation dynamics in large-scale structure, and uses a stochastic multiple scattering process to simulate the effects of statistically-described small-scale structure. Radiative3D simulates realistic source events described by moment tensor elements, allowing it to be used to simulate a variety of focal mechanisms, including explosions, double-couple earthquakes, CLVD's, etc.

## Build Process

Radiative3D builds with GCC on MacOS (OS X), Linux, and Raspbian.  (And perhaps also Windows.)

```
$ git clone https://github.com/christophersanborn/Radiative3D.git
$ cd Radiative3D
$ make
```

Results in a binary named `main`.

Run with:

```
$ ./main [args]
```

User Manual here: [Radiative3D Manual Page](doc/MANUAL.md)

There are also supporting scripts (e.g. `do-crustpinch.sh`) to help with managing command line options and organizing the various output files and post-processing of data.  The "do-scripts" are writtin in BASH and may depend on the installation of additional command line tools.  (See Manual Page.)
