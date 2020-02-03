# Radiative3D
Radiative transport in 3D Earth models

Radiative3D is a 3D radiative transport software tool being developed by Christopher Sanborn and the Solid Earth Geophysics Research Group at the University of Connecticut. Radiative3D can be used to produce synthetic waveforms, travel-time curves, or volumetric visualizations of energy propagation through three-dimensional Earth models. Radiative3D uses ray tracing to simulate propagation dynamics in large-scale structure, and uses a stochastic multiple scattering process to simulate the effects of statistically-described small-scale structure. Radiative3D simulates realistic source events described by moment tensor elements, allowing it to be used to simulate a variety of focal mechanisms, including explosions, double-couple earthquakes, CLVD's, etc.

#### Publications:

* [Modelling Lg blockage, GJI, 2018](https://academic.oup.com/gji/article/214/2/1426/5000173)
* [Simulations with Radiative3D, 2017](https://opencommons.uconn.edu/dissertations/1460/)
* [Combined effects of deterministic and statistical structure, GJI, 2017](https://academic.oup.com/gji/article/210/2/1143/3833065)

## What's New:

##### Janaury 2020:

* Added support for model architecture composed of concentric spherical shells, in which the elastic velocity profiles vary quadraticaly with the radial coordinate (_v<sub>P,S</sub> = a r^2 + c_).  The grid defines the velocities at the top and bottom of each spherical layer. This allows for whole-Earth models, efficiently defined in terms of a reasonably small number of layers.

* Raspbian Buster on the Raspberry Pi 4 is now my primary development and testing environment.  Radiative3D is fast, efficient, and runs on inexpensive hardware!  Of course, it still runs great on Intel and AMD-based workstations as well.

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
