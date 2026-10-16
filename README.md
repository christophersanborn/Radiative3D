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

## Multithreaded simulation performance

Radiative3D supports automatic phonon worker selection and an explicit
`--threads N` override. Use `--threads 1` for a serial reference run. The
simulation-only Lop Nor benchmark can be run with:

```sh
LOPNOR_BENCH_PHONONS=210000000 \
LOPNOR_BENCH_OUTPUT_DIR=/path/to/benchmark-output \
tests/benchmark_lopnor_threads.sh
```

This benchmark does not invoke Octave or generate figures. It runs the same
fixed-seed Lop Nor scenario with one worker, 12 workers, and automatic worker
selection. Results below were measured on:

| Date | Host | Compiler | Workload |
|---|---|---|---|
| 2026-10-16 | macOS host, 24 logical CPUs | Apple Clang/G++ C++11 build | 210,000,000 phonons |

| Mode | Workers | Phonons | Seconds | Phonons/s | Speedup |
|---|---:|---:|---:|---:|---:|
| Serial | 1 | 210,000,000 | 221.917 | 946,300 | 1.00x |
| Explicit | 12 | 210,000,000 | 27.623 | 7,602,360 | 8.03x |
| Automatic | 15 | 210,000,000 | 28.122 | 7,467,463 | 7.89x |

The benchmark was rerun after removing unnecessary report-mutex acquisition
when reporting is disabled. The original implementation serialized these
frequent report calls even with `--reports=ALL_OFF`; the corrected results
show the expected strong scaling up to the useful core range of this host.
Automatic selection uses the 16 performance logical CPUs reported by macOS and
reserves one; the explicit 12-worker run was slightly faster in this sample.

There are also supporting scripts (e.g. `do-crustpinch.sh`) to help with managing command line options and organizing the various output files and post-processing of data.  The "do-scripts" are writtin in BASH and may depend on the installation of additional command line tools.  (See Manual Page.)
