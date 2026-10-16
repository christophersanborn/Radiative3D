# Multithreaded Simulation Handoff

## Branch state

The implementation is on `multithread-srf`. The branch is intended to remain
available for testing and review and is based on the multithreaded development
line.

## User controls

- Multithreading is enabled automatically.
- `--threads 1` forces a serial reference run.
- `--threads N` explicitly selects a positive worker count.
- `--seed N` selects a reproducible non-negative seed.
- If no seed is supplied, the generated effective seed is printed at startup.

Automatic selection is performance-core aware where supported. On macOS it
uses `hw.perflevel0.logicalcpu`, reserves one performance logical CPU, and
clamps to the requested phonon count. Other platforms use
`std::thread::hardware_concurrency()` with the same one-thread reservation
policy. Explicit `--threads` values remain authoritative.

## Performance finding

The main scaling bottleneck was unnecessary global report-mutex acquisition.
The Lop Nor benchmark uses `--reports=ALL_OFF`, but several hot-path report
methods locked the shared mutex before checking whether reporting was enabled.
Moving the lock inside the enabled-report branch restored strong scaling.

Measured on a macOS host with 24 logical CPUs, including 16 performance
logical CPUs, using 210,000,000 phonons:

| Mode | Workers | Seconds | Phonons/s | Speedup |
|---|---:|---:|---:|---:|
| Serial | 1 | 221.917 | 946,300 | 1.00x |
| Explicit | 12 | 27.623 | 7,602,360 | 8.03x |
| Automatic | 15 | 28.122 | 7,467,463 | 7.89x |

The explicit 12-worker run was slightly faster than automatic selection in
this sample; automatic selection is hardware-aware, not a promise to choose
the optimum count for every workload.

## Reproduction

Build and run the full regression checks:

```sh
make cleanall && make
tests/test_simulation_options.sh
tests/test_worker_policy.sh
tests/test_seed_reproducibility.sh
tests/test_diagnostic_rng.sh
tests/test_multithread_accumulation.sh
tests/test_multithread_output.sh
tests/test_lopnor_thread_equivalence.sh
```

Run the simulation-only large benchmark without Octave:

```sh
LOPNOR_BENCH_PHONONS=210000000 \
LOPNOR_BENCH_OUTPUT_DIR=/path/to/benchmark-output \
tests/benchmark_lopnor_threads.sh
```

The benchmark runs serial, explicit 12-worker, and automatic modes with the
same fixed-seed Lop Nor scenario. Detailed multi-worker report ordering is not
guaranteed, but trace formats and aggregate outputs remain compatible.
