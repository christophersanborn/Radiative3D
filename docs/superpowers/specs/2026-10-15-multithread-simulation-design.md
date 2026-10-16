# Multithreaded Simulation Design

## Purpose

Enable Radiative3D to simulate more phonons per unit time and produce more
robust synthetic envelopes by safely parallelizing independent phonon
propagation tasks. Multithreading should be enabled automatically, while
remaining controllable for debugging, benchmarking, and reproducibility.

## Goals

- Select a sensible worker count automatically from available hardware and
  requested phonon count.
- Provide `--threads N`, including `--threads 1` for a single-thread reference
  run.
- Provide `--seed N` for reproducible random streams while retaining automatic
  seeding when no seed is supplied.
- Preserve existing propagation algorithms, envelope normalization, trace
  formats, and aggregate diagnostics.
- Avoid locks in the phonon-propagation hot path wherever practical.
- Ensure worker failures are reported on the main thread and do not produce a
  false successful-completion result.

## Non-goals

- No changes to the physical propagation algorithms or Earth-model geometry.
- No process-level distributed execution.
- No changes to the public seismometer trace file format.
- No broad redesign of the existing global/static architecture beyond state
  required for safe simulation workers.

## Architecture

`Model::RunSimulation()` will coordinate a dynamically sized worker pool.
Workers will claim phonon IDs from shared atomic work state, generate and
propagate one phonon per claimed ID, and repeat until no work remains.

The default worker count will be based on the available performance CPUs when
the platform exposes that information, the number of requested phonons, and a
conservative oversubscription policy. On macOS, automatic selection uses
`hw.perflevel0.logicalcpu` and reserves one performance logical CPU when more
than two are available. Other platforms fall back to
`std::thread::hardware_concurrency()` and reserve one hardware thread when
more than two are available. The effective count will never exceed the
requested phonon count. An explicit `--threads N` setting will override
automatic selection; positive values are valid, including values that
oversubscribe hardware.

Workers will not write progress output directly. The coordinator will report
progress from shared completion state at controlled intervals to prevent
interleaving and excessive synchronization.

Read-only model geometry and material data will be shared. Mutable state will
be isolated per worker where possible:

- random engine and random context;
- event-source generation state;
- seismometer accumulation bins;
- lost, timeout, and invalid counters;
- invalid-diagnostic bitfields.

After all workers have joined, the coordinator will merge worker results in a
fixed order and emit the existing final summaries and traces once.

## Randomness and reproducibility

The simulation path will stop relying on process-global `rand()`/`srand()`
state. A base seed will be accepted through `--seed N`; when omitted, it will
be generated from a high-resolution clock and process entropy. The effective
seed will be printed at simulation startup.

Each phonon will receive a stable random stream derived from the base seed and
phonon ID, preventing thread scheduling from changing the random inputs for a
given phonon. `--threads 1` will therefore provide a stable reference sequence
for a given seed. Aggregate floating-point results may still vary between
worker counts due to summation order; worker-result merging will nevertheless
be deterministic by using a fixed worker order.

## Accumulation and reporting

Each worker will use private seismometer accumulation state while propagating
phonons. Caught-phonon energy and count updates will not contend on shared
seismometer bins. The coordinator will merge energy bins and counts after the
workers finish, then use the existing output path.

Aggregate lost, timeout, and invalid counts and diagnostic flags will be
merged similarly.

Detailed event reports will remain compatible with existing report settings,
but concurrent report lines must be serialized or buffered so individual lines
cannot interleave. Their order is not guaranteed when multiple workers are
active. Report mutexes are acquired only when the corresponding report is
enabled, so the default suppressed-report simulation path does not pay a
global mutex cost for disabled reports.

## Error handling

Worker exceptions will be caught and stored for the coordinator. A worker
failure will request cancellation, allow all workers to join, and then rethrow
a descriptive exception on the main thread. The coordinator will not emit the
successful simulation-complete marker or incomplete final traces after a
worker failure.

## Command-line interface

Recognize:

```text
--threads N
--seed N
```

`--threads N` requires a positive integer. `--seed N` requires a non-negative
integer. Existing command lines remain valid when neither option is supplied.
The simulation configuration summary will include the requested phonon count,
effective worker count, and effective random seed. Seed and thread settings
will also be included in model-parameter output where that output currently
records simulation parameters.

## Testing and verification

The implementation will be verified with the existing C++11/pthread build
system and no new external dependency:

1. Build the application with the existing Makefile.
2. Exercise automatic workers, `--threads 1`, explicit worker counts, invalid
   thread values, and explicit/automatic seeds.
3. Run small seeded simulations in single- and multi-worker modes and verify
   completion, phonon accounting, diagnostics, and valid trace output.
4. Repeat a seeded run with the same worker count and verify stable output.
5. Run a workload large enough to exercise all workers and inspect report
   output for corruption or interleaving.
6. Use thread sanitization where supported by the local compiler.
7. Compare single- and multi-worker wall-clock behavior as a performance smoke
   test, without imposing a fixed speedup requirement.

## Compatibility

The existing startup/completion markers, aggregate summaries, seismometer
trace format, report labels, and physical calculations will remain intact.
Only report ordering in multi-worker mode is intentionally unspecified.
