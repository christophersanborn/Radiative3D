# Multithreaded Simulation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Safely parallelize independent phonon propagation with automatic worker selection, reproducible seeds, and merged per-worker envelope results.

**Architecture:** Add explicit simulation controls to `ModelParams` and command-line parsing. Replace shared `rand()`/fixed-four-thread execution with a worker pool whose phonon work is claimed atomically and whose random, accumulation, and diagnostics state is isolated per worker; merge results in a fixed order before existing output. Keep model geometry shared read-only and preserve existing trace formats.

**Tech Stack:** C++11, `std::thread`, `std::atomic`, `std::mt19937_64`, existing Makefile/pthread build, shell-based verification.

**Spec:** `docs/superpowers/specs/2026-10-15-multithread-simulation-design.md`

## Global Constraints

- Use the existing C++11/pthread build system and no new external dependency.
- Multithreading is automatic by default; `--threads 1` must force single-thread execution.
- `--threads N` requires a positive integer; explicit values may oversubscribe hardware.
- `--seed N` requires a non-negative integer; omitted seeds are generated automatically and printed.
- Do not change physical propagation algorithms, Earth-model geometry, envelope normalization, or public trace format.
- Preserve existing startup/completion markers, aggregate summaries, and report labels.
- Worker-local results must be merged in a fixed order; detailed report ordering is unspecified with multiple workers.

## Review Focus

- Zero, negative, malformed, or overflowing thread values must fail with a clear command-line error — covered in Task 1 CLI validation tests.
- Zero, negative, malformed, or overflowing seed values must fail clearly while valid unsigned seeds survive parsing — covered in Task 1 seed tests.
- Hardware concurrency reporting may return zero and very small phonon counts must not create zero workers — covered in Task 2 worker-count tests.
- A worker exception must cancel/join the pool and suppress successful completion output — covered in Task 4 failure-path test.
- Concurrent reporting must not corrupt individual report lines or final traces — covered in Task 5 multi-worker output smoke test.

---

### Task 1: Add simulation controls and deterministic random context

**Files:**
- Modify: `cmdline.hpp` — add option IDs and registrations for `--threads` and `--seed`.
- Modify: `main.cpp` — parse and validate the two options into `ModelParams`.
- Modify: `model.hpp` — add thread/seed parameters and accessors used by `Model`.
- Modify: `model.cpp` — include effective seed initialization and parameter reporting.
- Modify: `probability.hpp`, `probability.cpp`, `sources.hpp`, `sources.cpp`, `scatterers.hpp`, `scatterers.cpp`, `phonons.hpp`, `phonons.cpp` — route simulation random draws through an explicit per-phonon/per-worker random context instead of global `rand()` calls.

**Interfaces:**
- Produces `ModelParams::NumThreads` as a positive override value, with zero meaning automatic selection.
- Produces `ModelParams::RandomSeed` as an unsigned 64-bit value plus a boolean indicating whether the user supplied it.
- Produces an explicit random-context interface used by source generation, probability selection, scattering, and propagation.

- [ ] **Step 1: Add failing CLI tests/checks** for `--threads 1`, valid `--seed`, malformed values, zero/negative values, and values that exceed the signed integer range. Since the repository has no test framework, use a small shell test script under `tests/` that builds/runs the executable and asserts exit status and diagnostic text.
- [ ] **Step 2: Run the new checks** and confirm they fail because the options are not recognized or validated.
- [ ] **Step 3: Implement option IDs, token registration, and `process_option()` handling** with positive-thread and non-negative-seed validation; store values in `ModelParams` without changing existing defaults.
- [ ] **Step 4: Add the random context interface** with a per-context `std::mt19937_64` engine and helpers for uniform real/index selection, then update all simulation-path random consumers to accept or retrieve that context explicitly. Preserve probability behavior and avoid global mutable RNG state.
- [ ] **Step 5: Make phonon generation identify its random context** so each claimed phonon can derive a stream from the base seed and phonon ID. Keep non-simulation diagnostic APIs source-compatible where practical.
- [ ] **Step 6: Run the CLI checks** and build the application; expect all validation checks to pass and no global-RNG compile errors.
- [ ] **Step 7: Commit** with `git add` of the CLI, parameter, RNG, and test files and message `Add reproducible simulation controls`.

### Task 2: Implement worker-count policy and work claiming

**Files:**
- Modify: `model.hpp` — define worker-pool state and worker method signatures.
- Modify: `model.cpp` — implement automatic worker selection, explicit override handling, atomic phonon-ID claiming, progress coordination, and pool lifecycle.
- Test: `tests/test_worker_policy.sh` — exercise automatic selection indirectly through startup output and small workloads.

**Interfaces:**
- Consumes `ModelParams::NumThreads` and `ModelParams::RandomSeed` from Task 1.
- Produces `Model::ResolveWorkerCount()` behavior: explicit positive override wins; otherwise reserve one hardware thread when hardware concurrency is greater than two; clamp to requested phonon count and minimum one worker for positive workloads.
- Produces worker callback `SimulationThread(unsigned worker_index)` that claims stable phonon IDs from an atomic counter.

- [ ] **Step 1: Add failing worker-policy checks** for `--threads 1`, an explicit multi-worker value, automatic mode, and a workload with fewer phonons than hardware threads; assert the printed effective worker count.
- [ ] **Step 2: Run the checks** and confirm the current fixed-four-thread implementation fails at least the requested-count and explicit-count assertions.
- [ ] **Step 3: Implement `ResolveWorkerCount()`** using `hardware_concurrency()` with a zero result fallback, automatic reservation policy, explicit override, and phonon-count clamping.
- [ ] **Step 4: Replace the four hard-coded `std::thread` objects** with a vector of workers created from the resolved count. Use an atomic next-phonon counter or remaining counter that yields each ID exactly once.
- [ ] **Step 5: Move progress reporting to the coordinator/shared counters** so worker threads do not write progress output directly or interleave progress lines.
- [ ] **Step 6: Run worker-policy checks and a build**; verify exact phonon accounting for one and multiple workers.
- [ ] **Step 7: Commit** with message `Add automatic simulation worker pool`.

### Task 3: Isolate worker random and event-source state

**Files:**
- Modify: `events.hpp`, `events.cpp`, `sources.hpp`, `sources.cpp` — ensure event generation can use worker-local random context without mutating shared generator state.
- Modify: `model.cpp` — derive each phonon stream from the effective seed and stable phonon ID and pass it into event generation/propagation.
- Modify: `phonons.hpp`, `phonons.cpp`, `media.hpp`, `media.cpp`, `scatterers.hpp`, `scatterers.cpp` — update propagation/scattering random calls to use the active explicit context.
- Test: `tests/test_seed_reproducibility.sh` — compare repeated seeded outputs.

**Interfaces:**
- Consumes worker IDs from Task 2 and the random-context API from Task 1.
- Produces deterministic per-phonon random initialization from `(base_seed, phonon_id)` independent of worker scheduling.

- [ ] **Step 1: Add failing reproducibility checks** that run the same small model twice with the same seed and `--threads 1`, then compare generated traces and summaries.
- [ ] **Step 2: Run the checks** and confirm current time-based/global RNG seeding makes the result non-reproducible or cannot accept `--seed`.
- [ ] **Step 3: Implement stable stream derivation** using a documented 64-bit mixing function or `std::seed_seq` inputs from the base seed and phonon ID; construct the context before each phonon and ensure all random consumers use it.
- [ ] **Step 4: Remove simulation dependence on `srand()`/`rand()`** while leaving unrelated legacy APIs either converted or explicitly isolated from the simulation path.
- [ ] **Step 5: Run single-thread reproducibility checks** and verify identical trace files for repeated seeded runs.
- [ ] **Step 6: Run a multi-worker seeded smoke test** and verify all requested phonons complete; document if aggregate floating-point output differs by worker count while per-phonon stream assignment remains stable.
- [ ] **Step 7: Commit** with message `Isolate phonon random streams`.

### Task 4: Add worker-local accumulation and failure handling

**Files:**
- Modify: `dataout.hpp`, `dataout.cpp` — expose a worker-local accumulation snapshot/clone/merge interface for seismometer bins and aggregate counters without changing final output format.
- Modify: `model.hpp`, `model.cpp` — allocate one accumulation context per worker, merge contexts in worker-index order, capture exceptions, request cancellation, join all workers, and rethrow on the coordinator thread.
- Modify: `phonons.cpp` — report collection/lost/timeout/invalid results to the active worker context rather than shared mutable output state.
- Test: `tests/test_worker_failure.sh` — exercise an invalid/pathological configuration and assert no false completion marker.

**Interfaces:**
- Produces `DataReporter::WorkerAccumulator` (or equivalent focused type) containing per-seismometer bins, counts, lost/timeout/invalid totals, and invalid diagnostic flags.
- Produces `DataReporter::MergeWorkerAccumulator(const WorkerAccumulator&)` with deterministic merge order controlled by `Model`.
- Produces cancellation state checked between phonons and at safe propagation/reporting boundaries.

- [ ] **Step 1: Add failing accumulation checks** that run a small multi-worker simulation and verify trace output is valid and aggregate counts are present; add a failure-path assertion for absent completion marker.
- [ ] **Step 2: Implement worker accumulator storage** that mirrors only mutable seismometer and diagnostic data needed by the simulation, leaving metadata and final file output centralized.
- [ ] **Step 3: Make collection and aggregate report methods target the active worker accumulator** and keep detailed event reporting out of the worker-local aggregate path.
- [ ] **Step 4: Merge accumulators in fixed worker-index order** before invoking existing post-simulation trace/summary output.
- [ ] **Step 5: Add worker exception capture and cancellation** using `std::exception_ptr`, atomic cancellation, guaranteed joins, and coordinator-thread rethrow.
- [ ] **Step 6: Run accumulation, failure, and build checks**; verify no successful completion marker or final traces are emitted after a worker error.
- [ ] **Step 7: Commit** with message `Merge thread-local simulation results safely`.

### Task 5: Serialize detailed reporting and update documentation

**Files:**
- Modify: `dataout.hpp`, `dataout.cpp` — protect detailed report line emission with a mutex or equivalent atomic-line buffer; keep report labels and formats unchanged.
- Modify: `main.cpp`, `model.cpp`, `model.hpp` — print effective phonon count, worker count, and seed; include thread/seed values in parameter output.
- Modify: `README.md`, `doc/MANUAL.md` — document automatic workers, `--threads`, `--seed`, reproducibility expectations, and report-order behavior.
- Test: `tests/test_multithread_output.sh` — run a multi-worker report-enabled workload and check for complete, non-interleaved lines and valid final traces.

**Interfaces:**
- Consumes merged results and effective settings from Tasks 2–4.
- Preserves existing `DataReporter` report methods and output identifiers.

- [ ] **Step 1: Add failing output checks** for startup configuration fields and report-enabled multi-worker output.
- [ ] **Step 2: Add synchronized detailed-report emission** so each report record is written as one protected unit; do not serialize aggregate bin updates.
- [ ] **Step 3: Print effective simulation configuration** at simulation start and include seed/thread settings in existing parameter serialization without changing trace schema.
- [ ] **Step 4: Document CLI usage and caveats** including automatic seed display, `--threads 1` debugging, non-guaranteed event ordering, and possible floating-point merge differences across worker counts.
- [ ] **Step 5: Run the output test and build**; inspect generated report and trace files for malformed/interleaved records.
- [ ] **Step 6: Commit** with message `Document and harden multithreaded output`.

### Task 6: Full verification and performance smoke test

**Files:**
- Modify: `tests/` only if earlier checks need consolidation or portability fixes.

- [ ] **Step 1: Build from a clean object state** with `make cleanall` followed by `make`.
- [ ] **Step 2: Run all focused shell checks** for CLI validation, worker policy, reproducibility, accumulation/failure, and output integrity.
- [ ] **Step 3: Run a larger workload** with `--threads 1` and automatic/multiple workers; verify exact requested phonon accounting, completion, valid traces, and no data races observable in output.
- [ ] **Step 4: Run compiler thread sanitization if supported** and record any platform limitations; fix confirmed races before proceeding.
- [ ] **Step 5: Measure a representative workload** with one worker and automatic workers, reporting observed wall-clock behavior without asserting a fixed speedup.
- [ ] **Step 6: Review `git diff`, `git diff --check`, and branch status**, then commit any final test-only fixes with an explanatory message.
