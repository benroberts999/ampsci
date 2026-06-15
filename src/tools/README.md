# src/tools

Helper scripts. Run from the **project root**, e.g. `./src/tools/profile_hotspots.sh in.in`.

| script | purpose |
|---|---|
| `profile_hotspots.sh` | where CPU time goes: hot functions + call graph (perf) |
| `profile_alloc.sh` | heap allocation counts and sources (heaptrack) |
| `profile_cache.sh` | cache-miss rates + false-sharing check (perf stat / c2c) |
| `omp_scaling.sh` | OpenMP thread-scaling sweep: time vs thread count + efficiency |
| `profile_common.sh` | shared helpers, sourced by the `profile_*` scripts (not run directly) |
| `compile.sh` | build helper |
| `clang_format.sh`, `clang_format_all.sh` | apply clang-format |
| `clang-tidy.sh` | run clang-tidy |
| `lcov.sh` | coverage report |
| `bisect.sh` | git bisect helper |
| `build_includes.sh`, `license_year.sh` | maintenance |

---

## Profiling workflow

### Step 0: build a profiling binary (once)

```
make MODE=profile
```
- Same speed as a release build (`-O3 -DNDEBUG`), plus debug symbols.
- Lives in its own `build/profile/` tree; doesn't touch your normal `make` build.

### Step 1: is parallel scaling the problem? -- `omp_scaling.sh`

```
./src/tools/omp_scaling.sh --threads 24 --inc 4 --exe ./ampsci in.in
```
- Prints run time and efficiency vs thread count (efficiency = t1 / (tN * N), 1.0 = perfect).
- Efficiency drops off as threads increase -> the parallel part scales poorly.
- Stays high but time is still too long -> it's not parallelism; go to Step 2.

### Step 2: where is the time going? -- `profile_hotspots.sh`

```
./src/tools/profile_hotspots.sh in.in              # uses all cores
./src/tools/profile_hotspots.sh --threads 1 in.in  # 1 thread: cleanest hotspots
```
- One-time setup: `perf` needs `kernel.perf_event_paranoid <= 2`. The script tells
  you the exact `sudo sysctl ...` line if it's not set.
- Writes small reports into `./profile/hotspots/`; start with `SUMMARY.txt`:
  - `1_hot_functions.txt` - time spent *in* each function (self time)
  - `2_callgraph.txt` - top-down tree, time *under* each call path
  - `3_hot_callers.txt` - who calls the hottest function
- Reading it: `malloc`/`free`/`futex` near the top -> allocation/locking is the
  cost (Step 3a). Your own compute on top -> real work (Step 3b).

### Step 3a: is it allocations? -- `profile_alloc.sh`

```
./src/tools/profile_alloc.sh in.in        # needs: sudo apt install heaptrack
```
- `ALLOCATIONS.txt`: allocation counts, peak heap, and the worst call sites.
- A big count with small peak heap = many tiny short-lived allocations; only
  matters if `profile_hotspots.sh` shows malloc/free hot.
- Quick fix to test the hypothesis (no rebuild):
  `LD_PRELOAD=libtcmalloc_minimal.so.4 ./ampsci in.in`

### Step 3b: is it cache/memory bound? -- `profile_cache.sh`

```
./src/tools/profile_cache.sh in.in              # IPC + miss rates + false sharing
OMP_PROC_BIND=close OMP_PLACES=cores \
  ./src/tools/profile_cache.sh in.in            # pin to P-cores (hybrid CPUs)
```
- `CACHE.txt`: IPC and L1/LLC miss rates (low IPC + high LLC-miss = bandwidth-bound;
  more threads won't help -- that's usually the scaling ceiling).
- `FALSE_SHARING.txt`: cache lines bouncing between cores (`perf c2c`). Use
  `--no-c2c` to skip that (slower) second run.
- On Intel P/E hybrid CPUs, cache events show `<not supported>` unless pinned to
  P-cores (`OMP_PLACES=cores`).
