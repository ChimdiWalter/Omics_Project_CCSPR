# Bootstrap Stability — Arabidopsis All-Backend Attempt Log

**Date**: 2026-04-05
**Context**: After confirming 500-cell n_boot=50 is feasible at 39.45 s/iter (euclidean), user requested all 5 backends.

## Parallel Runs

### Run 1: euclidean + dtm (original)
- **Started**: 2026-04-05 ~20:31 UTC
- **Config**: `--max-cells 500 --n-boot 50 --skip-cv --backends euclidean dtm`
- **Output dir**: `results/arabidopsis_bootstrap/`
- **Timing probe**: euclidean = 39.45 s/iter → ~33 min for n_boot=50
- **Status**: IN PROGRESS

### Run 2: ricci + diffusion + phate_like (additional)
- **Started**: 2026-04-05 ~20:45 UTC (estimated)
- **Config**: `--max-cells 500 --n-boot 50 --skip-cv --backends ricci diffusion phate_like`
- **Output dir**: `results/arabidopsis_bootstrap_extra/`
- **Note**: ricci uses iters=2 (Ollivier-Ricci curvature iterations), may be slower per iteration
- **Status**: IN PROGRESS

## Fallback Policy

Per user instruction: **if Run 2 (ricci/diffusion/phate_like) fails, stick to euclidean + dtm results only**. This is consistent with the original Arabidopsis benchmark which only evaluated euclidean + dtm backends.

## Post-Run Actions

1. Merge per-backend checkpoint JSONs from both output directories
2. Create combined results CSV for manuscript reporting
3. Update feasibility log with final timings and results
4. All results, processes, and outputs to be used directly in manuscript updates

## Update: 2026-04-06 17:17 UTC

**Run 1 (euclidean+dtm)** has been running for ~21 hours with no completion logged. The euclidean backend is stuck on a single bootstrap iteration far longer than the 39 s/iter timing probe estimated. Per-iteration cost varies dramatically based on CVXPY solver convergence on each random subsample — pathological subsamples can take orders of magnitude longer.

**Run 2 (ricci+diffusion+phate_like)** showed similar issues:
- ricci: 160.11 s/iter probe (vs ~5 s on GSE161711) — would take ~133 min
- diffusion: 203.64 s/iter probe (vs ~0.9 s on GSE161711) — would take ~170 min
- phate_like: probe never completed
- **Killed at 2026-04-06 17:17 UTC** to free resources for Run 1

**Per user fallback policy**: Accepting euclidean + dtm only. If Run 1 also fails or runs indefinitely, fall back to GSE161711 results only and report Arabidopsis bootstrap stability as infeasible at the requested n_boot.

**Decision pending**: Wait for Run 1 to either complete or hit a wall. If no progress within next ~6 hours, kill and reduce n_boot further (e.g., n_boot=10 to match existing benchmark).

## Update: 2026-04-06 17:20 UTC — Run 1 killed

Run 1 (PID 1340352) had accumulated **1248 minutes (20.8 hours) of pure CPU time** on the euclidean backend without completing. The script does not log per-iteration progress within a backend, so there was no visibility into how many of the 50 bootstrap iterations had finished.

**Killed Run 1 manually**. No checkpoint files were saved (the script only writes per-backend checkpoints AFTER a backend completes its full n_boot iterations).

**Final state of n_boot=50, 500 cells attempt**:
- ZERO completed backends from either Run 1 or Run 2
- ZERO checkpoint files saved
- ZERO TIP vectors saved
- Only timing probe results captured in log files

## Root Cause Analysis: Why Arabidopsis is Pathologically Slow

The 50-feature PCA representation of Arabidopsis cells produces a **dense distance matrix** where every cell pair is non-trivially close. This causes:

1. **Larger Rips complex**: Many edges and triangles enter the filtration at low scales
2. **Larger boundary matrix**: The H1 boundary matrix passed to CVXPY is much larger than for sparse high-dimensional data like GSE161711
3. **Worse solver conditioning**: The CVXPY warning "Solution may be inaccurate" indicates the elastic-net problem is ill-conditioned, requiring many solver iterations
4. **Random per-iteration variance**: Some bootstrap subsamples produce particularly bad geometries that take orders of magnitude longer to solve

**Comparison**:
| Backend    | GSE161711 (96 samples, 5000 genes) | Arabidopsis (500 cells, 50 PCA) | Slowdown |
|------------|-----------------------------------|----------------------------------|----------|
| euclidean  | 19 s/iter                          | 39 s/iter (probe) → much higher in run | ~10–600× |
| ricci      | 13 s/iter                          | 160 s/iter (probe)               | ~12×     |
| diffusion  | 0.9 s/iter                         | 204 s/iter (probe)               | ~227×    |
| phate_like | 0.6 s/iter                         | not measured                      | —        |
| dtm        | 5.6 s/iter                         | 50 s/iter (probe)                 | ~9×      |

The slowdown is NOT explained by sample count alone (96 vs 500 = 5×). It is dominated by the **density of the distance matrix** in the PCA representation versus the sparser gene-space matrix in GSE161711.

## Update: 2026-04-06 17:21 UTC — Final Fallback Run

**Goal**: Get ANY Arabidopsis bootstrap stability result, even at n_boot=10 (matching existing benchmark).

**Run 3a (with timeout)**: `--max-cells 500 --n-boot 10 --backends euclidean dtm` wrapped in `timeout 3600`
- Started: 2026-04-06 17:21:18 UTC
- **Killed by user instruction** at 17:33 UTC: user wants no hard timeout if progress is being made
- The user reasoning: "dont time out at 1hr if results and progress is made"

**Run 3b (no timeout)**: Same config, no `timeout` wrapper
- Started: 2026-04-06 17:46 UTC
- **Completed**: 2026-04-06 23:53 UTC (6h 07m wall time)
- **Result**: SUCCESS — both backends 10/10 iterations

| Backend   | n_ok  | Support | Entropy | Gini  | Elapsed          | s/iter (probe → actual) |
|-----------|-------|---------|---------|-------|------------------|-------------------------|
| euclidean | 10/10 | 46      | 3.655   | 0.323 | 10,189 s (2.83h) | 31 → 1,019 (33×)        |
| dtm       | 10/10 | 48      | 3.733   | 0.289 | 11,609 s (3.22h) | 43 → 1,161 (27×)        |

The gap between probe timing and actual run timing confirms per-iteration variance from CVXPY solver conditioning on pathological subsample geometries. The probe used a different seed (SEED+9999) than the run (SEED+8000), so the subsamples differ.

**Stability observations**:
- Both backends highly concentrated: support sizes 46–48 out of 50 PCA features (~92–96% of features receive at least one vote across 10 bootstraps)
- Low Gini (0.29–0.32) and moderate entropy (3.66–3.73) indicate TIP mass is distributed across most features, which is expected given only 50 PCA components are available
- This is NOT a case where CC-SPR is selecting a tight subset — the PCA representation already compressed gene-level variation, so bootstrap-to-bootstrap feature rankings fluctuate across the 50 components

## Manuscript Reporting Plan (Pending Final Run)

Regardless of Run 3b outcome, the manuscript will report:

1. **GSE161711**: n_boot=1000, all 5 backends, complete success (8.3 hours total runtime)
2. **Arabidopsis**: One of the following, in order of preference:
   - **Best case**: Run 3b completes — n_boot=10, 500 cells, euclidean + dtm
   - **Fallback case**: Run 3b fails or runs indefinitely — report only existing benchmark (n_boot=10, 2999 cells from `run_arabidopsis.py`) with disclosure that higher n_boot was attempted but infeasible

The feasibility limitations themselves are a manuscript-worthy finding — they motivate future work on:
- Approximate persistent homology (sparse Rips, witness complexes)
- Better-conditioned cycle solvers (e.g., L1 regression instead of elastic net)
- Subsampling-aware variance reduction techniques
