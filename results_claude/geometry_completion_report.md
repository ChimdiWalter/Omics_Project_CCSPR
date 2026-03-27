# CC-SPR Geometry-Family Completion Report

## Summary

All five geometry backends (Euclidean, Ollivier-Ricci, Diffusion, PHATE-like, DTM) are now fully implemented, benchmarked, and documented. The manuscript has been upgraded from a two-backend paper to a complete geometry-family empirical study using exclusively real experimental data.

## What was done

### 1. Backend implementation (`src/ccspr/geometry/distance.py`)
- Added `_build_knn_affinity()`: shared kNN affinity with adaptive Gaussian kernel
- Added `_diffusion_distance()`: truncated eigendecomposition of Markov transition matrix (50 eigenvectors)
- Added `_phate_like_distance()`: potential distance via `-log(P^t + eps)`, lightweight PHATE surrogate
- Added `_dtm_distance()`: distance-to-measure weighted distance: `sqrt(d^2 + dtm_i^2 + dtm_j^2)`
- Updated `build_distance()` dispatcher to handle all 5 modes

### 2. Experiment runner (`experiments/run_geometry_family.py`)
- Phase 1: LUAD (max_samples=60, 2-split, n_boot=5)
- Phase 2: GSE161711 (n=96, 5-fold CV, n_boot=5)
- Phase 3: GSE165087 (320 bounded cells, 1-split, n_boot=1)
- Collects: F1, H1 diagnostics, TIP entropy/Gini/support, lambda ablation
- Checkpoint after each phase for fault tolerance

### 3. Slurm execution (job 12883024)
- Script: `slurm/hellbender_geometry_family.sbatch`
- Resources: 48 GB RAM, 8 CPUs, 12h walltime
- Status: COMPLETED
- All 15 backend x dataset combinations succeeded
- Slowest: DTM on GSE165087 (~65 min due to dense filtration on 320 samples)

### 4. Results (all real data, no placeholders)

**Benchmark F1 (mean +/- std):**

| Dataset | Euclidean | Ricci | Diffusion | PHATE-like | DTM |
|---------|-----------|-------|-----------|------------|-----|
| LUAD (n=60, 2-split) | 0.434+/-0.039 | 0.501+/-0.107 | 0.500+/-0.100 | 0.355+/-0.088 | 0.532+/-0.077 |
| GSE161711 (n=96, 5-fold) | 0.748+/-0.083 | 0.594+/-0.085 | 0.555+/-0.072 | 0.540+/-0.063 | 0.618+/-0.133 |
| GSE165087 (n=320, 1-split) | 0.925 | 0.925 | 0.925 | 0.925 | 0.925 |

**Key findings:**
- Backend ranking is dataset-dependent (central result)
- LUAD: DTM > Ricci > Diffusion > Euclidean > PHATE-like
- GSE161711: Euclidean > DTM > Ricci > Diffusion > PHATE-like
- GSE165087: ceiling effect (all identical at 0.925)
- Lambda ablation: LUAD insensitive; GSE161711 mildly sensitive
- Curvature backends produce higher TIP entropy, broader support, lower Gini

### 5. Manuscript outputs

| File | Pages | Size |
|------|-------|------|
| `manuscript_claude/cc_spr_final_geometry_complete.tex` | 29 pp | 1390 lines |
| `manuscript_claude/cc_spr_final_geometry_complete.pdf` | 29 pp | 859 KB |
| `manuscript_claude/cc_spr_final_geometry_complete_supplement.tex` | 9 pp | ~300 lines |
| `manuscript_claude/cc_spr_final_geometry_complete_supplement.pdf` | 9 pp | 198 KB |

Previous files (`cc_spr_final_merged.*`) are untouched.

### 6. Data artifacts

| File | Description |
|------|-------------|
| `results/geometry_family/geometry_family_results.json` | Full results JSON |
| `results/geometry_family/checkpoint.json` | Phase completion checkpoint |
| `results/geometry_family/tables/geometry_family_benchmark.csv` | F1 by backend x dataset |
| `results/geometry_family/tables/geometry_family_diagnostics.csv` | H1 diagnostics |
| `results/geometry_family/tables/geometry_family_ablation.csv` | Lambda ablation |
| `results/geometry_family/tables/geometry_family_timing.csv` | Wall-clock timing |

## What changed from previous manuscript

1. **Coverage**: 2 backends (Euclidean, Ricci) -> 5 backends (+ Diffusion, PHATE-like, DTM)
2. **Framing**: "mathematical extensions" -> "fully benchmarked backends"
3. **Tables**: All benchmark, diagnostics, ablation, and timing tables contain real numbers
4. **Honest reporting**: GSE165087 ceiling effect acknowledged; dataset-dependent ranking emphasized
5. **New supplement**: Dedicated to geometry-family results with full diagnostics and timing

## Integrity notes

- Every number in the manuscript traces to a CSV in `results/geometry_family/tables/`
- No placeholder values anywhere in the manuscript
- GSE165087 ceiling (all F1 = 0.925) reported honestly as a bounded feasibility check
- The paper maintains the framing that harmonic PH is the parent method and CC-SPR is the sparse extension
