# CC-SPR Geometry-Family Completion Plan

**Date**: 2026-03-21
**Objective**: Upgrade the merged manuscript into a genuinely completed geometry-family empirical study by implementing and benchmarking diffusion, PHATE-like, and DTM backends alongside the existing Euclidean and Ricci backends.

---

## 1. Current Empirical Coverage by Backend × Dataset

| Backend | LUAD | GSE161711 | GSE165087 (bounded) | BRCA |
|---------|------|-----------|---------------------|------|
| **Euclidean** | ✓ Full (2-split + 18-seed + ablation) | ✓ Full (5-fold + ablation) | ✓ Bounded (1-split) | ✓ Partial (2-fold) |
| **Ricci** | ✓ Full (2-split + 18-seed + ablation) | ✓ Full (5-fold + ablation) | ✓ Bounded (1-split) | — |
| **Diffusion** | ✗ Not implemented | ✗ Not implemented | ✗ Not implemented | — |
| **PHATE-like** | ✗ Not implemented | ✗ Not implemented | ✗ Not implemented | — |
| **DTM** | ✗ Not implemented | ✗ Not implemented | ✗ Not implemented | — |

## 2. Missing Backend × Dataset Cells

| Cell | Priority | Feasibility | Plan |
|------|----------|-------------|------|
| Diffusion × LUAD | **P1** | Easy — n=275, lightweight | Implement + run 2-split benchmark |
| Diffusion × GSE161711 | **P2** | Easy — n=96, lightweight | Implement + run 5-fold CV |
| Diffusion × GSE165087 | **P3** | Feasible — bounded 320 cells | Implement + run 1-split bounded |
| PHATE-like × LUAD | **P1** | Easy — derive from diffusion operator | Implement + run 2-split benchmark |
| PHATE-like × GSE161711 | **P2** | Easy — derive from diffusion operator | Implement + run 5-fold CV |
| PHATE-like × GSE165087 | **P3** | Feasible — bounded | Implement + run 1-split bounded |
| DTM × LUAD | **P1** | Easy — direct computation | Implement + run 2-split benchmark |
| DTM × GSE161711 | **P2** | Easy — direct computation | Implement + run 5-fold CV |
| DTM × GSE165087 | **P3** | Feasible — bounded | Implement + run 1-split bounded |

**All 9 cells can be filled.** None require heavy compute.

## 3. Backend Implementation Strategy

### 3.1 Diffusion Geometry
- Build shared kNN affinity graph (reuse sklearn `kneighbors_graph`)
- Compute Gaussian kernel weights: `W_ij = exp(-d_ij^2 / (2 * sigma_i * sigma_j))` where `sigma_i` = distance to k-th nearest neighbor
- Row-normalize to transition matrix P
- Compute diffusion distance at time t: `D_diff(i,j)^2 = sum_l (lambda_l^t)^2 * (psi_l(i) - psi_l(j))^2`
- Use truncated eigendecomposition via `scipy.sparse.linalg.eigsh` (top 50 eigenvectors)
- Lightweight: only matrix-vector products on sparse matrix, O(n * k * n_eig)

### 3.2 PHATE-like (Potential Distance)
- Derive from same diffusion operator P
- Compute potential: `U_t(i) = -log(P^t(i, :) + epsilon)` (information-geometric embedding)
- Potential distance: `D_pot(i,j) = ||U_t(i) - U_t(j)||_2`
- This is a lightweight surrogate for full PHATE (avoids MDS embedding step)
- Honest in manuscript: described as "PHATE-like potential distance" not full PHATE

### 3.3 DTM (Distance-to-Measure)
- For each point i, compute DTM weight: `d_DTM(i) = (1/k * sum_{j in kNN(i)} d(i,j)^2)^{1/2}`
- Weighted distance: `D_DTM(i,j) = sqrt(d(i,j)^2 + d_DTM(i)^2 + d_DTM(j)^2)`
- Straightforward, no eigendecomposition needed
- Robust to outliers by construction

### 3.4 Shared Infrastructure
- All new backends share the same kNN graph construction
- All produce n×n distance matrices compatible with existing pipeline
- All use the same caching infrastructure (SHA1-keyed)
- No changes needed to topology, solver, TIP, or evaluation code

## 4. Experimental Design

### 4.1 LUAD (Anchor Benchmark)
- **Splits**: 2 (same splits as existing, seed=42)
- **Variable genes**: 5000 (same as existing)
- **Bootstrap**: n_boot=5, subsample_frac=0.8
- **Backends**: euclidean, ricci, diffusion, phate_like, dtm
- **Solver**: elastic (CC-SPR default)
- **Parameters**: k=10, lambda=1e-6, top_k=20, normalize=std
- **Ricci iters**: 5 (for ricci only)
- **Diffusion t**: 2 (short time)
- **DTM k_dtm**: 10 (same as kNN k)
- **Collect**: F1, H1 lifetime, H1 prominence, TIP scores, support size

### 4.2 GSE161711 (Bulk Validation)
- **Splits**: 5 (5-fold CV, same as existing)
- **Variable genes**: 5000
- **Bootstrap**: n_boot=5, subsample_frac=0.8
- **Backends**: all 5
- **Same parameters as LUAD**
- **Collect**: same metrics

### 4.3 GSE165087 (Bounded Single-Cell)
- **Max cells**: 320
- **n_hvg**: 800, n_pcs: 20
- **Splits**: 1
- **Bootstrap**: n_boot=1
- **Backends**: all 5
- **Same parameters** (adjusted for bounded data)
- **Collect**: F1 only (no multi-split statistics)

### 4.4 Lightweight Ablations (LUAD + GSE161711 only)
- Backend comparison (primary): all 5 backends with fixed params
- Lambda sensitivity: [1e-6, 1e-5] on all 5 backends
- k sensitivity: [10, 15] on all 5 backends
- top_k sensitivity: [10, 20] on all 5 backends

**Total configurations**: 5 backends × 3 ablation axes × 2 values × 2 datasets = 60 lightweight configs
Plus 5 main configs × 2 datasets + 5 × 1 bounded = 15 main configs
**Grand total**: ~75 lightweight experiment configs

## 5. Resource Estimates

| Dataset | Per-backend (est.) | Total 5 backends |
|---------|-------------------|-----------------|
| LUAD (2-split) | ~2-5 min | ~15-25 min |
| GSE161711 (5-fold) | ~3-8 min | ~20-40 min |
| GSE165087 (1-split) | ~1-3 min | ~5-15 min |
| Ablations (LUAD) | ~1-2 min each | ~30-60 min |
| Ablations (GSE161711) | ~1-2 min each | ~30-60 min |
| **Total** | | **~2-3 hours** |

Memory: Peak <8 GB (LUAD n=275, GSE161711 n=96, GSE165087 n=320 bounded)
This is well within safe resource limits.

## 6. Analyses Intentionally Left for Future Work

1. Full-scale scRNA evaluation beyond bounded configurations
2. Multi-split scRNA statistics
3. Per-bootstrap entropy/Gini caching
4. Pathway enrichment analysis (KEGG, Reactome, MSigDB)
5. Survival stratification
6. Full PHATE pipeline (with MDS) — we use lightweight potential distance surrogate
7. maTilDA/R harmonic PH parity testing
8. Higher homology dimensions (H2+)
9. Large bootstrap counts for TIP convergence
10. Cross-backend × parameter Cartesian product sweeps on scRNA

## 7. Final Manuscript Structure

### Main Paper
1. Introduction (5 subsections, same as merged)
2. Background & Relationship to Harmonic PH (expanded comparison table)
3. Mathematical Formulation (chain complexes, objectives, stability)
4. Geometry Backend Family (5 subsections, now all empirically validated)
5. Topological Influence Profile
6. Experimental Design (now with 5-backend geometry-family protocol)
7. Results (organized by scientific question, now with 5-backend data)
8. Resource-Bounded Single-Cell Validation (5-backend bounded results)
9. Biological Interpretation (strengthened with geometry-family insights)
10. Limitations & Future Directions
11. Related Work
12. Conclusion

### Key New Tables
1. **Geometry-Family Benchmark Table**: 5 backends × 3 datasets (+ standard/harmonic baselines)
2. **Paired Backend Comparison on LUAD**: mean differences + Wilcoxon p-values
3. **Representative Diagnostics by Backend**: support size, H1 lifetime, H1 prominence per backend
4. **Backend Ablation Sensitivity**: lambda, k, top_k sensitivity per backend
5. **Harmonic PH vs CC-SPR Comparison Table**: updated with all backends
6. **Runtime/Memory Summary**: per backend per dataset

### Key New Figures
1. F1 by backend (grouped bar chart) — LUAD and GSE161711
2. H1 lifetime/prominence by backend (boxplots)
3. TIP concentration by backend
4. Backend sensitivity (lambda/k/top_k) line plots
5. Architecture diagram (if quality sufficient)

### Supplement
- Full ablation tables
- Per-dataset per-backend detailed metrics
- GSE165087 bounded diagnostics
- Cross-dataset divergence analysis
- Complete runtime logs

## 8. Output Files

| File | Status |
|------|--------|
| `src/ccspr/geometry/distance.py` | To update (add 3 backends) |
| `experiments/run_geometry_family.py` | To create |
| `results/geometry_family/` | To create (all new results) |
| `manuscript_claude/cc_spr_final_geometry_complete.tex` | To create |
| `manuscript_claude/cc_spr_final_geometry_complete_supplement.tex` | To create |
| `manuscript_claude/cc_spr_final_geometry_complete.pdf` | To compile |
| `manuscript_claude/cc_spr_final_geometry_complete_supplement.pdf` | To compile |
| `results_claude/geometry_completion_plan.md` | This file |
| `results_claude/geometry_completion_report.md` | To create |
