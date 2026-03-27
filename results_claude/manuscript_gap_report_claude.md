# Claude-Specific Manuscript Gap Report

Date: 2026-03-20
Repository: `/cluster/VAST/kazict-lab/e/lesion_phes/code/ccspr`
Output isolation: all new Claude outputs in `results_claude/` and `manuscript_claude/`

## 1. Already complete and reusable (no new computation needed)

### LUAD anchor benchmark
- `results/metrics.csv`: 2-split main benchmark (standard 0.767, harmonic 0.404, ccspr 0.313)
- `results/metrics_luad_multiseed.csv`: 18-pair multiseed statistics
- `results/paired_stats_luad_multiseed.csv`: paired t-test and Wilcoxon results
- Full ablation sweeps over: distance_mode, iters, k, lambda, top_k, normalize
- Diagnostic figures: TIP, H1 lifetime, H1 prominence, lambda robustness, paired deltas
- Best ablation summary table

### GSE161711 bulk disease-aligned validation
- `results/metrics_cll_venetoclax_gse161711.csv`: 5-split main benchmark (standard 0.971, harmonic 0.585, ccspr 0.603)
- Full ablation sweeps over all hyperparameter axes
- Diagnostic figures: TIP, H1 lifetime, H1 prominence, lambda robustness
- CC-SPR modestly exceeds harmonic (+0.018) while trailing standard

### GSE165087 bounded single-cell validation
- `results/metrics_cll_rs_scrna_gse165087.csv`: 1-split bounded (all models 0.925)
- Completed manuscript-safe run (job 12874303, 7min, 2.9 GB RSS)
- Diagnostic figures: TIP, H1 lifetime, H1 prominence, F1 bars

### BRCA supporting evidence
- `results/brca/metrics.csv`: 2-split (standard 0.611, harmonic 0.479, ccspr 0.410)

### Geometry study artifacts
- `results/tables/geometry_backend_summary.csv`
- `results/tables/geometry_diagnostics_summary.csv`
- `results/tables/representative_diagnostics_summary.csv`
- `results/tables/geometry_runtime_cost_summary.csv`
- `results/figures/geometry_backend_summary.png`
- All ablation profile tables and figures

### Existing manuscript
- `manuscript/cc_spr_full_paper.tex` (~357 lines, ~10 pages)
- `manuscript/cc_spr_full_supplement.tex` (~189 lines)
- Both compiled to PDF

## 2. What the paper still needs (addressable without new computation)

1. **Mathematical depth**: The existing paper has correct but brief math. A stronger submission needs:
   - Formal definition environment for the filtration stability proposition
   - Deeper treatment of the elastic net connection to LASSO/ridge
   - Explicit algorithm pseudocode
   - Convergence/optimality remark for the convex sparse objective
   - More rigorous backend perturbation analysis

2. **Harmonic PH comparison**: Needs expansion beyond the single table. Should discuss:
   - Computational complexity comparison
   - When harmonic representatives are actually sufficient
   - What sparse representatives uniquely enable

3. **Geometry backend sections**: Currently 1-2 paragraphs each. Need:
   - Diffusion: heat kernel definition, connection to graph Laplacian, multiscale property
   - PHATE-like: potential distance definition, information-geometric interpretation
   - DTM: formal definition, robustness guarantee sketch, connection to Wasserstein

4. **Results narrative**: Needs more careful interpretation:
   - Why standard baseline dominance does not invalidate the method
   - Backend-sensitive interpretation across datasets
   - Clearer ablation discussion connecting hyperparameter sensitivity to method design

5. **Bibliography**: Currently only 3 entries. Needs real references.

6. **Author/affiliation**: Currently placeholder.

## 3. What is missing and too expensive

- Diffusion, PHATE-like, DTM backend implementations and experiments
- Large-scale scRNA ablation sweeps (any config > manuscript-safe)
- Survival/enrichment analyses
- Cross-tool harmonic PH comparison (maTilDA parity)
- Per-bootstrap support entropy/Gini recomputation

## 4. What can be strengthened without new computation

- All mathematical sections (pure writing)
- All comparison tables (pure writing)
- Figure captions and interpretation (pure writing)
- Bibliography expansion (pure writing)
- Limitations and honest claims section (pure writing)
- Supplement execution envelope documentation (pure writing)
- Algorithm pseudocode (pure writing)

## Decision

The Claude-specific manuscript upgrade will focus entirely on writing quality and mathematical depth using existing results. No new computation is needed. All outputs go to `manuscript_claude/` and `results_claude/`.
