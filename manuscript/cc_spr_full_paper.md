# Geometry-Sensitive Persistent Representatives for Omics

Full manuscript source: `manuscript/cc_spr_full_paper.tex`.

This version positions harmonic persistent homology as the parent framework and CC-SPR as a sparse representative extension with a geometry backend family:
- Euclidean
- Ollivier-Ricci flow
- diffusion geometry
- PHATE-like potential geometry
- DTM-robust geometry

## Added in this upgrade

- Formal method expansion with notation, objectives, and proposition-style statements.
- Explicit sections:
  - `Diffusion geometry, heat kernels, and manifold learning`
  - `Information geometry and PHATE-like distances`
  - `Robust topology via DTM`
- LUAD multi-seed paired benchmark:
  - `results/tables/luad_multiseed_main_benchmark.csv`
  - `results/tables/luad_multiseed_paired_stats.csv`
  - `results/figures/luad_multiseed_benchmark.png`
- LUAD biological diagnostics (pilot):
  - `results/tables/luad_pathway_coherence.csv`
  - `results/tables/luad_survival_stratification.csv`
  - `results/figures/luad_biological_diagnostics.png`
  - `results/figures/luad_survival_km_by_method.png`
- Runtime/cost profiling:
  - `results/tables/runtime_cost_geometry_overall.csv`
  - `results/figures/runtime_cost_geometry.png`
  - `results/figures/runtime_cost_geometry_heatmap.png`

## Compile

```bash
pdflatex -interaction=nonstopmode -halt-on-error -output-directory manuscript manuscript/cc_spr_full_paper.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory manuscript manuscript/cc_spr_full_paper.tex
```
