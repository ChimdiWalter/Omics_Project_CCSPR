# Resource-Safe Execution Note (Claude)

Date: 2026-03-20
Repository: `/cluster/VAST/kazict-lab/e/lesion_phes/code/ccspr`

## Heavy runs that failed (not repeated)

| Job ID | Config | Outcome | Elapsed | Mem request | Peak RSS |
|--------|--------|---------|---------|-------------|----------|
| 12833143 | cll_rs_scrna.yaml | OUT_OF_MEMORY | 16:43:25 | 96G | ~100.3 GB |
| 12844071 | cll_rs_scrna_fasttrack.yaml | TIMEOUT | 20:00:06 | 96G | ~75.2 GB |
| 12849765 | cll_rs_scrna_lowmem.yaml | OUT_OF_MEMORY | 11:25:59 | 48G | ~49.7 GB |

These configurations were intentionally not repeated because they established a consistent failure envelope.

## Successful bounded run (reused)

| Job ID | Config | Outcome | Elapsed | Mem request | Peak RSS |
|--------|--------|---------|---------|-------------|----------|
| 12874303 | cll_rs_scrna_manuscript_safe.yaml | COMPLETED | 00:07:49 | 32G | ~2.9 GB |

## Geometry results reused (no recomputation)

All of the following were reused from existing repository outputs:
- LUAD 2-split main benchmark (standard 0.767, harmonic 0.404, ccspr 0.313)
- LUAD 18-pair multiseed paired statistics (Ricci vs Eu sparse: +0.046, p=0.308)
- LUAD full ablation sweeps (6 axes: distance_mode, iters, k, lambda, normalize, top_k)
- GSE161711 5-fold main benchmark (standard 0.971, harmonic 0.585, ccspr 0.603)
- GSE161711 full ablation sweeps (6 axes)
- GSE165087 bounded 1-split (all models 0.925)
- BRCA 2-split (standard 0.611, harmonic 0.479, ccspr 0.410)
- All existing figures (18 PNG files under results/figures/)
- All existing tables (9 CSV files under results/tables/)
- Geometry backend summary, diagnostics, and runtime tables

## High-cost analyses intentionally not repeated

1. Any scRNA config larger than manuscript-safe
2. Diffusion / PHATE-like / DTM backend implementation and experiments
3. Survival/enrichment analyses
4. Per-bootstrap support entropy/Gini recomputation
5. Cross-tool harmonic PH comparison (maTilDA parity)
6. Large single-cell ablation sweeps

## Why the manuscript is scientifically defensible

1. LUAD is a complete anchor benchmark with multiseed paired statistics
2. GSE161711 provides disease-aligned bulk validation with full ablation
3. BRCA provides supplementary consistency
4. Bounded GSE165087 establishes single-cell feasibility
5. All empirical claims are restricted to implemented backends (Euclidean, Ricci)
6. Diffusion/PHATE/DTM are presented as mathematical extensions, not empirical claims
7. Resource failures are documented transparently
8. The Ricci benefit is reported honestly as non-monotone and dataset-dependent

## Claude-specific output isolation

All new Claude outputs are written to:
- `manuscript_claude/` (LaTeX sources and PDFs)
- `results_claude/` (gap report, execution note)

No existing files in `manuscript/`, `results/`, or any other directory were modified.
