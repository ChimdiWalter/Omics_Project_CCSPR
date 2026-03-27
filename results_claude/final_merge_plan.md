# CC-SPR Final Manuscript Merge Plan

**Date**: 2026-03-21
**Objective**: Produce one unified final manuscript by merging the strongest parts of all existing manuscript families.

---

## 1. Manuscript Inventory

| Manuscript | Path | Lines | Role |
|---|---|---|---|
| Claude full paper | `manuscript_claude/cc_spr_full_paper_claude.tex` | 871 | **Base shell** — strongest math, geometry-family framing, biology |
| Claude supplement | `manuscript_claude/cc_spr_full_supplement_claude.tex` | 381 | Strongest supplement — execution envelope, code-level correction, quantitative summaries |
| Hellbender paper | `manuscript/cc_spr_full_paper.tex` | 357 | Best resource-bounded language, honest execution-envelope reporting |
| Hellbender supplement | `manuscript/cc_spr_full_supplement.tex` | 188 | Concise operational supplement |

## 2. Base Shell Decision

**The Claude manuscript (`cc_spr_full_paper_claude.tex`) is the base shell.**

Rationale:
- Heaviest mathematics (definitions, propositions, proofs, algorithm box)
- Formal geometry-backend definitions (Euclidean, Ricci, Diffusion, PHATE, DTM) with proper math
- Filtration-perturbation stability proposition with proof
- Sparsity guarantee proposition
- Convexity proposition with full proof
- Algorithm pseudocode
- Computational complexity analysis
- Systematic comparison table (Table 5)
- Related work section with real citations
- Biological interpretation paragraphs for each dataset

## 3. Section-by-Section Merge Decisions

### Introduction
- **Source**: Claude manuscript (Sections 1–1.2)
- **Changes**: Expand the biological motivation paragraph. Add explicit language about omics as geometry-rich, noisy, high-dimensional. Strengthen the "from prediction to interpretation" framing.
- **From Hellbender**: Use the "execution policy" paragraph tone for resource-aware language.

### Background and Relationship to Harmonic PH
- **Source**: Claude manuscript (Section 2)
- **Changes**: Expand into a full section with dedicated subsections. Add the comparison subsection/table with all 8 comparison axes (metric assumption, representative criterion, sparsity, interpretability, stability, computational cost, backend flexibility, biological alignment, implementation status).
- **From Hellbender**: Table structure from Hellbender Table 6 but with expanded axes.

### Mathematical Formulation
- **Source**: Claude manuscript (Section 3)
- **Changes**: This section is already the strongest. Add:
  - Entropy / Gini / support-size formal definitions
  - Representative support concentration remark
  - Remark on bridge-edge / bottleneck-sensitive backends sharpening persistent signal
  - Filtration perturbation and barcode sensitivity remark (already present as Proposition 2)

### Geometry Backends
- **Source**: Claude manuscript (Sections 4.1–4.5)
- **Changes**: Retitle subsections per instructions:
  - "Diffusion geometry, heat kernels, and manifold learning"
  - "Information geometry and PHATE-like distances"
  - "Robust topology via DTM"
  - "Ricci curvature, bridge edges, and bottleneck sharpening"
- Expand each with biological/omics relevance paragraphs.
- Add explicit "implementation status" note per subsection.

### TIP Section
- **Source**: Claude manuscript (Section 5)
- **Keep as-is**: Already strong with Definition and two Remarks.

### Comparison with Harmonic PH
- **Source**: Claude manuscript (Section 6)
- **Changes**: Expand comparison table to include all requested axes. Add "when is harmonic sufficient" and "when does CC-SPR add value" paragraphs (already present).

### Experimental Design
- **Source**: Claude manuscript (Section 7)
- **From Hellbender**: Tighten execution-envelope language. Use Hellbender's wording for the bounded scRNA policy.

### Results
- **Source**: Claude manuscript (Section 8)
- **Changes**: Reorganize by question:
  1. Does sparsity help representative interpretability?
  2. When does geometry help?
  3. When does Ricci help, and when does it not?
  4. What do the topological diagnostics reveal?
  5. What is gained biologically even when predictive performance does not dominate?
- Keep all existing tables and figures.
- Add biological interpretation paragraphs for LUAD subtype structure, GSE161711 disease-aligned interpretation, bounded GSE165087 heterogeneity interpretation.

### Resource-Bounded Single-Cell Section
- **Source**: Hellbender manuscript (Section 6) — strongest language
- **Merge**: Claude manuscript's ablation-safe results (F1 = 0.823, lambda sensitivity)
- Add table of execution envelope (already in both)

### Biological Interpretation Section (NEW)
- **Source**: Strengthen from Claude manuscript's scattered biology paragraphs
- **New content**: Dedicated section on:
  - LUAD subtype structure and topology
  - GSE161711 disease-aligned bulk structure (venetoclax resistance, BCL-2 family)
  - Bounded GSE165087 cellular heterogeneity interpretation
  - Why topology-aware sparse representatives identify concentrated disease-relevant modules
  - How CC-SPR improves interpretability even when predictive gains are modest

### Limitations and Future Directions
- **Source**: Claude manuscript (Section 9)
- **From Hellbender**: Cleaner resource-bounded wording
- **Expand**: Add explicit future-work items for each limitation

### Related Work
- **Source**: Claude manuscript (Section 10)
- **Keep**: Already strong with real citations

### Conclusion
- **Source**: Claude manuscript (Section 11)
- **Strengthen**: Re-emphasize harmonic PH parent framework, geometry-aware backend family, biological interpretability, honest computational scope

## 4. Figures — Main Paper vs Supplement

### Main Paper Figures (promoted)
1. **LUAD anchor**: `f1_bars_ci.png` + `paired_deltas_luad_multiseed.png` (side-by-side)
2. **Geometry/ablation**: `ablation_profiles_luad_gse161711.png`
3. **H1 diagnostics**: `h1_lifetime_eu_vs_ricci.png` + `h1_prominence_eu_vs_ricci_cll_venetoclax_gse161711.png`
4. **TIP stability**: `tip_eu_vs_ricci.png` + `tip_eu_vs_ricci_cll_venetoclax_gse161711.png`
5. **Lambda robustness**: `lambda_robustness.png`
6. **Geometry backend summary**: `geometry_backend_summary.png`
7. **Architecture diagram**: `cc_spr_fullstack_v9.png` (if quality sufficient)

### Main Paper Tables
1. Main benchmark results (Table 1)
2. LUAD multiseed paired statistics (Table 2)
3. Best ablation settings (Table 3)
4. Harmonic-PH vs CC-SPR comparison (Table 4)
5. scRNA execution envelope (Table 5)

### Supplement Figures
- GSE165087-specific diagnostic plots (bounded feasibility)
- BRCA figures if available
- `ablation_best_summary.png`
- Per-dataset TIP/lifetime/prominence plots not in main paper

### Supplement Tables
- Runtime/cost summary from logs
- Full ablation profiles
- Representative diagnostics summary
- Geometry diagnostics summary
- Dataset metadata

## 5. Claims Kept, Softened, or Removed

### Claims Kept (with evidence)
- CC-SPR is a strict extension of harmonic PH (mathematical argument)
- Sparse elastic objective is convex (Proposition + proof)
- Sparsity guarantee (Proposition + proof sketch)
- Filtration stability under backend perturbation (Proposition + proof)
- DTM stability (Proposition, cited)
- LUAD: standard > harmonic > CC-SPR in raw F1 (Table 1)
- GSE161711: CC-SPR > harmonic by +0.018 (Table 1)
- GSE165087 ablation-safe: CC-SPR > harmonic by +0.025
- LUAD: Ricci sparse > Euclidean sparse by +0.046 (not significant, p=0.308)
- Ablation divergence between datasets (Table 3)
- TIP provides stability diagnostic (empirical figures)

### Claims Softened
- "Ricci helps" → "Ricci helps on LUAD but not GSE161711; dataset-dependent"
- "Geometry-aware is better" → "Geometry-awareness exposes dataset-dependent structure"
- Any single-cell claims → "bounded feasibility check, not comprehensive benchmark"

### Claims Removed
- No claims about diffusion/PHATE/DTM empirical performance (not implemented)
- No claims about pathway enrichment (not computed)
- No claims about survival stratification
- No cross-tool harmonic PH comparison (maTilDA parity incomplete)

## 6. Experiments Already Sufficient

- LUAD main benchmark (2-split + 18-seed paired stats) ✓
- GSE161711 5-fold CV with full ablation sweep ✓
- BRCA 2-fold supporting evidence ✓
- GSE165087 manuscript-safe + ablation-safe ✓
- LUAD + GSE161711 full 6-axis ablation sweep ✓
- All diagnostic figures (TIP, H1 lifetime, prominence, lambda robustness) ✓

## 7. Lightweight Experiments Worth Adding

**Priority 1 — LUAD multiseed paired statistics refresh**: Already completed (18 seeds). No new experiment needed.

**Priority 2 — GSE161711 bounded ablation refresh**: Already completed in cached sweeps. No rerun needed.

**Priority 3 — Lightweight diagnostics from existing outputs**:
- Support size, entropy, Gini, TIP concentration: **NOT available from cached outputs** (per `representative_diagnostics_summary.csv`). Would require rerunning topology loops → **skip, note as limitation**.
- H1 lifetime, H1 prominence: Available ✓
- Runtime/memory/cost table from logs: Available ✓

**Priority 4 — Optional tiny scRNA sensitivity**: Already completed via ablation-safe config (max_cells=400, n_pcs=25, lambda sensitivity). No additional run needed.

**Verdict: No new experiments required.** All evidence for the final manuscript is already cached.

## 8. Analyses Too Expensive → Future Work

1. Diffusion/PHATE/DTM backend implementation and benchmarking
2. Full-scale scRNA ablation grids (multi-split, multi-bootstrap)
3. Per-bootstrap representative entropy/support-size recomputation
4. Pathway enrichment analysis
5. Survival stratification
6. maTilDA/R harmonic PH parity testing
7. Larger bootstrap counts for TIP convergence analysis

## 9. Output Files

| File | Status |
|---|---|
| `manuscript_claude/cc_spr_final_merged.tex` | To create |
| `manuscript_claude/cc_spr_final_merged_supplement.tex` | To create |
| `manuscript_claude/cc_spr_final_merged.pdf` | To compile |
| `manuscript_claude/cc_spr_final_merged_supplement.pdf` | To compile |
| `results_claude/final_merge_plan.md` | This file |
| `results_claude/final_manuscript_upgrade_report.md` | To create |
