# Revision Gap Report (Claude)

Date: 2026-03-20
Target: `manuscript_claude/cc_spr_full_paper_claude.tex`

## 1. Sections that sound report-like rather than paper-like

- **Section 5.2 "Execution envelope and resource policy"**: The Slurm job table and "Code-level correction" paragraph read like an engineering log. The job failure table should be condensed or moved to supplement; the code fix paragraph should be removed from the main paper entirely.
- **Section 5.3 "Protocol"**: Very terse and list-like. Needs narrative flow.
- **Section 6.6 "Resource-bounded single-cell validation"**: Reads as a progress report ("does not crash the cluster"). Needs scientific framing.
- **Conclusion**: Reads as a summary list rather than a forward-looking synthesis.

## 2. Claims broader than completed experiments

- **Abstract line**: "We evaluate CC-SPR on four cancer genomics datasets" — technically correct but the four datasets have very different evidence levels (18-pair anchor vs single bounded split). The abstract should signal this asymmetry.
- **Backend sections (diffusion, PHATE, DTM)**: Each is ~0.5–1 page of mathematical development for backends with zero empirical results. This length creates an implied promise. Should be shortened and clearly scoped.
- **Comparison table (Table 2)**: Lists 6 methods including 3 unimplemented ones. The "extension" label is clear but the visual weight implies parity.

## 3. Weak biological interpretation

- **LUAD**: No discussion of what genes or pathways the TIP profiles highlight. No biological hypothesis about why topology-derived features underperform — just "standard baselines win."
- **GSE161711**: No CLL biology. No discussion of venetoclax response biology or what the topological signal might correspond to.
- **GSE165087**: No heterogeneity interpretation. No discussion of what cell populations or transcriptional programs the sparse representatives capture.
- **General**: The paper treats datasets as benchmark slots rather than biological settings. A methods paper still needs qualitative biological interpretation to be convincing.

## 4. Notation / formatting / citation issues

- **Algorithm 1, line 1**: `\textsc{Backend}` — `\textsc` is not a valid LaTeX command; should be `\textsc` or `\text{Backend}` or use `\texttt`.
- **harmonicPH reference**: No authors listed — just "Probing omics data via harmonic persistent homology. bioRxiv / preprint, 2024."
- **Inconsistent ℓ notation**: Uses `$\ell_2$` in text but `\|.\|_2` in math. Should be consistent.
- **Missing `\label` on some tables** in supplement.
- **Author field**: Empty — `\author{}`.

## 5. Lightweight analyses that could help

- **Integrate ablation-safe GSE165087 results**: The run completed (standard 0.888, harmonic 0.798, CC-SPR 0.823) but is not in the manuscript. This is more informative than the manuscript-safe tie result.
- **Add a lightweight runtime comparison table**: Already available from `results/tables/geometry_runtime_cost_summary.csv`.

## 6. Figure/table promotion decisions

- **Promote**: ablation-safe GSE165087 results into main text (more informative than the tie)
- **Move to supplement**: Slurm job history table (Table 3), code-level correction paragraph
- **Keep in main**: All diagnostic figures, ablation table, paired stats table, comparison table
- **Supplement improvement**: Convert from file-inventory style to structured scientific supplement with interpretation

## 7. Summary of planned revisions

1. Reframe paper as methods/interpretability contribution, not benchmark paper
2. Move engineering details to supplement
3. Add qualitative biological interpretation for each dataset
4. Integrate ablation-safe GSE165087 results (CC-SPR 0.823 > harmonic 0.798)
5. Shorten unimplemented backend sections
6. Sharpen harmonic PH comparison
7. Fix notation and references
8. Improve conclusion
