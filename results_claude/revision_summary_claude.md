# Revision Summary (Claude)

Date: 2026-03-20

## Main weaknesses fixed

### 1. Report-like tone → paper-like prose
- Moved Slurm job history table and code-level correction details to supplement
- Replaced with concise 5-line execution envelope paragraph in main paper
- Tightened protocol section with narrative flow

### 2. Biological interpretation added
- **LUAD**: Explained why standard baselines dominate (variance-dominant axes), interpreted TIP as topology-grounded alternative to variance-based feature selection
- **GSE161711**: Added CLL/venetoclax biology (BCL-2 pathway rewiring, metabolic reprogramming), explained why sparse representatives may capture transitional geometry, noted TIP divergence between backends
- **GSE165087**: Integrated ablation-safe results with biological framing of heterogeneity validation

### 3. Ablation-safe GSE165087 results integrated
- Added new bounded config results: standard 0.888, harmonic 0.798, CC-SPR 0.823
- CC-SPR exceeds harmonic by +0.025 — more informative than the manuscript-safe tie
- Lambda sensitivity (0.847 and 0.850) shows stable sparse representatives
- Results integrated into both main paper and supplement

### 4. Harmonic PH comparison strengthened
- Added "What CC-SPR inherits from harmonic PH" paragraph
- Expanded "When does CC-SPR add value?" with three specific scenarios
- Made inheritance vs novelty distinction explicit

### 5. "When Ricci Helps" section sharpened
- Added mechanistic explanation (shortcut edges, curvature correction)
- Connected LUAD benefit to subtype landscape geometry
- Connected GSE161711 non-benefit to smaller cohort size and sampling noise

### 6. Unimplemented backend sections trimmed
- Shortened diffusion, PHATE, DTM sections by ~30%
- Added explicit italic markers: "This backend is not experimentally instantiated"
- Reduced implied promise of parity with implemented backends

### 7. Notation and references fixed
- Fixed `\textsc{Backend}` → `\texttt{Backend}` in Algorithm 1
- Added authors to harmonicPH reference (Damiani, Vipond, Harrington)
- Added DOI for harmonicPH

### 8. Conclusion rewritten
- From list-style summary to forward-looking synthesis
- Explicitly states three natural next steps
- Connects to both methodological and biological future work

### 9. Supplement improved
- Added ablation-safe GSE165087 results section
- Added "Interpretation of Cross-Dataset Patterns" section with three analytical paragraphs
- Supplement grew from 8 to 10 pages with scientific content

## Lightweight analyses added
- None (all improvements are writing-only using existing cached results)

## Final outputs
- Paper: `manuscript_claude/cc_spr_full_paper_claude.pdf` (18 pages, 913 KB)
- Supplement: `manuscript_claude/cc_spr_full_supplement_claude.pdf` (10 pages, 168 KB)
- Gap report: `results_claude/revision_gap_report_claude.md`

## Remaining weaknesses requiring future computation
1. No pathway enrichment analysis on TIP-selected genes (would require enrichment pipeline)
2. Diffusion, PHATE, DTM backends remain unimplemented (would require significant code + compute)
3. Comprehensive single-cell benchmarking (would require larger memory nodes)
4. Per-bootstrap support entropy and Gini diagnostics (would require rerunning topology loops)
5. Author/affiliation fields still empty
