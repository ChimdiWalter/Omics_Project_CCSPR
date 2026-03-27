# Manuscript Upgrade Report (Claude)

Date: 2026-03-20

## Upgrades relative to the existing manuscript

### Mathematical depth
- Added formal Definition environments (persistence barcode, geometry backend, diffusion distance, PHATE potential distance, DTM, TIP)
- Added Proposition with proof (convexity, sparsity guarantee, filtration stability, DTM stability)
- Added Remark environments (dense support, elastic net connection, non-monotonicity, multiscale property, information geometry, TIP stability, TIP across geometries)
- Added Algorithm pseudocode (Algorithm 1: CC-SPR Pipeline)
- Added computational complexity analysis
- Added connection to elastic net regression literature

### Geometry backend sections (substantially expanded)
- Diffusion geometry: formal definition, heat kernel connection, multiscale property, omics relevance (~1 page)
- PHATE-like distances: formal definition, information-geometric interpretation, omics relevance (~0.5 page)
- DTM: formal definition, stability proposition, DTM-filtration construction, robustness rationale (~0.75 page)
- Ricci: expanded with formal Ollivier definition, flow update equation, non-monotonicity remark (~0.75 page)

### Harmonic PH comparison
- Expanded comparison table (7 columns including implementation status)
- Added paragraph on when harmonic representatives are sufficient
- More precise scope statements

### Results and interpretation
- Added BRCA to main benchmark table
- Improved figure captions with specific quantitative interpretations
- Expanded cross-dataset ablation divergence discussion
- Dedicated subsection on when Ricci helps vs does not

### Bibliography
- Expanded from 3 entries to 23 proper references
- Covers: TDA foundations, biology applications, Ricci curvature, diffusion maps, PHATE, DTM, elastic net, sparse topology

### Structure
- Added formal Background section (Section 2)
- Separated Related Work into its own section
- Better section hierarchy and flow
- Added table of contents to supplement

## Page count
- Main paper: 18 pages (up from ~10)
- Supplement: 8 pages (up from ~5)

## No new computation was performed
All numerical results, figures, and tables are reused from existing repository outputs.
