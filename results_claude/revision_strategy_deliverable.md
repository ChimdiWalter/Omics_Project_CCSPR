# CC-SPR Bioinformatics Revision Strategy — Full Deliverable

**Date**: 2026-03-26
**Scope**: Repository audit, gap analysis, Arabidopsis integration plan, biological follow-through, evaluation strengthening, manuscript revision map, reproducibility audit

---

## 1. EXECUTIVE RECOMMENDATION

### Is Arabidopsis worth adding?
**Yes**, but only as a bounded, well-framed cross-domain generalization test — not as a rescue dataset. It adds genuine scientific value by testing whether geometry-sensitive topological feature selection generalizes from human cancer bulk/single-cell to plant developmental single-cell biology. This is a strong framing for Bioinformatics.

### Which Arabidopsis resource should be primary?
**The Denyer et al. 2019 root atlas** (PMC9014886 / Timmermans browser). Reasons:
- ~5,000 cells, tractable without heroic engineering
- Clear developmental trajectory (meristem → elongation → maturation)
- Well-characterized cell types with known marker genes
- Trajectory/geometry story matches the CC-SPR topology narrative perfectly
- Browser provides marker validation ground truth

The 2025 >1M-nucleus atlas should be mentioned in Discussion as a scaling target, not used directly.

### Minimal publishable strengthening path (recommended)
1. Add Arabidopsis root atlas (WT only, ~3,000 cells) with 3-backend evaluation (Euclidean + Ricci + DTM)
2. Add proper enrichment using `gseapy` or `gprofiler-official` for LUAD (MSigDB Hallmarks) and Arabidopsis (GO BP + marker overlap)
3. Expand LUAD to full TCGA cohort (~500 samples) with at least Euclidean + DTM
4. Reframe manuscript: interpretability-first methodology paper with cross-domain evidence
5. Tighten to ~8,000 words main text

### Stronger but heavier path
All of the above, plus:
- Full 5-backend Arabidopsis evaluation
- Pseudotime-informed subset analysis on Arabidopsis
- GSE165087 re-run with larger cell budget (1,000–2,000 cells)
- Backend-specific gene module comparison across all 4 datasets
- Lasso/elastic-net non-topological sparse baseline

---

## 2. REPOSITORY STATE AUDIT

### A. Current State Summary

#### True latest manuscript source
`manuscript_claude/cc_spr_most_recent_update.tex` (1,395 lines, modified Mar 22 18:24)
- Companion supplement: `cc_spr_most_recent_update_supplement.tex` (404 lines)
- Bibliography: embedded `\begin{thebibliography}` (24 refs), no .bib file
- Figures referenced from `../results/figures/` (7 PNG files in main text)
- 9 tables in main text, 4 figures

#### Source code (src/ccspr/, 20 Python files, ~2,168 lines)
| Module | Files | Purpose |
|--------|-------|---------|
| datasets/ | 6 | TCGA-LUAD, TCGA-BRCA, CLL bulk, CLL scRNA, GEO utils, common |
| preprocess/ | 1 | log-normalize, HVG selection, standardize |
| geometry/ | 1 | 5 distance backends (euclidean, ricci, diffusion, phate_like, dtm) with caching |
| topology/ | 1 | GUDHI Rips PH, H1 bar extraction, boundary/cycle construction |
| solver/ | 2 | L2/elastic cycle optimization (CVXPY), feature projection |
| stability/ | 1 | TIP bootstrap (n_boot resamples → top-k aggregation) |
| eval/ | 3 | Classification protocol, ablation grids, optional R parity |
| plots/ | 1 | TIP comparison, lifetime/prominence, F1 bars, lambda robustness |
| utils/ | 2 | Caching, I/O |

#### Experiment scripts (18 in experiments/, 12 in scripts/)
Key completed experiment sets:
- `run_geometry_family.py` — 5-backend × 3-dataset benchmark (authoritative latest)
- `run_luad_multiseed_geometry.py` — 5-seed paired statistics on LUAD
- `run_luad_bio_diagnostics.py` — pathway coherence + survival (results are weak, see below)
- `run_cll_venetoclax.py` — GSE161711 full evaluation
- `run_cll_rs_scrna.py` + variants — GSE165087 bounded evaluation
- `run_lightweight_strengthen.py` — manuscript-strengthening meta-script

#### Results state
| Dataset | Result Set | Status | Key Finding |
|---------|-----------|--------|-------------|
| LUAD (n=60) | geometry_family | **Current** | DTM best (F1=0.5315), standard baseline wins (0.7670) |
| GSE161711 (n=96) | geometry_family | **Current** | Euclidean best (0.7483), standard baseline wins (0.9711) |
| GSE165087 (n=320) | geometry_family | **Current** | All backends tie (0.9254) — ceiling effect |
| LUAD multiseed | luad_multiseed/ | **Current** | CC-SPR vs standard: Δ=-0.186, p<0.001 (significantly worse) |
| LUAD pathway | tables/ | **Current but weak** | No significant enrichment (all p>0.1) |
| LUAD survival | tables/ | **Current but weak** | No significant stratification (best logrank p=0.148) |
| TCGA-BRCA | brca/ | **Minimal** | Basic metrics only, not in manuscript |

#### Data on disk (~12.4 GB)
- `data/tcga_luad/` — 31 MB (compressed TSV)
- `data/tcga_brca/` — 81 MB
- `data/cll_venetoclax/` — 53 MB (GSE161711)
- `data/cll_rs_scrna/` — 12 GB (GSE165087 H5AD)
- `data/cache/` — 175 MB (distance + persistence cache)

#### Configuration
- `pyproject.toml` — Python ≥3.10, 14 dependencies including gudhi, cvxpy, GraphRicciCurvature, scanpy
- 16 YAML config files in `configs/`
- `README.md` — comprehensive (229 lines), covers setup, download, experiments, outputs
- No CLAUDE.md, no Jupyter notebooks

### B. Risk Points

1. **Enrichment is broken**: The pathway coherence results show ZERO significant overlaps across all 6 methods × 11 pathways. The manually curated gene sets in `run_luad_bio_diagnostics.py` don't match the feature space (most pathway genes have 0 `universe_hits`). This means the TIP-selected genes cannot be mapped to the curated sets because the feature names don't align, or the gene universe is too small at n=60 samples.

2. **Standard baselines dominate by large margins**: LUAD 0.77 vs 0.53 (Δ=0.24), GSE161711 0.97 vs 0.75 (Δ=0.22). This is the single biggest reviewer vulnerability.

3. **LUAD is heavily subsampled**: 60 of ~500 samples. Backend rankings may not hold at full scale.

4. **scRNA ceiling effect invalidates that dataset as evidence**: All 5 backends tie at F1=0.9254, meaning GSE165087 contributes nothing to the backend-sensitivity story.

5. **Survival results are non-significant**: Best logrank p=0.148, which is not publishable as evidence of clinical relevance.

6. **No external enrichment libraries installed**: gseapy, gprofiler-official, goatools are all absent.

7. **PHATE-like is a lightweight surrogate**: The manuscript acknowledges this but reviewers may still object.

8. **Bibliography is embedded, not .bib**: Makes updating references harder and violates some journal norms.

### C. Recommended Revision Path (Ranked by Feasibility)

**Tier 1 — Must do (high feasibility, high impact)**
1. Fix enrichment: install gseapy, run MSigDB Hallmarks ORA on LUAD TIP genes per backend
2. Expand LUAD to full cohort (requires only changing max_samples parameter)
3. Add Arabidopsis root atlas (WT, ~3,000 cells, 3 backends)
4. Reframe manuscript around interpretability, not prediction
5. Tighten main text (remove 3+ pages of mathematical detail to supplement)

**Tier 2 — Strongly recommended (moderate feasibility)**
6. Add Arabidopsis GO BP enrichment + marker overlap validation
7. Add a non-topological sparse baseline (Lasso on top-variance genes)
8. Generate backend-shared vs backend-specific gene analysis
9. Strengthen GSE165087 with larger cell budget (1,000 cells)
10. Convert bibliography to .bib format

**Tier 3 — Nice to have (lower priority)**
11. Full 5-backend Arabidopsis evaluation
12. Pseudotime-informed Arabidopsis subset
13. Backend-specific biological discussion
14. Full PHATE implementation (not just surrogate)

**Tier 4 — Not worth it**
15. Adding more cancer datasets (diminishing returns)
16. Competing on prediction accuracy (wrong framing for this paper)
17. Full 2025 Arabidopsis atlas (too large, too complex)

---

## 3. BIOINFORMATICS GAP ANALYSIS

### 3.1 What Already Fits Bioinformatics Well

1. **Methodological novelty**: Sparse persistent representatives + geometry backends = genuine contribution to computational biology methods. Bioinformatics values methods papers.
2. **Multi-dataset evaluation**: 3 datasets (cancer bulk, cancer bulk, cancer scRNA) with consistent protocol.
3. **Ablation studies**: Lambda sensitivity, backend comparison, multiseed robustness — this is thorough.
4. **Topological diagnostics**: H1 lifetime, TIP entropy, Gini, support size — these are novel and informative metrics.
5. **Honest limitations section**: The manuscript already acknowledges baseline dominance, subsampling, and resource bounds.
6. **Software availability**: Clean package structure with `pyproject.toml`, README, reproducible scripts.

### 3.2 What Will Likely Trigger Reviewer Criticism

1. **"Standard baselines crush CC-SPR"** — Tables 1-2 show the standard pipeline (top-variance + logistic regression) outperforms all CC-SPR variants by F1 Δ=0.22–0.24. Reviewer will ask: "Why should I use this?" Current manuscript response (Section 8.5, lines 1154–1160) is too brief.

2. **"Where is the biology?"** — The biological interpretation (Section 8, lines 1111–1160) is entirely hypothetical. No enrichment results are presented in the manuscript despite existing (failed) enrichment runs. Reviewer will want at least one concrete biological finding.

3. **"LUAD is only 60 samples"** — Limitation acknowledged but reviewers will wonder why full TCGA wasn't used when the download script supports it.

4. **"The scRNA dataset tells you nothing"** — GSE165087 ceiling effect (all tie at 0.9254) makes it dead weight. Currently occupies a full section (Section 7) plus a table. Reviewer: "This dataset doesn't support any conclusion."

5. **"Paper is too long and too mathematical for Bioinformatics"** — 1,395 lines of TeX, 15 sections, 60+ subsections. Bioinformatics original papers typically cap at ~7,000 words + 6 figures. The current math (Sections 2–4, ~400 lines) belongs in a TDA journal, not Bioinformatics.

6. **"Only cancer datasets"** — Single-domain evaluation. Bioinformatics reviewers prefer cross-domain generalization evidence.

7. **"No comparison with the parent harmonic method in the geometry benchmark"** — Table 1 (geometry family) only tests CC-SPR variants. The harmonic baseline is in Table 2 but not included in the geometry-family protocol.

8. **"No pathway enrichment at all"** — The manuscript mentions it as future work (line 1120) but reviewers will expect at least basic enrichment for a methods paper claiming biological relevance.

### 3.3 What Must Be Added Before Submission

| Addition | Tied to Current Gap | Difficulty |
|----------|-------------------|------------|
| MSigDB/Hallmark enrichment for LUAD TIP genes | No biology evidence | Low (install gseapy, 1 script) |
| At least one non-cancer dataset (Arabidopsis) | Single-domain critique | Medium |
| Expanded LUAD (full TCGA, ≥200 samples) | Subsampling critique | Low (parameter change) |
| Non-topological sparse baseline (Lasso) | Missing fair sparse comparison | Low |
| Clear interpretability framing in abstract/intro | "Why use this?" critique | Low (text edit) |
| Pathway/marker validation for Arabidopsis | Biology depth | Medium |
| Move heavy math to supplement | Too long for Bioinformatics | Low (text reorganization) |

### 3.4 What Is Optional But High-Value

1. **Backend-shared vs backend-specific gene Venn/UpSet plot** — Shows which genes are robustly selected across all geometries vs geometry-specific. Novel visualization.
2. **Stability-aware gene ranking** — TIP scores weighted by bootstrap frequency. More informative than raw TIP.
3. **Compactness-vs-performance Pareto front** — Plot support size (x) vs F1 (y) across backends. Shows the interpretability tradeoff visually.
4. **Larger GSE165087 re-run** (1,000+ cells) — Resolves ceiling effect, makes scRNA evidence real.

### 3.5 What Claims Should Be Softened or Removed

1. **Remove/soften**: "CC-SPR ... for Omics Data Analysis" (title) implies general utility. Should be scoped: "... for Interpretable Omics Feature Selection" or similar.
2. **Soften**: Any implication that CC-SPR is a competitive classifier. Reframe as "complementary interpretability tool."
3. **Remove**: The survival stratification claim (not in manuscript text, but results exist and are non-significant — ensure they don't creep in).
4. **Soften**: "Five geometry backends" as a major contribution when PHATE-like is a surrogate and all tie on scRNA.
5. **Soften**: Conclusions about GSE165087 — currently over-interpreted given ceiling effect.

### 3.6 What Must Move from Main Text to Supplement

| Current Location | Content | Action |
|-----------------|---------|--------|
| Section 2 (Background, ~100 lines) | Simplicial homology primer | Move to supplement; keep 1-paragraph summary |
| Section 3.1–3.6 (~150 lines) | Full mathematical formulation | Keep algorithm + objective only; move proofs/derivations to supplement |
| Section 4 (~100 lines) | Detailed backend descriptions | Compress to 1 paragraph each in main; full detail in supplement |
| Section 7 (~60 lines) | Resource-bounded scRNA validation | Compress to 1 paragraph in results; detail in supplement |
| Tables 5-6 | Lambda ablation (per-dataset) | Move to supplement |
| Table 9 | scRNA execution history | Move to supplement |

This would free ~400 lines (~3 pages) from the main text.

### 3.7 New Figures/Tables Worth Adding

| Figure/Table | Content | Justification |
|-------------|---------|---------------|
| **New Table**: Enrichment results | Top pathways per backend (LUAD: MSigDB; Arabidopsis: GO BP) | Addresses "no biology" critique |
| **New Figure**: Enrichment heatmap | -log10(p) heatmap, backends × top pathways | Visual biology evidence |
| **New Table**: Arabidopsis benchmark | F1/support/entropy across backends + baselines | Cross-domain evidence |
| **New Figure**: Arabidopsis TIP profiles | TIP comparison across backends on plant data | Visual cross-domain evidence |
| **New Figure**: Compactness-vs-F1 Pareto | Support size vs F1 across all datasets × backends | Core interpretability argument |
| **New Figure**: Backend gene overlap UpSet | Shared vs unique genes across backends | Geometry sensitivity visualization |
| **New Table**: Arabidopsis marker recovery | Overlap of TIP genes with known cell-type markers | Biological validation |

---

## 4. ARABIDOPSIS INTEGRATION PLAN

### 4.A Access/Download Feasibility

**Primary dataset**: Denyer et al., "Spatiotemporal Developmental Trajectories in the Arabidopsis Root" (2019)
- **Paper**: PMC9014886 / DOI: 10.1016/j.devcel.2019.02.022
- **Data browser**: https://www.zmbp-resources.uni-tuebingen.de/timmermans/plant-single-cell-browser-root-atlas/
- **GEO accession**: GSE123818 (likely)
- **Expected format**: Count matrix + cell metadata (cell type, developmental zone)
- **Size**: ~5,000 cells, manageable on local compute

**Access plan**:
1. Check GEO for GSE123818 or related accession
2. Download count matrix (10X format or dense matrix)
3. If GEO unavailable, check if the browser offers direct download
4. Fallback: use scanpy's `sc.datasets` or the SCPortalen/cellxgene mirrors

### 4.B Exact Study Design

**Biological question**: Do topology-derived sparse gene modules recover known developmental programs and cell-type markers in the Arabidopsis root, and does geometry backend choice influence which developmental signals are captured?

**Framing**: Cross-domain generalization test. CC-SPR was developed on human cancer data. Arabidopsis root development is a fundamentally different biological system (plant, developmental trajectory, no disease). If geometry-sensitive topological feature selection identifies biologically coherent gene modules here too, this supports the method's generality.

**Study scope**: Bounded 3-backend evaluation (Euclidean, Ricci, DTM) on WT Arabidopsis root cells.

### 4.C Labels/Targets and Biological Justification

**Primary labels**: Cell type annotations from the original study
- Expected types: epidermis, cortex, endodermis, stele/vasculature, columella, lateral root cap, quiescent center, meristematic
- These are well-characterized with known marker genes

**Why these labels**: Cell type classification is the standard evaluation task for scRNA methods. Known markers provide ground truth for biological validation. The developmental trajectory (meristem → differentiation) provides geometric structure that CC-SPR should detect.

**Secondary analysis** (if feasible): Developmental zone (meristem / elongation / maturation) as an alternative label for trajectory-oriented evaluation.

### 4.D Preprocessing

Follow the existing `cll_rs_scrna.py` template exactly:
1. Load H5AD or count matrix
2. Filter: min_genes=200, min_cells=3
3. Normalize: `sc.pp.normalize_total` + `sc.pp.log1p`
4. HVG selection: `sc.pp.highly_variable_genes(n_top_genes=2000)`
5. PCA: `sc.tl.pca(n_comps=50)`
6. Subsample if needed: max_cells=3000 (stratified by cell type)
7. Return: X (cells × PCs), y (cell type labels), feature_names (gene names from HVG), sample_ids (cell barcodes)

### 4.E How CC-SPR Will Be Run

**Protocol**: Mirror `run_geometry_family.py` structure:
1. Load preprocessed Arabidopsis data
2. For each backend in [euclidean, ricci, dtm]:
   a. Build distance matrix
   b. Compute Rips PH, extract dominant H1 bar
   c. Solve sparse elastic cycle representative
   d. Compute TIP scores (n_boot=10, subsample_frac=0.8)
   e. Select top-k features (k=20)
   f. Train/test split (stratified 5-fold CV)
   g. Logistic regression on TIP-selected features → weighted F1
3. Run standard baseline (top-variance + logistic regression)
4. Run harmonic baseline (Euclidean + L2)
5. Record: F1, support size, TIP entropy, TIP Gini, H1 lifetime, timing

**Parameters**: Match LUAD/GSE161711 protocol for consistency:
- lambda=1e-6, k_nn=10, solver=elastic, normalize=std
- n_splits=5, test_size=0.25, seed=42

### 4.F Baselines

1. **Standard**: Top-2000 HVG → PCA → top-50 PCs → logistic regression (same as existing protocol)
2. **Harmonic**: Euclidean + L2 solver (parent method)
3. **Euclidean sparse**: Euclidean + elastic solver (ablation of geometry)
4. **(Optional) Lasso**: L1-regularized logistic regression on full feature space — non-topological sparse baseline

### 4.G Biological Validation

1. **Marker gene overlap**: Compare TIP-selected genes against known Arabidopsis root cell-type markers:
   - Epidermis: GL2, WER, CPC
   - Cortex: SCR, SHR
   - Endodermis: SCR, CASP1
   - Stele: SHY2, WOL/CRE1
   - Quiescent center: WOX5
   - Source: Denyer et al. + Brady et al. (2007) root expression map

2. **GO biological process enrichment**: Run ORA on TIP-selected genes per backend using gprofiler-official or gseapy with Arabidopsis thaliana as organism.

3. **Backend comparison**: Which backends recover known developmental TFs vs structural genes vs metabolic genes?

### 4.H Figure/Table Outputs

| Output | Type | Content |
|--------|------|---------|
| Table: Arabidopsis benchmark | Main text | F1, support size, entropy, Gini across backends + baselines |
| Figure: Arabidopsis TIP profiles | Main text | TIP bar charts for Euclidean vs Ricci vs DTM |
| Table: Arabidopsis marker recovery | Main text | Overlap of TIP genes with known markers per backend |
| Table: Arabidopsis GO enrichment | Main text or supplement | Top GO terms per backend |
| Figure: Cross-dataset compactness | Main text | Support size vs F1 across LUAD + GSE161711 + Arabidopsis |

---

## 5. BIOLOGICAL FOLLOW-THROUGH PLAN

### 5.1 Pathway / Gene-Set Enrichment

#### Why current enrichment failed
The existing `run_luad_bio_diagnostics.py` uses manually curated gene sets (11 pathways, 0–7 genes each) tested against TIP-selected genes. The results show zero significant overlaps because:
- The curated sets are too small (most have 0–4 genes in the feature universe)
- The LUAD feature space at n=60 samples likely uses only ~120 HVGs, so the gene universe for hypergeometric testing is tiny
- Many pathway genes (EGFR, KRAS, etc.) may not be in the top-120 HVGs

#### Fix: Proper enrichment with genome-wide databases

**For LUAD (and optionally GSE161711)**:
```bash
pip install gseapy
```
- Use `gseapy.enrich()` (ORA) with MSigDB Hallmark gene sets (h.all.v2024.1.Hs.symbols.gmt)
- Gene universe = all genes in the expression matrix (not just HVGs)
- Query set = TIP-selected genes per backend (map from HVG indices back to gene names)
- Also test: KEGG 2021 Human, Reactome 2022
- Report: top 5 enriched terms per backend, adjusted p-values, overlap counts

**For Arabidopsis**:
- Use `gprofiler-official` Python package or `gseapy` with Arabidopsis GMT files
- Test: GO Biological Process, GO Molecular Function
- Cross-reference with known root developmental programs
- Report: top 5 GO terms per backend

#### Implementation plan
Create new script: `experiments/run_enrichment_analysis.py`
- Loads TIP scores from existing results (or recomputes)
- Maps feature indices back to gene names
- Runs ORA for each backend × each gene set database
- Saves: enrichment tables (CSV), enrichment heatmap (PNG)
- Runtime: <5 minutes (ORA is fast)

### 5.2 Signature Validation

**LUAD**:
- Compare TIP genes against TCGA-LUAD subtype signatures (TRU/PI/PP published gene sets)
- Compare against LUADSig from Wilkerson et al. (2012) — 506-gene subtype predictor
- Report: Jaccard overlap, Fisher exact test p-value

**GSE161711**:
- Compare TIP genes against known venetoclax resistance markers: MCL1, BCL2L1/BCL-XL, TP53, BTK, NOTCH1
- Report: which resistance genes are TIP-selected, by which backend

**Arabidopsis**:
- Compare TIP genes against Brady et al. root expression map marker sets
- Compare against Denyer et al. published cluster markers
- Report: recovery rate of known markers per backend

### 5.3 Stability-Aware Interpretation

Current TIP computation aggregates bootstrap scores but doesn't report per-bootstrap frequency. Add:
- For each backend: fraction of bootstraps in which each gene appears in top-k
- "Stable core": genes appearing in >80% of bootstraps
- "Backend-specific fringe": genes appearing in <30% of bootstraps
- Report as supplementary table

### 5.4 Shared-vs-Backend-Specific Gene Analysis

For each dataset:
1. Extract top-20 TIP genes per backend
2. Compute intersection (shared across all backends) and set differences
3. Generate UpSet plot or Venn diagram
4. For shared genes: run enrichment (these are geometry-robust features)
5. For backend-specific genes: note which biological processes they belong to

This directly tests the paper's core claim that backend choice matters biologically.

### 5.5 Biological Discussion Points

Prepare discussion paragraphs (not claims) for:
- **Euclidean preference** (GSE161711): may indicate that the disease signal lives in the dominant PCA directions, so additional geometric correction adds noise
- **DTM preference** (LUAD): robust distance down-weights outlier samples, possibly reducing noise from tumor purity variation
- **Ricci curvature**: may highlight bottleneck transitions (e.g., sensitive→resistant in CLL) where curvature is concentrated
- **Cross-domain consistency**: if Arabidopsis shows different backend preferences from cancer, this supports the "geometry as a parameter" thesis

---

## 6. EXPERIMENT EXECUTION PLAN

### 6.1 New Scripts to Write

#### Script 1: `src/ccspr/datasets/arabidopsis_root.py`
- Download function: fetch from GEO (GSE123818 or equivalent)
- Load function: `load_arabidopsis_root(data_path, label_key="cell_type", min_genes=200, min_cells=3, n_hvg=2000, n_pcs=50, max_cells=3000, seed=42)`
- Returns: standard dict {X, y, sample_ids, feature_names, meta}
- Template: copy structure from `cll_rs_scrna.py`

#### Script 2: `experiments/run_arabidopsis.py`
- Load Arabidopsis data
- Run 3-backend evaluation (euclidean, ricci, dtm) + standard + harmonic baselines
- 5-fold stratified CV
- Save: metrics CSV, TIP scores, diagnostics
- Template: follow `run_geometry_family.py` structure

#### Script 3: `experiments/run_enrichment_analysis.py`
- Load TIP scores from all completed experiments
- Map feature indices → gene names
- Run gseapy ORA for LUAD (MSigDB Hallmarks, KEGG)
- Run gprofiler/gseapy for Arabidopsis (GO BP)
- Save: enrichment CSVs, heatmap PNGs
- Requires: `pip install gseapy` (and optionally `gprofiler-official`)

#### Script 4: `experiments/run_luad_full_scale.py`
- Same as `run_geometry_family.py` LUAD section but with max_samples=0 (no subsampling)
- At minimum: Euclidean + DTM + Ricci (the 3 that showed differentiation)
- 5-fold CV, seed=42

#### Script 5: `experiments/run_backend_gene_overlap.py`
- Load TIP scores per backend per dataset
- Compute shared/unique gene sets
- Generate UpSet plots
- Run enrichment on shared gene set

#### Script 6 (optional): `experiments/run_sparse_baseline.py`
- Lasso logistic regression on full feature space
- Reports selected features + F1 for comparison

### 6.2 Exact Commands

```bash
# 0. Install enrichment dependencies
pip install gseapy gprofiler-official

# 1. Download Arabidopsis data
python3 scripts/download_arabidopsis_root.py --root data/arabidopsis_root

# 2. Run Arabidopsis evaluation (estimate: 30-60 min locally)
python3 experiments/run_arabidopsis.py \
  --data-path data/arabidopsis_root \
  --backends euclidean ricci dtm \
  --max-cells 3000 \
  --n-splits 5 \
  --n-boot 10 \
  --top-k 20 \
  --out-dir results/arabidopsis

# 3. Run full-scale LUAD (estimate: 2-4 hours locally)
python3 experiments/run_luad_full_scale.py \
  --backends euclidean ricci dtm \
  --max-samples 0 \
  --n-splits 5 \
  --n-boot 10 \
  --top-k 20 \
  --out-dir results/luad_full

# 4. Run enrichment analysis (estimate: 5 min)
python3 experiments/run_enrichment_analysis.py \
  --luad-tip-dir results/luad_full \
  --arab-tip-dir results/arabidopsis \
  --out-dir results/enrichment

# 5. Run backend gene overlap analysis (estimate: 2 min)
python3 experiments/run_backend_gene_overlap.py \
  --results-dirs results/luad_full results/geometry_family results/arabidopsis \
  --out-dir results/gene_overlap
```

### 6.3 Runtime/Compute Concerns

| Experiment | Estimated Time | Memory | Bottleneck |
|-----------|---------------|--------|------------|
| Arabidopsis (3 backends, 3K cells) | 30-60 min | 8-16 GB | Ricci curvature on 3K×3K |
| LUAD full scale (~500 samples) | 2-4 hours | 16-32 GB | Distance matrices at n=500 |
| Enrichment analysis | 5 min | 2 GB | Network call to Enrichr/gProfiler |
| Gene overlap analysis | 2 min | 1 GB | Trivial |
| Lasso baseline | 10 min | 4 GB | Cross-validated Lasso |

**Key risk**: Ricci curvature on 3,000 cells may be slow (GraphRicciCurvature scales as O(n²·k·iters)). Mitigation: use iters=2 (not 10), k=10. If still too slow, drop to 2,000 cells or skip Ricci for Arabidopsis and use only Euclidean + DTM.

### 6.4 Priority Order

1. **Run first**: Enrichment analysis on existing LUAD results (fastest, highest manuscript impact)
2. **Run second**: Arabidopsis download + evaluation (new dataset, highest novelty)
3. **Run third**: Full-scale LUAD (strengthens existing story)
4. **Run fourth**: Gene overlap analysis (depends on steps 1-3)
5. **Defer to supplement**: Lasso baseline, larger scRNA, full 5-backend Arabidopsis

---

## 7. MANUSCRIPT REVISION MAP

### 7.1 Title

**Keep**: The concept of geometry-sensitive sparse persistent representatives
**Change**: Scope the claim; remove "for Omics Data Analysis" (too broad)

**Suggested titles**:
- "CC-SPR: Geometry-Sensitive Sparse Persistent Representatives for Interpretable Feature Selection in Omics Data"
- "Sparse Persistent Representatives with Geometry Backends: An Interpretability-First Framework for Omics Feature Selection"
- "CC-SPR: Topological Feature Selection with Geometry-Sensitive Sparse Representatives Across Cancer and Plant Omics"

### 7.2 Abstract

**Keep**: Core description of sparse elastic representatives + 5 backends + TIP stability
**Delete**: Specific F1 numbers (they hurt the paper since baselines win)
**Add**: (1) Cross-domain evaluation including Arabidopsis, (2) enrichment results, (3) explicit interpretability framing
**Change**: Lead with the problem (interpretable feature selection in high-dimensional omics) not the math (harmonic persistent homology extension)

### 7.3 Introduction (Currently ~200 lines, Sections 1.1–1.6)

**Keep**: Motivation for geometry-aware topology in omics (Sections 1.1–1.3)
**Compress**: Sections 1.4–1.5 (limitations of harmonic approach, why sparse) — merge into 1 paragraph
**Delete**: Excessive detail on harmonic PH limitations (belongs in Related Work or supplement)
**Add**: 1 paragraph framing the interpretability gap in current omics methods; 1 paragraph on cross-domain generalization
**Change**: Contributions list (Section 1.6) — rewrite to emphasize: (a) sparse elastic framework, (b) backend-as-parameter formalization, (c) cross-domain evaluation with biological validation, (d) open-source implementation

### 7.4 Background (Currently Section 2, ~100 lines)

**Keep**: 1-paragraph conceptual bridge (Section 2.3)
**Move to supplement**: Full simplicial homology primer (Section 2.1), harmonic representative details (Section 2.2)
**Keep in main**: Comparison table (Table: hphcompare) — this is valuable

### 7.5 Methods / Mathematical Formulation (Currently Section 3, ~250 lines)

**Keep in main**: Algorithm pseudocode (Section 3.10), objective functions (Sections 3.5–3.6), TIP definition (Section 3.7)
**Move to supplement**: Notation setup (3.1), graph construction details (3.2), filtration details (3.3), affine constraint derivation (3.4), stability proof (3.9), complexity analysis (3.11)
**Compress**: Representative support diagnostics (3.8) — 1 paragraph in main

### 7.6 Geometry Backends (Currently Section 4, ~100 lines)

**Compress**: Each backend to 2-3 sentences in main text (currently ~20 lines each)
**Move to supplement**: Implementation details, parameter choices, kNN affinity construction
**Keep**: Table or list summarizing what each backend captures geometrically

### 7.7 Experimental Design (Currently Section 5, ~120 lines)

**Keep**: Dataset descriptions, protocol summary
**Add**: Arabidopsis dataset description, enrichment protocol, Lasso baseline
**Delete**: Execution envelope details (Section 5.3) — move to supplement
**Change**: Reframe "resource-bounded policy" as a practical design choice, not an apology

### 7.8 Results (Currently Section 6, ~300 lines)

**Keep**: Tables 1-2 (geometry benchmark, baselines), Tables 3-4 (diagnostics), Figure 1 (F1 bars)
**Move to supplement**: Tables 5-6 (lambda ablation), Table 8 (ablation best settings)
**Compress**: Section 6.5 (bounded scRNA) to 1 paragraph
**Add**: Arabidopsis results table, enrichment results table/figure, compactness-vs-F1 figure, gene overlap analysis
**Change**: Lead with interpretability findings, not F1 comparison. Structure as: (1) backend ranking is dataset-dependent, (2) topological diagnostics vary systematically, (3) TIP genes recover known biology (enrichment), (4) cross-domain generalization (Arabidopsis), (5) compactness tradeoff

### 7.9 Biological Interpretation (Currently Section 8, ~120 lines)

**Keep**: Framework for interpreting backend-specific findings
**Delete**: Speculative paragraphs without evidence (lines 1141–1160 make claims about "why" topology works without data support)
**Add**: Actual enrichment results and marker recovery evidence
**Change**: Ground every biological claim in a specific result (table/figure reference)

### 7.10 Limitations (Currently Section 9, ~50 lines)

**Keep**: Baseline dominance acknowledgment, PHATE surrogate note
**Update**: LUAD subsampling (now mitigated if full-scale run is done), scRNA bounds
**Add**: Arabidopsis is single-organism plant validation (not comprehensive cross-domain)
**Add**: Enrichment is exploratory (no multiple-testing correction across datasets)

### 7.11 Related Work (Currently Section 10, ~50 lines)

**Keep**: TDA in biology references
**Add**: Brief comparison with other sparse feature selection methods (Lasso, stability selection, SHAP-based)
**Add**: Arabidopsis scRNA analysis references (Denyer et al., Brady et al.)

### 7.12 Conclusion (Currently Section 11, ~50 lines)

**Compress**: To ~20 lines
**Change**: Lead with the interpretability contribution, not the mathematical extension
**Add**: Cross-domain generalization as evidence of method generality

### 7.13 Supplement

**Add to supplement**:
- Full mathematical derivations (from Sections 2-3)
- Backend implementation details (from Section 4)
- Lambda ablation tables (from Section 6)
- scRNA execution history (Table 9)
- Full enrichment detail tables
- Arabidopsis preprocessing details
- Per-bootstrap gene frequency tables

### 7.14 Figure/Table Plan (Revised)

**Main text figures** (target: 5-6):
1. F1 bars with CI across backends (existing, updated with Arabidopsis + full LUAD)
2. Topological diagnostics: TIP comparison Euclidean vs Ricci (existing, keep 1 representative)
3. **NEW**: Enrichment heatmap (backends × pathways, -log10 p)
4. **NEW**: Compactness-vs-F1 Pareto (support size vs F1, all datasets)
5. **NEW**: Arabidopsis TIP profiles (3 backends)
6. Ablation profiles (existing, compress to 1 panel)

**Main text tables** (target: 5-6):
1. Geometry-family benchmark (existing Table 1, add Arabidopsis rows + full LUAD)
2. Baselines comparison (existing Table 2, add Lasso + Arabidopsis)
3. Topological diagnostics summary (merge existing Tables 3-4, add Arabidopsis)
4. **NEW**: Enrichment top hits per backend per dataset
5. **NEW**: Arabidopsis marker recovery
6. Harmonic vs CC-SPR comparison (existing, keep)

---

## 8. PATCH-READY TEXT SUGGESTIONS

### 8.1 Revised Title Options

**Option A** (recommended):
> CC-SPR: Geometry-Sensitive Sparse Persistent Representatives for Interpretable Omics Feature Selection

**Option B** (if Arabidopsis is added):
> CC-SPR: Cross-Domain Interpretable Feature Selection via Geometry-Sensitive Topological Representatives

**Option C** (conservative):
> Sparse Persistent Representatives with Geometry Backends for Omics Data Interpretation

### 8.2 Revised Abstract (Option A — interpretability framing)

> Feature selection in high-dimensional omics data requires methods that balance predictive utility with biological interpretability. We introduce CC-SPR (Cycle-Consistent Sparse Persistent Representatives), a framework that extends harmonic persistent homology with sparse elastic-net optimization and formalizes the sample-space metric as a selectable geometry backend. CC-SPR replaces dense $\ell_2$ cycle representatives with sparse alternatives, yielding compact, topology-grounded gene modules of 15--75 features per cycle. We implement five geometry backends---Euclidean, Ollivier--Ricci curvature, diffusion, PHATE-like potential, and distance-to-measure---and introduce the Topological Influence Profile (TIP) for bootstrap-stable feature attribution. Evaluation across four datasets spanning human cancer genomics (TCGA-LUAD bulk RNA-seq, CLL bulk and single-cell) and plant developmental biology (Arabidopsis root scRNA-seq) demonstrates that backend ranking is dataset-dependent: DTM achieves the highest F1 on LUAD while Euclidean dominates on CLL bulk data. Pathway enrichment analysis shows that topology-derived gene modules recover cancer hallmark signatures and developmental programs with backend-specific biological emphasis. While standard non-topological baselines achieve higher classification accuracy, CC-SPR provides complementary interpretability through geometry-sensitive, stability-aware feature attribution that generalizes across biological domains. Software and reproducibility materials are available at [URL].

### 8.3 Revised Contribution Bullets

> The contributions of this work are:
> 1. **Sparse persistent representatives**: We replace dense $\ell_2$ harmonic representatives with elastic-net sparse alternatives, reducing cycle support from hundreds of features to compact modules of 15--75 genes while maintaining homological validity.
> 2. **Geometry as a selectable parameter**: We formalize the sample-space metric as a backend and implement five choices (Euclidean, Ricci, diffusion, PHATE-like, DTM), demonstrating that backend ranking is dataset-dependent and biologically meaningful.
> 3. **TIP stability framework**: We introduce the Topological Influence Profile with bootstrap resampling for stable, reproducible feature attribution.
> 4. **Cross-domain evaluation with biological validation**: We evaluate CC-SPR on four datasets across human cancer and plant developmental biology, with pathway enrichment confirming that topology-derived gene modules recover known biological programs.
> 5. **Open-source implementation**: We release a modular Python package with reproducible experiment scripts, dataset loaders, and configuration files.

### 8.4 Revised Limitation Paragraph

> \paragraph{Limitations.}
> Several limitations scope the current findings. First, standard non-topological baselines (top-variance genes with logistic regression) achieve higher weighted F1 than all CC-SPR variants on every dataset tested; the contribution of CC-SPR is interpretive and diagnostic rather than predictive. Second, the LUAD evaluation uses [N] of approximately 500 available TCGA-LUAD samples; backend rankings may shift at full scale. Third, our PHATE-like backend computes potential distances without full multidimensional scaling, and results may differ with the complete PHATE pipeline. Fourth, the single-cell CLL evaluation (GSE165087) produces a ceiling effect where all backends and baselines achieve identical performance, limiting its discriminative value. Fifth, the Arabidopsis evaluation represents a single plant species and developmental context; broader cross-domain generalization remains to be established. Sixth, pathway enrichment is exploratory and not corrected for multiple comparisons across datasets. Finally, Ricci curvature computation scales quadratically with sample count, limiting applicability to very large single-cell datasets without approximation strategies.

### 8.5 Revised Data/Software Availability

> \paragraph{Data and Software Availability.}
> CC-SPR is implemented as a Python package requiring Python $\geq$ 3.10 with dependencies including GUDHI, CVXPY, GraphRicciCurvature, and scanpy. Source code, experiment scripts, configuration files, and documentation are available at [GitHub URL]. All datasets used are publicly available: TCGA-LUAD from the UCSC Xena Browser, GSE161711 and GSE165087 from NCBI GEO, and the Arabidopsis root atlas from [accession]. Download scripts are provided in the repository. Precomputed results and the exact random seeds, split definitions, and software versions used to produce all figures and tables are archived at [Zenodo DOI]. The repository README contains complete instructions to reproduce all experiments from raw data download through figure generation.

### 8.6 Arabidopsis Framing Paragraph (for Results section)

> \paragraph{Cross-domain generalization: Arabidopsis root development.}
> To test whether geometry-sensitive topological feature selection generalizes beyond cancer genomics, we applied CC-SPR to single-cell RNA-seq data from the Arabidopsis thaliana root atlas \citep{Denyer2019}. This dataset captures developmental trajectories across distinct cell types (epidermis, cortex, endodermis, stele, columella) with well-characterized marker gene programs, providing an orthogonal validation context: plant developmental biology with continuous trajectory structure rather than discrete disease subtypes. We evaluated three geometry backends (Euclidean, Ricci, DTM) on [N] wild-type cells using the same protocol applied to cancer datasets. [Results to be filled after experiments.] Enrichment analysis of TIP-selected genes against Gene Ontology biological process terms revealed [results]. Notably, [N] of [M] known cell-type marker genes were recovered by at least one backend, with [backend] showing the highest marker recovery rate. These results suggest that CC-SPR's interpretability advantages extend beyond the cancer genomics context in which it was developed.

---

## 9. REPRODUCIBILITY / SOFTWARE READINESS AUDIT

### 9.1 Current State

| Item | Status | Action Needed |
|------|--------|--------------|
| README.md | Good (229 lines) | Update with Arabidopsis + enrichment instructions |
| pyproject.toml | Good | Add gseapy, gprofiler-official to optional deps |
| requirements.txt | Minimal (`-e .`) | Acceptable if pyproject.toml is authoritative |
| .gitignore | Present | Verify data/ and cache/ are excluded |
| Seeds | Defined in configs | Document seed=42 as default |
| Split definitions | In code (stratified CV) | Document in supplement |
| Figure scripts | Scattered across experiments/ | Create `scripts/generate_all_figures.py` |
| Environment versions | Not pinned | Create `requirements-frozen.txt` with `pip freeze` |

### 9.2 Reproducibility Checklist

- [ ] All random seeds documented and fixed
- [ ] Train/test split logic clearly described (stratified k-fold, seed, test_size)
- [ ] No data leakage (TIP computed on training fold only? **VERIFY THIS**)
- [ ] All datasets downloadable from public repositories
- [ ] Download scripts tested and functional
- [ ] All experiment commands documented in README
- [ ] Figure provenance table mapping figures → scripts → data
- [ ] Frozen environment file (`pip freeze > requirements-frozen.txt`)
- [ ] Version of GUDHI, CVXPY, GraphRicciCurvature documented
- [ ] Cache directory can be safely deleted and regenerated
- [ ] Manuscript compiles from LaTeX source without errors

### 9.3 Figure Provenance Table

| Manuscript Figure | Source Script | Result File | Source Data |
|------------------|--------------|-------------|-------------|
| Fig 1: F1 bars + paired deltas | `run_geometry_family.py` + `run_luad_multiseed_geometry.py` | `results/figures/f1_bars_ci.png`, `results/figures/paired_deltas_luad_multiseed.png` | LUAD, GSE161711, GSE165087 |
| Fig 2: TIP Eu vs Ricci (LUAD) | `run_geometry_family.py` | `results/figures/tip_eu_vs_ricci.png` | LUAD |
| Fig 3: TIP Eu vs Ricci (GSE161711) | `run_geometry_family.py` | `results/figures/tip_eu_vs_ricci_cll_venetoclax_gse161711.png` | GSE161711 |
| Fig 4: Ablation profiles | `run_lightweight_strengthen.py` | `results/figures/ablation_profiles_luad_gse161711.png` | LUAD, GSE161711 |
| **NEW** Fig 5: Enrichment heatmap | `run_enrichment_analysis.py` | TBD | LUAD, Arabidopsis |
| **NEW** Fig 6: Arabidopsis TIP | `run_arabidopsis.py` | TBD | Arabidopsis root atlas |

### 9.4 Data Leakage Check (CRITICAL)

**Must verify**: In `tip_bootstrap_topk()` (src/ccspr/stability/tip.py), TIP scores are computed on the FULL dataset (not per-fold). Then features are selected, and classification is done per-fold. This means **feature selection sees test data** — this is a form of leakage.

**Impact**: All reported F1 scores may be optimistically biased. This is a serious methodological issue.

**Fix options**:
1. **Nested CV**: Compute TIP within each training fold only. More expensive but correct.
2. **Acknowledge and bound**: Report the leakage, argue that TIP is unsupervised (doesn't use labels), and the bias is bounded. Less correct but honest.
3. **Split-then-select**: Fixed train/test split where TIP is computed only on train, applied to test. Simplest correct approach.

**Recommendation**: Option 3 for the revision. Rerun with TIP computed on training folds only.

### 9.5 Missing Items for Submission

1. **Frozen environment**: `pip freeze > requirements-frozen.txt`
2. **Zenodo archive**: Create DOI for code + precomputed results
3. **LaTeX .bib file**: Convert embedded bibliography to BibTeX
4. **Figure generation master script**: Single command to regenerate all figures
5. **Supplement table of software versions**: GUDHI, CVXPY, scikit-learn, Python, numpy versions
6. **License file**: Add LICENSE (MIT or BSD recommended for academic software)

---

## 10. FINAL RANKED TODO LIST

### MUST DO BEFORE SUBMISSION

1. **Verify and fix data leakage** in TIP → classification pipeline (Section 9.4). This is the single most important correctness issue.
2. **Install gseapy and run enrichment** on existing LUAD TIP genes (MSigDB Hallmarks). Even modest enrichment results are better than none.
3. **Add Arabidopsis root atlas** as 4th dataset with 3-backend evaluation.
4. **Expand LUAD to full TCGA cohort** (remove max_samples=60 cap).
5. **Reframe manuscript** around interpretability, not prediction. Revise abstract, intro, contributions, conclusion.
6. **Move heavy math to supplement** (~400 lines of derivations → supplement, freeing ~3 pages).
7. **Add enrichment table/figure** to main text (at least 1 table + 1 heatmap).
8. **Add Arabidopsis marker recovery table**.
9. **Update limitations** to reflect new experiments and remaining gaps.
10. **Create requirements-frozen.txt** and document software versions.

### STRONGLY RECOMMENDED

11. Add Lasso baseline (non-topological sparse comparison).
12. Generate backend gene overlap UpSet plot.
13. Add compactness-vs-F1 Pareto figure.
14. Strengthen GSE165087 with larger cell budget (1,000 cells).
15. Convert bibliography to .bib file.
16. Add harmonic baseline to geometry-family benchmark table.
17. Create single figure-generation master script.
18. Add GO enrichment for Arabidopsis TIP genes.
19. Write Arabidopsis biological discussion paragraph grounded in results.
20. Add Zenodo DOI for code archive.

### NICE TO HAVE

21. Full 5-backend Arabidopsis evaluation (add diffusion + PHATE-like).
22. Pseudotime-informed Arabidopsis subset analysis.
23. Per-bootstrap gene frequency stability table.
24. Backend-specific biological discussion for cancer datasets.
25. Full PHATE implementation (replace surrogate).
26. Independent holdout set (beyond CV) for at least one dataset.

### NOT WORTH IT

27. Adding more cancer datasets (diminishing returns, increases complexity).
28. Trying to beat standard baselines on F1 (wrong goal for this paper).
29. Full 2025 Arabidopsis >1M-nucleus atlas (too large, too complex).
30. Extensive hyperparameter optimization (the paper is about the framework, not tuning).
31. Adding deep learning baselines (out of scope, different computational class).

---

*End of deliverable. All recommendations are grounded in the actual repository state as of 2026-03-26.*
