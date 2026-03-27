# Full-Scale LUAD Run Log — 2026-03-27

## Hardware
- CPU: Intel i9-9920X @ 3.50GHz, 24 threads
- RAM: 125GB
- GPU: 2x RTX 2080 Ti (not used — CVXPY is CPU-bound)

## Dataset
- TCGA-LUAD: 275 samples, 20530 features, 3 classes (Bronchioid/Magnoid/Squamoid)
- Variable genes: 2000
- Standard baseline F1 (full n=275): 0.7501 ± 0.0294

## Run 1: n=200 subsampled, 3-fold CV, 5 bootstraps (SUPERSEDED)

These ran first but were killed/superseded by full n=275 runs.

| Backend | F1 | Std | n_samples | n_splits | n_boot | Status |
|---------|-----|-----|-----------|----------|--------|--------|
| standard | 0.7298 | 0.0761 | 200 | 3 | — | ✓ complete |
| diffusion | 0.5535 | 0.0462 | 200 | 3 | 5 | ✓ complete |
| ricci | 0.6169 | 0.0361 | 119 | 3 | 3 | ✓ complete |
| phate_like | 0.5643 | 0.0432 | 200 | 3 | 5 | ✓ complete |
| euclidean | — | — | 200 | 3 | 5 | ✗ killed mid-run |
| dtm | — | — | 200 | 3 | 5 | ✗ never started |

### TIP Top Genes (n=200 runs)
- diffusion: SYCE1, CYP1A1, MAGEA8, MAGEA10, PASD1
- ricci (n=119): CALB1, KLK14, CALML3, UGT3A1, MAGEB2
- phate_like: NOVA1, CRISP3, CELF3, FAM40B, ORM1

## Run 2: n=275 full scale, 5-fold CV, 5 bootstraps (ACTIVE)

All 5 backends launched 2026-03-27 ~15:59 UTC.

| Backend | Process | Threads | Nice | Task ID | Status |
|---------|---------|---------|------|---------|--------|
| euclidean + dtm | 1 | 8 | 10 | btuf9c98c | running |
| ricci | 2 | 6 | 12 | bzfg256ug | running |
| diffusion | 3 | 5 | 14 | boyeh86z5 | running |
| phate_like | 4 | 5 | 14 | b5z1vxas2 | running |

Parameters: --max-samples 0 --n-boot 5 --n-splits 5 --var-genes 2000 --top-k 20

### Output directories
- euclidean + dtm: results/luad_full/
- ricci: results/luad_full_ricci/
- diffusion: results/luad_full_diffusion/
- phate_like: results/luad_full_phate/

## Manuscript Status

File: manuscript_claude/cc_spr_revised_full.tex (DO NOT OVERWRITE — create new version)

### Already done (previous conversations):
- Arabidopsis section updated with TIP gene results (euclidean + dtm)
- Table tab:arab_genes added (top 5 TIP genes per backend)
- Cross-domain consistency section expanded
- All 25 references verified and corrected:
  - harmonicPH: Gurnari/Guzman-Saenz/Utro/Bose/Basu/Parida, Sci Reports 2025
  - rizvi2017: "Kandber" → "Kandror"
  - vipond2021: "Sherland" → "Macklin"
  - gurnari2024 duplicate removed

## Run 2 Final Results (n=275, all complete)

| Backend | F1 | Std | Top TIP Genes | Walltime |
|---------|-----|-----|---------------|----------|
| standard | 0.7501 | 0.0294 | — | — |
| ricci | 0.6339 | 0.0503 | DEFB4A, LOC254559, SPRR3, TMED6, PAK3 | ~2.5h |
| dtm | 0.5814 | 0.0611 | TDRD9, LOC220594, GJB7, MKRN3, MMP23A | ~1.8h |
| euclidean | 0.5655 | 0.0361 | GDF5, PCDHA6, GABRB3, UCA1, CDKL2 | ~42min |
| diffusion | 0.5522 | 0.0609 | MAGEA10, LBP, PNMA5, PAGE2, CTAG1B | ~2.2h |
| phate_like | 0.5445 | 0.0640 | ATP12A, SST, L1CAM, LOC100190940, MAGEA8 | ~35min |

## Manuscript Updates (COMPLETED 2026-03-27)
- ✓ Added full-scale LUAD row to geobench table (Table 1)
- ✓ Added full-scale LUAD ranking paragraph with TIP gene descriptions
- ✓ Added Table: full-scale TIP genes (tab:fullscale_genes)
- ✓ Updated abstract with scale-sensitive ranking finding
- ✓ Updated datasets table and description (pilot + full)
- ✓ Updated replication protocol (added LUAD full line)
- ✓ Updated limitations: replaced "n=60 may change" with "rankings are scale-sensitive"
- ✓ Updated Ricci and DTM backend discussion paragraphs
- ✓ Updated conclusion with full-scale numbers
- ✓ All references verified (0 orphans)
- ✓ Compiles cleanly: 36 pages PDF
