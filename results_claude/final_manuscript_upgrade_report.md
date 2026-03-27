# CC-SPR Final Manuscript Upgrade Report

**Date**: 2026-03-21
**Output files**:
- `manuscript_claude/cc_spr_final_merged.tex` (main paper)
- `manuscript_claude/cc_spr_final_merged_supplement.tex` (supplement)

---

## 1. Exact Merge Decisions Made

### Base shell
The Claude manuscript (`cc_spr_full_paper_claude.tex`, 871 lines) was used as the base shell because it had the strongest mathematics, most complete geometry-family development, and best biological interpretation.

### Section-by-section decisions

| Section | Source | Changes |
|---|---|---|
| Introduction | Claude manuscript | Expanded with 5 subsections: omics as geometry-rich setting, why PH matters, why harmonic PH is important, limitations, why sparse geometry-aware reps matter biologically. Added explicit "What CC-SPR is not" paragraph. |
| Background & HPH comparison | Claude manuscript | Expanded comparison table to 9 axes. Added "when is harmonic sufficient" and "when does CC-SPR add value" paragraphs. |
| Mathematical formulation | Claude manuscript | Added: chain groups/boundary maps definition, entropy/Gini/support-size formal definitions, support concentration remark, bridge-edge/bottleneck remark, per-bootstrap diagnostics note (honest about what's not cached). |
| Geometry backends | Claude manuscript | Retitled per instructions. Added biology paragraphs for each backend. Added implementation status notes. Expanded Ricci section with bottleneck-sharpening explanation. |
| TIP | Claude manuscript | Kept as-is (already strong). |
| Experimental design | Claude + Hellbender | Used Hellbender resource-policy language. Added three-way separation (completed/bounded/conceptual). |
| Results | Claude manuscript | Reorganized by question (5 subsections). Added "what is gained biologically" subsection. |
| Resource-bounded scRNA | Hellbender manuscript (primary) + Claude | Used Hellbender's failure analysis language. Added ablation-safe results from Claude version. Added code-level correction. |
| Biological interpretation | **NEW section** | Dedicated section with 5 subsections on LUAD subtype structure, GSE161711 disease biology, bounded GSE165087 heterogeneity, concentrated disease-relevant modules, interpretability value independent of prediction. |
| Limitations & future | Claude + Hellbender | Expanded to 6 numbered limitations with explicit future-work items for each. Added broader future directions paragraph. |
| Related work | Claude manuscript | Kept as-is (already has real citations). |
| Conclusion | Claude manuscript | Strengthened with explicit bullet-point summary of what the evaluation confirms. |

### Supplement decisions
- Used Claude supplement as base, expanded with:
  - All scRNA diagnostic figures (manuscript-safe + ablation-safe)
  - GSE161711 extended diagnostics (F1 + lambda robustness)
  - Full ablation details with cross-dataset divergence analysis
  - Cross-dataset pattern interpretation section
  - Complete runtime/resource summary table

## 2. Lightweight Experiments/Ablations Added

**None required.** All evidence for the final manuscript was already cached in the repository:
- LUAD: 2-split benchmark + 18-seed paired stats + full 6-axis ablation ✓
- GSE161711: 5-fold CV + full 6-axis ablation ✓
- BRCA: 2-fold CV ✓
- GSE165087: manuscript-safe + ablation-safe with lambda sensitivity ✓
- All diagnostic figures (TIP, H1 lifetime, prominence, lambda robustness) ✓
- All derived tables (geometry backend summary, ablation profiles, runtime costs) ✓

## 3. Final Manuscript PDF Path

`manuscript_claude/cc_spr_final_merged.pdf`

## 4. Final Supplement PDF Path

`manuscript_claude/cc_spr_final_merged_supplement.pdf`

## 5. Remaining Limitations and Future Work

### Limitations stated in paper
1. Standard baselines win on prediction (shared with parent harmonic PH framework)
2. Backend family exceeds implementation (diffusion/PHATE/DTM not instantiated)
3. Ricci benefit is non-monotone and dataset-dependent
4. Single-cell evidence is resource-bounded (3 failed configs, 2 bounded successes)
5. Public substitutes replace controlled-access datasets
6. Representative diagnostics are partial (no cached entropy/Gini)

### Future work identified
- Implement and benchmark diffusion, PHATE-like, DTM backends
- Scale single-cell evaluation beyond bounded configurations
- Connect TIP-selected gene sets to pathway enrichment (KEGG, Reactome, MSigDB)
- Survival stratification using TIP-derived features
- Cross-tool harmonic PH comparison (maTilDA parity)
- Larger bootstrap counts for TIP convergence analysis
- Per-bootstrap support entropy/Gini caching for post hoc analysis
- Extension to higher homology dimensions ($H_2$+)
- Integration of topology-derived features as inputs to conventional classifiers

## 6. Key Improvements Over Previous Manuscripts

| Improvement | Details |
|---|---|
| Mathematical depth | Added chain group/boundary map definitions, entropy/Gini definitions, support concentration remark, bridge-edge remark |
| Biological interpretation | New dedicated section (5 subsections) with LUAD subtype, GSE161711 disease, GSE165087 heterogeneity interpretation |
| Honest claims | Three-way separation of completed/bounded/conceptual work; explicit implementation status per backend |
| Results organization | Reorganized by scientific question rather than file dump |
| Comparison table | Expanded to 9 axes (was 7) with clearer formatting |
| Introduction structure | 5 subsections building from general motivation to specific contributions |
| Resource transparency | Integrated Hellbender failure analysis language throughout |
| Paper length | ~1100 lines main paper + ~400 lines supplement (unconstrained by page limit) |
