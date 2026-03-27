# CC-SPR Explainer for Geneticists: Production Notes

## Purpose

This document records decisions, sources, and caveats for the standalone explainer
`manuscript_claude/cc_spr_explainer_for_geneticists.tex`.

## Source material

All content in the explainer is derived from a single source:

- **Main manuscript**: `manuscript_claude/cc_spr_most_recent_update.tex` (1395 lines)
- **Compiled PDF**: `manuscript_claude/cc_spr_most_recent_update.pdf`

No external papers were read beyond references cited within the manuscript itself.
The Gurnari et al. (2025) paper is described based on the manuscript's own summary
(Sections 1.3, 2, and the comparison table).

## Sections created

| # | Section | Content source |
|---|---------|---------------|
| 1 | Executive Summary | Abstract + Section 1.6 (Contributions) + Section 1.7 (What CC-SPR is not) |
| 2 | Background: What Gurnari et al. Did | Sections 1.2-1.3, 2.1-2.2 |
| 3 | What Problem Remained | Sections 1.4-1.5 (Limitations 1-2) |
| 4 | What CC-SPR Is Trying to Solve | Section 1.5-1.6 |
| 5 | Intuitive Explanation of the Method | Algorithm 1, Section 3 pipeline description |
| 6 | Minimal Math Intuition | Sections 2.2, 2.5-2.6, 3.5 (simplified) |
| 7 | Geometry Backends for Biologists | Sections 4.1-4.5, Tables 5-6 diagnostics |
| 8 | Datasets and Biological Questions | Section 5.1 |
| 9 | What Experiments Were Done | Sections 5.2-5.4, 7 |
| 10 | What the Ablations Mean | Section 6.3, Tables 7-8 |
| 11 | Main Results | Section 6 (Tables 3-4, 9-10) |
| 12 | Biological Interpretation | Section 8 |
| 13 | Comparison to the Gurnari Paper | Table 1 + Section 2.3-2.4 |
| 14 | Limitations | Section 9 |
| 15 | Why This Matters | Section 8.4-8.5, Conclusion |
| A1 | Glossary | Definitions from throughout the paper |
| A2 | Acronyms | Collected from manuscript |
| A3 | Math-to-Biology Cheat Sheet | Original translation table |

## Figures reused

No figures from the manuscript were included in the explainer.
The explainer is a text-only document to avoid dependency on figure files and to keep
the document self-contained as a standalone PDF. Figures are described in text where
relevant (e.g., TIP profiles, ablation plots).

## Key editorial decisions

1. **No new biology was invented.** All biological interpretations (LUAD subtypes,
   venetoclax resistance mechanisms, CLL heterogeneity) come directly from the
   manuscript's Section 8.

2. **Math was minimized.** The "Minimal Math Intuition" section uses analogies
   (rubber band, stability selection) rather than equations. Formal definitions
   were moved to the glossary.

3. **Honest about limitations.** The explainer repeats all six limitations from the
   manuscript and adds context for biologists about what "methods paper stage" means.

4. **No overhyping.** The explainer explicitly states that standard baselines win on
   prediction, that single-cell evidence is bounded, and that pathway enrichment
   has not been performed.

5. **Gurnari comparison framed as extension.** The explainer emphasizes that CC-SPR
   is a strict extension of the Gurnari framework (special case at lambda_1=0,
   Euclidean), not a competing method.

## Key numbers referenced (all from manuscript tables)

- LUAD standard baseline F1: 0.767
- LUAD best CC-SPR (DTM): 0.532
- GSE161711 standard baseline F1: 0.971
- GSE161711 best CC-SPR (Euclidean): 0.748
- GSE165087 all backends: 0.925 (identical)
- DTM projected support on GSE161711: 24 features
- Euclidean projected support on GSE161711: 59 features
- LUAD multiseed paired delta (Ricci vs Euclidean): +0.046, p=0.308
- LUAD ablation best top-K=10: F1=0.667
- GSE161711 ablation best top-K=50: F1=0.806

## Caveats

1. The Gurnari et al. paper itself was not read independently; all descriptions
   of their work come from the CC-SPR manuscript's characterization.

2. The PHATE-like backend description includes a caveat about the lightweight
   implementation, as stated in the manuscript.

3. No independent verification of the numerical results was performed; all
   numbers are taken directly from the manuscript tables.

4. The biological context (venetoclax resistance mechanisms, LUAD subtype names)
   is drawn from the manuscript and general domain knowledge; no new literature
   review was conducted.
