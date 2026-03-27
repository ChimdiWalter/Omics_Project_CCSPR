# CC-SPR Reference Audit — Most Recent Update

**Audit date:** 2026-03-22
**Audited manuscript:** `manuscript_claude/cc_spr_final_geometry_complete.tex`
**Corrected manuscript:** `manuscript_claude/cc_spr_most_recent_update.tex`

---

## Reference-by-Reference Verification

### 1. `\bibitem{harmonicPH}` — **CORRECTED (MAJOR)**

**Old entry:**
> D. Damiani, O. Vipond, and H. A. Harrington. Probing omics data via harmonic persistent homology. *bioRxiv*, 2024. `doi:10.1101/2024.09.15.613148`.

**Verification finding:**
This entry conflates two distinct works:

- **The omics application paper:** Gurnari, D., Guzmán-Sáenz, A., Utro, F., Morovati, V., Lodi, S., and Harrington, H. A. "Probing omics data via harmonic persistent homology." *Scientific Reports* 15, 38836 (2025). arXiv preprint: 2311.06357. Published in the journal — NOT a bioRxiv preprint.
- **The mathematical foundation:** Basu, S. and Cox, N. "Harmonic Persistent Homology." *SIAM Journal on Applied Algebra and Geometry* 8(1):189–224 (2024). arXiv preprint: 2105.15170.

The authors Damiani, Vipond, and Harrington do not appear as authors of "Probing omics data via harmonic persistent homology." The DOI `10.1101/2024.09.15.613148` does not match either paper. This was a hallucinated/conflated reference.

**Action:** Split into two separate references: `\bibitem{gurnariHPH2025}` for the omics application and `\bibitem{basuCoxHPH2024}` for the mathematical foundation. All in-text citations updated to point to the appropriate reference depending on context.

**Authoritative source:** Nature/Scientific Reports page (https://www.nature.com/articles/s41598-025-12189-y); arXiv (https://arxiv.org/abs/2105.15170); SIAM J. Appl. Algebra Geom. (https://epubs.siam.org/doi/10.1137/22M1518761)

**Status:** CORRECTED

---

### 2. `\bibitem{edelsbrunner2010}`

**Entry:** H. Edelsbrunner and J. Harer. *Computational Topology: An Introduction*. AMS, 2010.

**Verification:** Correct. Published by AMS in 2010. ISBN 978-0-8218-4925-5.

**Authoritative source:** AMS Bookstore (https://bookstore.ams.org/view?ProductCode=MBK/69)

**Status:** VERIFIED

---

### 3. `\bibitem{carlsson2009}`

**Entry:** G. Carlsson. Topology and data. *Bulletin of the AMS*, 46(2):255–308, 2009.

**Verification:** Correct. Published in Bull. Amer. Math. Soc. 46(2), 255–308, 2009.

**Authoritative source:** AMS (https://www.ams.org/bull/2009-46-02/S0273-0979-09-01249-X/)

**Status:** VERIFIED

---

### 4. `\bibitem{rizvi2017}` — **CORRECTED (MINOR)**

**Old entry:**
> A. H. Rizvi, P. G. Camara, E. K. Kandber, et al. Single-cell topological RNA-seq analysis reveals insights into cellular differentiation and development. *Nature Biotechnology*, 35(6):551–560, 2017.

**Verification finding:** The third author's name is misspelled. It should be "E. K. Kandror" not "E. K. Kandber."

**Authoritative source:** Nature Biotechnology (https://www.nature.com/articles/nbt.3854); PubMed (https://pubmed.ncbi.nlm.nih.gov/28459448/)

**Action:** Corrected "Kandber" to "Kandror."

**Status:** CORRECTED

---

### 5. `\bibitem{xia2014}`

**Entry:** K. Xia and G.-W. Wei. Persistent homology analysis of protein structure, flexibility, and folding. *Int. J. Numer. Methods Biomed. Eng.*, 30(8):814–844, 2014.

**Verification:** Correct. DOI: 10.1002/cnm.2655.

**Authoritative source:** Wiley Online Library (https://onlinelibrary.wiley.com/doi/abs/10.1002/cnm.2655)

**Status:** VERIFIED

---

### 6. `\bibitem{vipond2021}` — **CORRECTED (MINOR)**

**Old entry:**
> O. Vipond, J. A. Bull, P. S. Sherland, et al. Multiparameter persistent homology landscapes identify immune cell spatial patterns in tumors. *PNAS*, 118(41), 2021.

**Verification finding:** The third author is Philip S. Macklin, not "P. S. Sherland." Full author list: Vipond, Bull, Macklin, Tillmann, Pugh, Byrne, Harrington.

**Authoritative source:** PNAS (https://www.pnas.org/doi/10.1073/pnas.2102166118); PubMed (https://pubmed.ncbi.nlm.nih.gov/34625491/)

**Action:** Corrected "P. S. Sherland" to "P. S. Macklin."

**Status:** CORRECTED

---

### 7. `\bibitem{moon2019phate}`

**Entry:** K. R. Moon, D. van Dijk, Z. Wang, et al. Visualizing structure and transitions in high-dimensional biological data. *Nature Biotechnology*, 37(12):1482–1492, 2019.

**Verification:** Correct. DOI: 10.1038/s41587-019-0336-3.

**Authoritative source:** Nature Biotechnology (https://www.nature.com/articles/s41587-019-0336-3)

**Status:** VERIFIED

---

### 8. `\bibitem{zomorodian2005}`

**Entry:** A. Zomorodian and G. Carlsson. Computing persistent homology. *Discrete & Computational Geometry*, 33(2):249–274, 2005.

**Verification:** Correct. DOI: 10.1007/s00454-004-1146-y.

**Authoritative source:** Springer (https://link.springer.com/article/10.1007/s00454-004-1146-y)

**Status:** VERIFIED

---

### 9. `\bibitem{friedman1998}`

**Entry:** J. Friedman. Computing Betti numbers via combinatorial Laplacians. *Algorithmica*, 21(4):331–346, 1998.

**Verification:** Correct. DOI: 10.1007/PL00009218.

**Authoritative source:** Springer (https://link.springer.com/article/10.1007/PL00009218)

**Status:** VERIFIED

---

### 10. `\bibitem{zou2005elasticnet}`

**Entry:** H. Zou and T. Hastie. Regularization and variable selection via the elastic net. *JRSS Series B*, 67(2):301–320, 2005.

**Verification:** Correct. Full journal name: Journal of the Royal Statistical Society: Series B (Statistical Methodology). DOI: 10.1111/j.1467-9868.2005.00503.x.

**Authoritative source:** Wiley (https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9868.2005.00503.x)

**Status:** VERIFIED

---

### 11. `\bibitem{chazal2014persistence}` — **CORRECTED (MINOR)**

**Old entry:**
> F. Chazal, V. de Silva, M. Glisse, and S. Oudot. The structure and stability of persistence modules. *Springer*, 2016.

**Verification finding:** This is a book (SpringerBriefs in Mathematics, 2016). The bibitem key says "2014" which is misleading (the arXiv preprint is from 2012). The content year 2016 in the text is correct. Added ISBN for completeness.

**Authoritative source:** Springer (https://link.springer.com/book/10.1007/978-3-319-42545-0)

**Action:** Minor formatting improvement. Added "SpringerBriefs in Mathematics" series name.

**Status:** CORRECTED (minor)

---

### 12. `\bibitem{ollivier2007}`

**Entry:** Y. Ollivier. Ricci curvature of metric spaces. *Comptes Rendus Mathematique*, 345(11):643–646, 2007.

**Verification:** Correct. DOI: 10.1016/j.crma.2007.10.041.

**Authoritative source:** Numdam (https://www.numdam.org/item/CRMATH_2007__345_11_643_0/)

**Status:** VERIFIED

---

### 13. `\bibitem{ollivier2009}`

**Entry:** Y. Ollivier. Ricci curvature of Markov chains on metric spaces. *Journal of Functional Analysis*, 256(3):810–864, 2009.

**Verification:** Correct. DOI: 10.1016/j.jfa.2008.11.001.

**Authoritative source:** ScienceDirect (https://www.sciencedirect.com/science/article/pii/S002212360800493X)

**Status:** VERIFIED

---

### 14. `\bibitem{ni2019community}`

**Entry:** C.-C. Ni, Y.-Y. Lin, F. Luo, and J. Gao. Community detection on networks with Ricci flow. *Scientific Reports*, 9:9984, 2019.

**Verification:** Correct. DOI: 10.1038/s41598-019-46380-9.

**Authoritative source:** Nature/Scientific Reports (https://www.nature.com/articles/s41598-019-46380-9)

**Status:** VERIFIED

---

### 15. `\bibitem{coifman2006diffusion}`

**Entry:** R. R. Coifman and S. Lafon. Diffusion maps. *Applied and Computational Harmonic Analysis*, 21(1):5–30, 2006.

**Verification:** Correct. DOI: 10.1016/j.acha.2006.04.006.

**Authoritative source:** ScienceDirect (https://www.sciencedirect.com/science/article/pii/S1063520306000546)

**Status:** VERIFIED

---

### 16. `\bibitem{haghverdi2016}`

**Entry:** L. Haghverdi, M. Büttner, F. A. Wolf, F. Buettner, and F. J. Theis. Diffusion pseudotime robustly reconstructs lineage branching. *Nature Methods*, 13(10):845–848, 2016.

**Verification:** Correct. DOI: 10.1038/nmeth.3971.

**Authoritative source:** Nature Methods (https://www.nature.com/articles/nmeth.3971)

**Status:** VERIFIED

---

### 17. `\bibitem{chazal2011dtm}`

**Entry:** F. Chazal, D. Cohen-Steiner, and Q. Mérigot. Geometric inference for probability measures. *Foundations of Computational Mathematics*, 11(6):733–751, 2011.

**Verification:** Correct. DOI: 10.1007/s10208-011-9098-0.

**Authoritative source:** Springer (https://link.springer.com/article/10.1007/s10208-011-9098-0)

**Status:** VERIFIED

---

### 18. `\bibitem{chazal2017robust}` — **CORRECTED (MINOR)**

**Old entry:**
> F. Chazal, B. Fasy, F. Lecci, B. Michel, A. Rinaldo, and L. Wasserman. Robust topological inference: distance to a measure and kernel distance. *JMLR*, 18(159):1–40, 2017.

**Verification finding:** The paper was published in JMLR volume 18, issue 159, pages 1–40, in **2018**, not 2017.

**Authoritative source:** JMLR (https://jmlr.org/papers/v18/15-484.html)

**Action:** Corrected year from 2017 to 2018.

**Status:** CORRECTED

---

### 19. `\bibitem{gudhi2014}` — **CORRECTED (MINOR)**

**Old entry:**
> The GUDHI Project. GUDHI User and Reference Manual, 2014–2024.

**Verification finding:** There is a formal citable publication: Maria, C., Boissonnat, J.-D., Glisse, M., and Yvinec, M. "The Gudhi Library: Simplicial Complexes and Persistent Homology." In *Mathematical Software — ICMS 2014*, Springer LNCS 8592, pp. 167–174, 2014.

**Authoritative source:** Springer (https://link.springer.com/chapter/10.1007/978-3-662-44199-2_28)

**Action:** Replaced vague manual reference with the formal ICMS 2014 publication.

**Status:** CORRECTED

---

### 20. `\bibitem{giusti2015}`

**Entry:** C. Giusti, R. Ghrist, and D. S. Bassett. Two's company, three (or more) is a simplex. *Journal of Computational Neuroscience*, 41(1):1–14, 2016.

**Verification:** The content is correct — J Comput Neurosci 41:1–14, 2016. DOI: 10.1007/s10827-016-0608-6. The bibitem key says "2015" but the text correctly says 2016. No change to displayed text needed.

**Authoritative source:** Springer (https://link.springer.com/article/10.1007/s10827-016-0608-6)

**Status:** VERIFIED (bibitem key is a misnomer but displayed text is correct)

---

### 21. `\bibitem{chan2013topology}`

**Entry:** J. M. Chan, G. Carlsson, and R. Rabadan. Topology of viral evolution. *PNAS*, 110(46):18566–18571, 2013.

**Verification:** Correct. DOI: 10.1073/pnas.1313480110.

**Authoritative source:** PNAS (https://www.pnas.org/content/110/46/18566)

**Status:** VERIFIED

---

### 22. `\bibitem{chen2011hardness}`

**Entry:** C. Chen and D. Freedman. Hardness results for homology localization. *Discrete & Computational Geometry*, 45(3):425–448, 2011.

**Verification:** Correct. DOI: 10.1007/s00454-010-9322-8.

**Authoritative source:** Springer (https://link.springer.com/article/10.1007/s00454-010-9322-8)

**Status:** VERIFIED

---

### 23. `\bibitem{dey2011optimal}`

**Entry:** T. K. Dey, J. Sun, and Y. Wang. Approximating loops in a shortest homology basis from point data. *SoCG*, 2010.

**Verification:** Correct. Published in Proceedings of the 26th Annual Symposium on Computational Geometry (SoCG '10), pp. 166–175, 2010. The bibitem key says "2011" but the text correctly says "2010."

**Authoritative source:** ACM DL (https://dl.acm.org/doi/10.1145/1810959.1810989)

**Status:** VERIFIED (bibitem key is a misnomer but displayed text is correct)

---

## Summary

| Status | Count | References |
|--------|-------|------------|
| VERIFIED (no change) | 16 | edelsbrunner2010, carlsson2009, xia2014, moon2019phate, zomorodian2005, friedman1998, zou2005elasticnet, ollivier2007, ollivier2009, ni2019community, coifman2006diffusion, haghverdi2016, chazal2011dtm, giusti2015, chan2013topology, chen2011hardness, dey2011optimal |
| CORRECTED (major) | 1 | harmonicPH → split into gurnariHPH2025 + basuCoxHPH2024 |
| CORRECTED (minor) | 5 | rizvi2017 (author typo), vipond2021 (wrong author), chazal2014persistence (formatting), chazal2017robust (year), gudhi2014 (proper citation) |
| REMOVED | 0 | — |
| UNRESOLVED | 0 | — |

## In-Text Citation Corrections Required

The `\cite{harmonicPH}` tag appears in multiple locations and each usage must be reassigned:

1. **Line 81** (Section 1.3 "Why harmonic persistent homology is important"): "Harmonic persistent homology \cite{harmonicPH} addresses this interpretability gap..." — This cites the **omics application** of harmonic PH. Should cite both: `\cite{basuCoxHPH2024, gurnariHPH2025}`.

2. **Line 172** (Section 2.2 "Harmonic representatives"): "...coinciding with the harmonic representative in the Hodge-theoretic sense when the chain complex is equipped with the standard inner product \cite{harmonicPH, friedman1998}." — This cites the **mathematical foundation**. Should be `\cite{basuCoxHPH2024, friedman1998}`.

3. **Line 1222** (Section 8 "Related Work"): "The harmonic PH framework \cite{harmonicPH} introduced the bridge from persistent classes to sample-level representatives..." — This describes the **omics application pipeline**. Should cite both: `\cite{basuCoxHPH2024, gurnariHPH2025}`.

### Body-text wording corrections

1. **Line 81**: The sentence "Harmonic persistent homology addresses this interpretability gap by computing a representative cycle..." is accurate for both the mathematical formulation (Basu & Cox) and the omics application (Gurnari et al.). Citation split is sufficient; no wording change needed.

2. **Line 172**: Citation split is sufficient. The Hodge-theoretic connection is from Basu & Cox.

3. **Line 1222**: Revised wording to distinguish the mathematical foundation from the omics application: "The harmonic PH framework, introduced mathematically by Basu and Cox and applied to omics data by Gurnari et al., provides the bridge..."

4. **Line 1231**: "Using diffusion or PHATE distances as persistent-homology backends is a natural extension that has been explored informally but, to our knowledge, not systematized within a sparse representative framework." — The "to our knowledge" claim is reasonable and can stand; no existing publication systematizes this within a sparse representative framework.
