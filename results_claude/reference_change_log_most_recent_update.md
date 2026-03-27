# CC-SPR Reference Change Log — Most Recent Update

**Date:** 2026-03-22
**Audited manuscript:** `manuscript_claude/cc_spr_final_geometry_complete.tex`
**Corrected manuscript:** `manuscript_claude/cc_spr_most_recent_update.tex`

---

## Change 1 (MAJOR): Split `harmonicPH` into two correct references

### Old reference text
```latex
\bibitem{harmonicPH}
D.~Damiani, O.~Vipond, and H.~A.~Harrington.
Probing omics data via harmonic persistent homology.
\emph{bioRxiv}, 2024.
\texttt{doi:10.1101/2024.09.15.613148}.
```

### New corrected references
```latex
\bibitem{basuCoxHPH2024}
S.~Basu and N.~Cox.
Harmonic persistent homology.
\emph{SIAM Journal on Applied Algebra and Geometry}, 8(1):189--224, 2024.

\bibitem{gurnariHPH2025}
D.~Gurnari, A.~Guzm\'{a}n-S\'{a}enz, F.~Utro, V.~Morovati, S.~Lodi, and H.~A.~Harrington.
Probing omics data via harmonic persistent homology.
\emph{Scientific Reports}, 15:38836, 2025.
```

### Reason
The old entry was a hallucinated reference: the authors (Damiani, Vipond, Harrington), venue (bioRxiv), year (2024), and DOI were all incorrect. "Probing omics data via harmonic persistent homology" was written by Gurnari et al. and published in Scientific Reports (2025). The mathematical foundation of harmonic persistent homology is by Basu & Cox in SIAM J. Appl. Algebra Geom. (2024). These are two distinct works that must be cited separately.

### Sections affected
- Section 1.3 (line ~81): `\cite{harmonicPH}` → `\cite{basuCoxHPH2024, gurnariHPH2025}`
- Section 2.2 (line ~172): `\cite{harmonicPH, friedman1998}` → `\cite{basuCoxHPH2024, friedman1998}`
- Section 8 "Related Work" (line ~1222): `\cite{harmonicPH}` → `\cite{basuCoxHPH2024, gurnariHPH2025}`
- Body wording in Section 8 revised to distinguish the two works

---

## Change 2 (MINOR): Fix author typo in `rizvi2017`

### Old reference text
```
A.~H.~Rizvi, P.~G.~Camara, E.~K.~Kandber, et al.
```

### New corrected reference text
```
A.~H.~Rizvi, P.~G.~Camara, E.~K.~Kandror, et al.
```

### Reason
"Kandber" is a misspelling of "Kandror" (Elena K. Kandror). Verified via PubMed and Nature Biotechnology.

### Sections affected
Bibliography only.

---

## Change 3 (MINOR): Fix author name in `vipond2021`

### Old reference text
```
O.~Vipond, J.~A.~Bull, P.~S.~Sherland, et al.
```

### New corrected reference text
```
O.~Vipond, J.~A.~Bull, P.~S.~Macklin, et al.
```

### Reason
"P. S. Sherland" is not an author of this paper. The third author is Philip S. Macklin. Verified via PNAS and PubMed.

### Sections affected
Bibliography only.

---

## Change 4 (MINOR): Fix year in `chazal2017robust`

### Old reference text
```
\emph{JMLR}, 18(159):1--40, 2017.
```

### New corrected reference text
```
\emph{JMLR}, 18(159):1--40, 2018.
```

### Reason
JMLR volume 18, number 159 was published in 2018, not 2017. Verified via the JMLR website (https://jmlr.org/papers/v18/15-484.html).

### Sections affected
Bibliography only.

---

## Change 5 (MINOR): Improve `chazal2014persistence` formatting

### Old reference text
```
\emph{Springer}, 2016.
```

### New corrected reference text
```
SpringerBriefs in Mathematics. Springer, 2016.
```

### Reason
Added series name for proper book citation formatting. The publication is a book in the SpringerBriefs in Mathematics series, not a journal article.

### Sections affected
Bibliography only.

---

## Change 6 (MINOR): Replace vague `gudhi2014` with formal publication

### Old reference text
```
\bibitem{gudhi2014}
The GUDHI Project.
GUDHI User and Reference Manual, 2014--2024.
```

### New corrected reference text
```
\bibitem{gudhi2014}
C.~Maria, J.-D.~Boissonnat, M.~Glisse, and M.~Yvinec.
The Gudhi library: simplicial complexes and persistent homology.
In \emph{Mathematical Software -- ICMS 2014}, Springer LNCS 8592, pp.~167--174, 2014.
```

### Reason
The original reference was a vague manual citation. The GUDHI library has a formal citable publication from ICMS 2014. Verified via Springer.

### Sections affected
Bibliography only.

---

## Change 7: Body-text wording corrections

### Section 8 "Related Work" — Harmonic PH paragraph

**Old text:**
```
The harmonic PH framework \cite{harmonicPH} introduced the bridge from persistent classes
to sample-level representatives, enabling feature attribution through projection.
CC-SPR builds directly on this work, preserving its scaffold while modifying the
representative objective and generalizing the metric.
```

**New text:**
```
The harmonic PH framework, introduced mathematically by Basu and Cox
\cite{basuCoxHPH2024} and applied to omics data by Gurnari et al.\ \cite{gurnariHPH2025},
provides the bridge from persistent classes to sample-level representatives, enabling
feature attribution through projection.
CC-SPR builds directly on this work, preserving its scaffold while modifying the
representative objective and generalizing the metric.
```

**Reason:** Distinguishes the mathematical foundation from the omics application, as required by the reference split.

---

## Summary of all changes

| # | Type | Reference | What changed |
|---|------|-----------|--------------|
| 1 | MAJOR | harmonicPH | Split into basuCoxHPH2024 + gurnariHPH2025; all in-text citations updated |
| 2 | MINOR | rizvi2017 | Author typo: Kandber → Kandror |
| 3 | MINOR | vipond2021 | Wrong author: Sherland → Macklin |
| 4 | MINOR | chazal2017robust | Year: 2017 → 2018 |
| 5 | MINOR | chazal2014persistence | Added series name |
| 6 | MINOR | gudhi2014 | Replaced vague manual ref with formal ICMS 2014 publication |
| 7 | WORDING | Section 8 body | Distinguished math foundation from omics application |
