# Mass-filter PDF review

2026-10-31. Five-page working note, Alan Day and Chengheng (AI research assistant footnote).

Executed `powershell -ExecutionPolicy Bypass -File code/build_mass_filter.ps1`: two pdfLaTeX passes, exit 0. Desktop compile returned success. Final log has no undefined references or overfull boxes. The final PDF was rendered with `pdftoppm -scale-to 1100 -png`; all five pages were inspected in the preceding work session, with no clipped equations, missing glyphs or table overflow.

Page coverage: (1) four filters and covariance; (2) average mass and normalization; (3) spectral slopes; (4) endpoint and explicitly assumed comparison model; (5) numerical evidence, limitations and references. The PDF reports 23 diagnostic checks; the subsequently completed 444-check raw arithmetic replay is supplementary evidence in audit.json, not an additional claim printed in the PDF.

64 conditional paired scrambled runs and 16 candidate mass-average runs are retained. Errors are empirical; neither positive scalar weights nor these checks establish the missing physical finite-kappa trace. Source, compiler log, reading PDF and this review are bound together by evidence.json. No earlier scientific record was rewritten.
