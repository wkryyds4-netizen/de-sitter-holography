# Global Racah note: build and visual review

2026-10-31. Five-page A4 research note, signed Alan Day and Chengheng with the AI research assistant footnote.

- Two pdfLaTeX passes completed through the dedicated build script. Final log has no undefined references, overfull boxes, rerun requests or fatal errors.
- The desktop LaTeX compiler returned success for the final source.
- Poppler rendered all five pages to PNG, maximum dimension 1050 pixels. Every page was visually inspected: authors and footnote, weights and measures, full scalar formula and mapping table, numerical evidence table, proposition and references. No clipping, collisions or missing glyphs were observed.
- The analytic positive-kernel argument is separated from finite-point numerical evidence. The PDF explicitly limits the tested dictionary and does not claim interval certification, physical finite-gravity sewing or a completed holographic dual.
- Numerical validation was already complete before document production; its 44 checks were not rerun merely for layout. Source code and numerical data are bound by hashes in the evidence record.

The first two attempted authoring calls failed at JavaScript parsing because TeX dollar-brace notation was interpreted as template interpolation. A subsequent patch was rejected for two operations on one path. These calls did not modify scientific files. A single update patch succeeded, followed by successful compilation.
