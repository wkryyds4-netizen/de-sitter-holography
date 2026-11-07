# Higher Wick descendants and fixed-charge ordinary fields

This snapshot collects two completed local stages after `b6b2a45b2d0e8e0f7d66edb0aa70c4910e880bd4`. Both use standard CAR, the same Gaussian disorder for all operators, half filling m=N/2, and the double-scaled limit p squared over N tending to fixed positive lambda. They concern fixed polynomial degrees and fixed finite Schatten orders.

| Stage | Reading copy and scope | Evidence |
| --- | --- | --- |
| Higher Wick descendants | [Six-page PDF](../output/pdf/de-sitter-susskind-slater-higher-wick.pdf): general disjoint whole-Wick D/L sources, the selected density functional and fixed-degree exact pinching for the charged tower. | [Manifest](../results/susskind-slater-higher-wick/local-evidence.json) |
| Ordinary charged fields | [Six-page PDF](../output/pdf/de-sitter-susskind-central-ordinary.pdf): exact charge Fourier projection, ordinary internal Gaussian pairs and the same-copy central-sector strong bridge to Wick descendants. | [Manifest](../results/susskind-central-ordinary/local-evidence.json) |

The higher-Wick stage extends the selected density identity to admissible fixed marked words. On the D-only charged tower the functional is weighted by f(h)=1-sqrt(q)+sqrt(q)h squared, with q=exp(-4 lambda). This is a statement about the specified moment algebra; it does not assert that every higher field is bounded or has its tracial marginal.

The ordinary stage retains the wrong-occupancy singleton terms created by internal Gaussian pairs. Exact finite-N Fourier conditioning, with separate near-angle and far-angle estimates, gives K_cen=K_full-2I(z) squared, where I(z) counts the internal D pairs. The resulting fixed-degree same-copy strong bridge transfers to exact energy pinching and the already controlled fixed polynomial preparations. The source kernel is defined in the accompanying derivation; setting a charge variable to zero in a weak limiting law is not used as a substitute for lattice conditioning.

Saved checks include 248 higher-Wick local layouts, 30 density examples and 18 earlier-interface comparisons; the ordinary stage includes 2,186 Fourier coefficients, 96 Gram equalities and a separate 5,220-case local-source check. Its first checker accidentally enumerated zero local layouts. That coverage gap and its original output are preserved; only the separate nonempty supplementary enumeration supplies the local coverage. High-precision Fourier quadrature checks are consistency tests with a stated cutoff-tail bound, not interval quadrature certificates.

Both PDFs were compiled twice and every page inspected at their milestones. This publication checks saved file hashes and PDF page counts without rerunning successful scientific drivers or changing sealed evidence. LaTeX, derivations, source notes, code, outputs and diagnostic records accompany the PDFs. Historical statements that these stages were local retain their time-of-record meaning.

The next task is a degree-uniform finite-N ordinary tail bound in the central charge sector. A candidate route uses row and column CAR gradients along normalized charge-changing CP channels; its identities and tagged-gradient estimates still require verification. Full-Fock Clifford Riesz estimates cannot simply be restricted to a fixed charge sector. The complete selected bilocal also needs its neutral and memory-return interfaces checked in that sector. Spatial/mass generators and finite-G dynamics remain later tasks.

Issue [#1](https://github.com/wkryyds4-netizen/de-sitter-holography/issues/1) concerns the earlier KS pentagon branch. There is no new KS reintegration or endpoint-discrepancy result in this snapshot. The [existing review](issue-1-review.md) and its qualifications remain applicable, and the issue remains open.

Research and preparation used AI assistance. Private context, journals, temporary files and third-party full texts remain outside this publication.
