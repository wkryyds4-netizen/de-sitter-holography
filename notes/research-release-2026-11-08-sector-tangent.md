# Fixed-charge tails and the covariance tangent

This snapshot follows commit 261a31ac806c3cf783226bbb6f05b6c1c123a821. It collects the completed central-sector ordinary-tail stage and a smaller covariance-tangent substage. Both concern the same static Gaussian p-body CAR model, normalized annealed sector Schatten norms, fixed energy normalization and p²/N in a compact subset of (0,infinity).

| Material | Scope | Saved evidence |
| --- | --- | --- |
| [Eight-page gradient-tail PDF](../output/pdf/de-sitter-susskind-sector-gradient-tail.pdf), [LaTeX](../manuscript/susskind-sector-gradient-tail/main.tex) | Adjacent-charge CP gradients, graded derivative shifts and actual ordinary degree-uniform tails in a central O(p) charge window. | [Manifest](../results/susskind-sector-gradient-tail/local-evidence.json) |
| [Four-page covariance-tangent PDF](../output/pdf/de-sitter-susskind-sector-covariance-tangent.pdf), [LaTeX](../manuscript/susskind-sector-covariance-tangent/main.tex) | Exact first-order channel expansion with an ordered two-deficit remainder; fixed-degree local-size tangent with O(p^-1) error. | [Manifest](../results/susskind-sector-covariance-tangent/local-evidence.json) |

The gradient stage records 8,692 exact graded CAR identities and inherits the saved channel checks. The tangent stage records 1,824 exact path identities and 1,824 normalized S2 remainder inequalities on the full neutral matrix-unit basis at N=2,3,4, including vanishing reference factors. Its 60 rational scalar rows check coefficient estimates. The general estimates are analytic arguments using the cited pure-Wick bounds; finite checks are not proofs of the Gaussian asymptotic limit.

The tangent retains finite lambda_N, C_m and theta_(v,m), and its error constant depends on the fixed Wick degree. It controls A_v=p(Phi-C_m theta_v)Z_v. Covariance filtering of the lower-chaos eta descendants, central exact pinching of the full neutral field, and actual memory/bilocal return remain the next interfaces. The complete selected bilocal and the bulk spatial/mass and finite-G dictionary are not claimed complete.

Both PDFs have been compiled twice and visually checked. Saved diagnostics, code, derivations and source notes accompany the evidence. This release verifies existing hashes without rerunning successful scientific drivers. The tangent science remains under results/susskind-sector-neutral-transfer because it is an input to that larger, unfinished stage.

Issue [#1](https://github.com/wkryyds4-netizen/de-sitter-holography/issues/1) is about the earlier KS pentagon branch. This snapshot adds no KS reintegration or diagnosis of the 9,890 endpoint-reproduction discrepancies; the [existing review](issue-1-review.md) retains its stated qualifications. The issue stays open. Research and release preparation used AI assistance. Private context, journals, temporary files and cached third-party full texts are excluded.
