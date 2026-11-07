# 2026-11-08 consolidated Slater-state research snapshot

This snapshot collects four completed stages after `560d683503beb9b6d85832e0bbe6b65f45f54d31`. The latest result is a finite-N orbit-density bound that transfers exact energy-pinching errors to a single selected excitation. The accompanying drafts contain the assumptions and proofs; the finite checks test specified identities and do not independently certify the general theorems.

| Stage | Result and reading copy | Evidence |
| --- | --- | --- |
| Stationary states | [Eight-page draft](../output/pdf/de-sitter-susskind-slater-stationary.pdf): stationary states of the existing h,M norm algebra, absence of a normal stationary density, finite-time averaging, and the covariant fixed-sector trace interface. | [Manifest](../results/susskind-slater-stationary/local-evidence.json) |
| Covariant channel | [Nine-page draft](../output/pdf/de-sitter-susskind-slater-covariant.pdf): the actual same-J size fluctuation u, fixed mixed words and compact-time convergence, and a stationary trace retaining nontrivial h,u dynamics. | [Manifest](../results/susskind-slater-covariant/local-evidence.json) |
| Neutral slope | [Eight-page draft](../output/pdf/de-sitter-susskind-slater-neutral-slope.pdf): the actual rank-to-size bridge, the full q-Fock amplification, the selected-state density, and finite-resolution dephasing. | [Manifest](../results/susskind-slater-neutral-slope/local-evidence.json) |
| Orbit density | [Six-page draft](../output/pdf/de-sitter-susskind-slater-orbit-density.pdf): uniform density moments for homogeneous linear stabilizer-covariant preparations and exact pinching inside a covariant block in the M-dagger prepared state. | [Manifest](../results/susskind-slater-orbit-density/local-evidence.json) |

The microscopic model uses standard CAR, overlapping Hermitian complex Gaussian p-body couplings and the same disorder for the Hamiltonian and selected operators. The asymptotic results use a half-filled, annealed Slater state with p squared over N tending to a fixed positive lambda; the linear orbit-density lemma itself holds at finite N in the stated sector. Energies and times use the normalizations specified in each draft.

The latest estimate is E tr_m(R_P^r) <= (2r-1)!! a_N^r for homogeneous linear preparations, with a_N = E ||P O||^2 and normalized sector trace. For P=M-dagger, a_N=1. Combined with the previously established tracial bounds, it transfers the actual finite-N exact-pinching error into that selected state. Covariant postdressing is included. General nonlinear preparations and words with several interlaced selected insertions are not proved by this estimate.

LaTeX sources, source comparisons, derivation notes, code, scientific outputs and build logs accompany all four PDFs. The PDFs were already compiled twice and visually inspected at their respective milestones. Publication preparation checks saved hashes and Git object bytes without rerunning successful scientific generators or altering sealed manifests. Historical statements that a stage was local describe its status when recorded.

The next calculation concerns quadratic preparations M-dagger h and (M-dagger)^2, their Gaussian-chaos contractions and the partial tensor flattenings needed for dimension-uniform density bounds. A possible extension to fixed polynomial degree remains an unverified working route and is not part of this release. The selected-state interface with J2 and the full bilocal, the bulk spatial and mass generators, and finite-G dynamics remain open. The previously completed full-trace bilocal four-point construction is retained.

Issue [#1](https://github.com/wkryyds4-netizen/de-sitter-holography/issues/1) concerns the earlier KS pentagon branch. This snapshot adds no independent reintegration of the 21,143 cells and does not resolve the 9,890 historical endpoint replay failures. The [existing review](issue-1-review.md) remains the relevant evidence; the new Slater-state results do not strengthen the pentagon claim. The issue remains open.

Research and preparation used AI assistance. Private project context, journals, temporary working material and downloaded third-party attachments remain outside the publication inputs.
