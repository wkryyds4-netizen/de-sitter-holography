"""Check research records and current result hashes; freeze a unique verified ZIP."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import re
import zipfile
from decimal import Decimal, localcontext

ROOT = Path(__file__).resolve().parent.parent

RECORDS = [
    ("results/study.json", "stage1-study", False),
    ("results/validation.json", "stage1", True),
    ("results/stage2/study.json", "stage2-study", False),
    ("results/stage2/validation.json", "stage2", True),
    ("results/stage2/operator-bridge-validation.json", "operator-bridge", True),
    ("results/stage3/study.json", "stage3-recoil-study", False),
    ("results/stage3/validation.json", "stage3-recoil", True),
    ("results/stage3/packets/study.json", "stage3-packets-study", False),
    ("results/stage3/packets/validation.json", "stage3-packets", True),
    ("results/stage3/fourpoint/uncrossed-study.json", "stage3-fourpoint-study", False),
    ("results/stage3/fourpoint/uncrossed-validation.json", "stage3-fourpoint", True),
    ("results/stage3/fourpoint/clock-limit-study.json", "stage3-clock-limit-study", False),
    ("results/stage3/fourpoint/clock-limit-validation.json", "stage3-clock-limit", True),
    ("results/stage3/fourpoint/crossing-validation.json", "stage3-crossing", True),
    ("results/stage3/clock-rate/study.json", "stage3-clock-rate-study", False),
    ("results/stage3/clock-rate/validation.json", "stage3-clock-rate", True),
    ("results/stage3/fourpoint/crossing-regular-validation.json", "stage3-crossing-regular", True),
    ("results/stage3/fourpoint/conditional-crossing-pilot.json", "conditional-crossing-pilot", False),
    ("results/stage3/fourpoint/conditional-crossing-refined.json", "conditional-crossing-refined", False),
    ("results/stage3/fourpoint/conditional-crossing-cartesian.json", "conditional-crossing-cartesian", False),
    ("results/stage3/fourpoint/conditional-validation.json", "conditional-crossing", True),
    ("results/stage3/packet-crossing/setup-validation.json", "packet-crossing-setup", True),
    ("results/stage3/packet-crossing/positive-validation.json", "packet-crossing-positive", True),
    ("results/stage3/packet-crossing/tail-validation.json", "packet-crossing-tail", True),
    ("results/stage3/packet-crossing/packet-crossing-study.json", "packet-crossing-study", False),
    ("results/stage3/packet-crossing/validation.json", "packet-crossing", True),
    ("results/stage4/study.json", "stage4-study", False),
    ("results/stage4/validation.json", "stage4", True),
    ("results/stage4/microscopic-comparison.json", "microscopic-comparison", True),
    ("results/stage4/time-corrections.json", "microscopic-time-corrections", True),
]

# Historical failures stay available and are included in the byte inventory.
# Their code hashes deliberately refer to the earlier implementations.
HISTORICAL_RECORDS = [
    "results/stage2/validation-initial-grid.json",
    "results/stage2/tails-initial-grid.json",
    "results/stage2/study-initial-grid.json",
    "results/stage2/operator-bridge-initial-diagnostic.json",
    "results/stage3/fourpoint/clock-limit-initial-unsplit-study.json",
    "results/stage3/fourpoint/clock-limit-initial-unsplit-validation.json",
    "results/stage3/fourpoint/crossing-regular-initial-validation.json",
]
CROSSING_PROVENANCE = {
    "code/observer_crossing.py", "code/validate_observer_crossing.py",
    "sources/kolchmeyer-liu-2411.08090v1.html",
    "literature/stage3-fourpoint-source-map.md",
}
COMPLETED_STUDIES = {
    "conditional-crossing-pilot", "conditional-crossing-refined", "conditional-crossing-cartesian",
    "packet-crossing-study",
}
ROOT_HASH_RECORDS = {
    "stage3-clock-rate-study", "stage3-clock-rate",
}
PACKET_RAW_ARCHIVE = "results/stage3/packet-crossing/packet-crossing-study-points.json.gz"
PACKET_VALIDATIONS = {
    "packet-crossing-setup", "packet-crossing-positive", "packet-crossing-tail", "packet-crossing",
}
EXPECTED_VALIDATION_COUNTS = {
    "packet-crossing-setup": 28, "packet-crossing-positive": 34, "packet-crossing-tail": 29,
    "packet-crossing": 93,
}
EXPECTED_INPUT_HASHES = {
    "packet-crossing-setup": {"sources/kolchmeyer-liu-2411.08090v1.html"},
    "packet-crossing-study": {
        "sources/kolchmeyer-liu-2411.08090v1.html", PACKET_RAW_ARCHIVE,
    },
    "packet-crossing-tail": {
        "results/stage3/packet-crossing/positive-validation.json",
        "notes/stage3-packet-tail-envelope.md",
    },
    "packet-crossing": {
        "results/stage3/packet-crossing/packet-crossing-study.json",
        "results/stage3/packet-crossing/setup-validation.json",
        "results/stage3/packet-crossing/positive-validation.json",
        "results/stage3/packet-crossing/tail-validation.json",
        "results/stage3/fourpoint/conditional-crossing-cartesian.json",
        PACKET_RAW_ARCHIVE,
    },
    "stage3-crossing-regular": {
        "sources/kolchmeyer-liu-2411.08090v1.html",
        "results/stage3/fourpoint/crossing-validation.json",
    },
    "conditional-crossing": {
        "results/stage3/fourpoint/conditional-crossing-pilot.json",
        "results/stage3/fourpoint/conditional-crossing-refined.json",
        "results/stage3/fourpoint/conditional-crossing-cartesian.json",
    },
}
EXPECTED_RATE_HASHES = {
    "code/validate_clock_rate.py", "code/quantized_observer.py",
    "code/semiclassical_clock.py", "notes/stage3-clock-rate.md",
}
# Require every archived source to remain represented in the checksum manifest.
# New future sources are allowed; dropping one of these would otherwise go unseen.
REQUIRED_SOURCE_ARCHIVES = {
    "sources/giulini-marolf-gr-qc9902045v2.html",
    "sources/dlmf-14.3.html",
    "sources/dlmf-5.11.html",
    "sources/han-hu-2610.14175v1.html",
    "sources/cui-kolchmeyer-2607.13665v1.html",
    "sources/harlow-zhao-2607.14215v2.html",
    "sources/collier-et-al-2501.01486v1.html",
    "sources/ali-et-al-2003.11046v2.html",
    "sources/witten-2303.02837v2.html",
    "sources/chen-stanford-tang-yang-2607.14042v1.html",
    "sources/narovlansky-2506.02109v2.html",
    "sources/pinamonti-0806.0803v2.pdf",
    "sources/neeb-2011-aif2660.pdf",
    "sources/kolchmeyer-liu-2411.08090v1.html",
    "sources/chandra-liska-post-2610.30247v1.html",
    "sources/bros-et-al-0901.4223v1.html",
    "sources/bros-et-al-0901.4223v1.pdf",
    "sources/dlmf-14.24.html",
    "sources/marolf-yu-2610.24073v1.html",
    "sources/goldberg-schlag-math0306108v2.html",
    "sources/germain-pusateri-rousset-2018.pdf",
    "sources/dlmf-15.8.html",
    "sources/dlmf-15.2.html",
    "sources/sun-2111.04591v2.html",
    "sources/dlmf-14.12.html",
    "sources/dlmf-14.20.html",
}
# These are evidence/data records, not passed-check summaries.
CG_AUDIT_DIR = "results/stage3/cg-domain-source-audit/"
CG_ROOT_HASH_RECORDS = {
    CG_AUDIT_DIR + "evidence.json": {
        "notes/stage3-cg-domain-investigation.md",
        CG_AUDIT_DIR + "kolchmeyer-liu-2411.08090v1.pdf",
        CG_AUDIT_DIR + "v1-page76.png",
    },
    CG_AUDIT_DIR + "vertex-analytic-reference.json": {
        "code/validate_vertex_gamma_reference.py", "code/validate_vertex_source.py",
        CG_AUDIT_DIR + "vertex-check.json",
    },
}
ORBIT_IDENTITY_RECORD = "results/stage3/orbit-core/orbit-core-identity-70dps.json"
CG_ROOT_HASH_RECORDS[ORBIT_IDENTITY_RECORD] = {
    "code/validate_orbit_core_identities.py",
    "sources/kolchmeyer-liu-2411.08090v1.html",
    "sources/dlmf-15.8.html", "sources/germain-pusateri-rousset-2018.pdf",
}
MANUSCRIPT_FILES = [
    "manuscript/" + name for name in (
        "main.tex", "orbit.tex", "spectral.tex", "physical.tex", "fields.tex",
        "fourpoint.tex", "normalization.tex", "vertex.tex", "references.tex", "references.bib",
        "position.tex", "annihilation.tex", "crosslimit.tex", "endpoint.tex", "assembly.tex", "logtime.tex", "preparation.tex", "movingcurve.tex", "densityendpoints.tex", "twobranch.tex", "growingbox.tex"
    )
]
CG_ROOT_HASH_RECORDS["results/manuscript/build-record.json"] = set(MANUSCRIPT_FILES) | {
    "code/build_manuscript.ps1", "output/pdf/de-sitter-observer-foundations.pdf",
    "results/manuscript/main.log", "results/manuscript/main.blg",
    "results/stage3/preparation/fixed-curve.pdf",
    "results/stage3/moving-curve/moving-curve.pdf",
    "results/stage3/density-endpoints/endpoint-curves.pdf",
    "results/stage3/two-branch/branch-window.pdf",
}
CG_ROOT_HASH_RECORDS["results/stage3/position-action/evidence.json"] = {
    "manuscript/position.tex", "notes/stage3-direct-position-creation.md",
    "literature/stage3-position-action-sources.md",
    "notes/stage3-positive-frequency-orbit-core.md",
    "notes/stage3-ambient-spectral-test-core.md",
    "notes/stage3-maximal-weak-frequency-space.md",
    "notes/stage3-all-n-generated-space.md",
    "notes/stage3-regulated-local-descent.md",
    "sources/kolchmeyer-liu-2411.08090v1.html",
    "sources/giulini-marolf-gr-qc9902045v2.html",
}

CG_ROOT_HASH_RECORDS["results/stage3/annihilation/evidence.json"] = {
    "manuscript/annihilation.tex", "notes/stage3-integrated-annihilation.md",
    "literature/stage3-annihilation-sources.md", "manuscript/position.tex",
    "manuscript/physical.tex", "manuscript/fields.tex",
    "notes/stage3-regulated-local-descent.md",
    "notes/stage3-product-seed-ga-density.md",
    "sources/kolchmeyer-liu-2411.08090v1.html",
}

CG_ROOT_HASH_RECORDS["results/stage3/cross-limit/comparison.json"] = {
    "code/cross_limit_comparison.py", "code/semiclassical_clock.py",
    "code/quantized_observer.py", "manuscript/crosslimit.tex",
    "sources/kolchmeyer-liu-2411.08090v1.html",
}
CG_ROOT_HASH_RECORDS["results/stage3/cross-limit/evidence.json"] = {
    "manuscript/crosslimit.tex", "notes/stage3-crossed-limit-reduction.md",
    "literature/stage3-crossed-large-mass-sources.md",
    "notes/stage3-crossing-integral-audit.md", "notes/stage3-clock-rate.md",
    "manuscript/fields.tex", "manuscript/fourpoint.tex",
    "results/stage3/cross-limit/comparison.json",
    "sources/kolchmeyer-liu-2411.08090v1.html",
}

CG_ROOT_HASH_RECORDS["results/stage3/endpoint/diagnostics.json"] = {
    "code/hyper_endpoint.py", "code/observer_crossing.py", "manuscript/endpoint.tex",
    "sources/kolchmeyer-liu-2411.08090v1.html", "sources/dlmf-15.8.html",
}
CG_ROOT_HASH_RECORDS["results/stage3/endpoint/evidence.json"] = {
    "manuscript/endpoint.tex", "notes/stage3-hypergeometric-endpoint.md",
    "literature/stage3-endpoint-sources.md", "results/stage3/endpoint/diagnostics.json",
    "manuscript/crosslimit.tex", "code/hyper_endpoint.py",
    "sources/kolchmeyer-liu-2411.08090v1.html", "sources/dlmf-15.8.html",
}
CG_ROOT_HASH_RECORDS["results/stage3/assembly/diagnostics.json"] = {
    "code/crossing_assembly.py", "code/observer_crossing.py", "code/observer_crossing_regular.py",
    "manuscript/assembly.tex", "manuscript/endpoint.tex", "sources/kolchmeyer-liu-2411.08090v1.html",
}
CG_ROOT_HASH_RECORDS["results/stage3/assembly/evidence.json"] = {
    "manuscript/assembly.tex", "manuscript/endpoint.tex", "manuscript/crosslimit.tex",
    "notes/stage3-assembled-fixed-time.md", "literature/stage3-assembly-sources.md",
    "notes/stage3-crossing-integral-audit.md", "notes/stage3-clock-rate.md",
    "results/stage3/assembly/diagnostics.json", "sources/kolchmeyer-liu-2411.08090v1.html",
}
CG_ROOT_HASH_RECORDS["results/stage3/logtime/diagnostics.json"] = {
    "code/logtime_limit.py", "code/crossing_assembly.py", "manuscript/logtime.tex",
    "manuscript/assembly.tex", "manuscript/crosslimit.tex", "sources/kolchmeyer-liu-2411.08090v1.html",
}
CG_ROOT_HASH_RECORDS["results/stage3/logtime/evidence.json"] = {
    "manuscript/logtime.tex", "manuscript/assembly.tex", "manuscript/crosslimit.tex",
    "manuscript/fourpoint.tex", "notes/stage3-logarithmic-time.md", "literature/stage3-logtime-sources.md",
    "results/stage3/logtime/diagnostics.json", "sources/kolchmeyer-liu-2411.08090v1.html",
}
CG_ROOT_HASH_RECORDS["results/stage3/preparation/diagnostics.json"] = {
    "code/preparation_asymptotics.py", "code/logtime_limit.py", "code/crossing_assembly.py",
    "manuscript/preparation.tex", "manuscript/logtime.tex",
    "sources/kolchmeyer-liu-2411.08090v1.html", "sources/dlmf-5.11.html",
}
CG_ROOT_HASH_RECORDS["results/stage3/preparation/evidence.json"] = {
    "manuscript/preparation.tex", "manuscript/logtime.tex", "manuscript/assembly.tex",
    "notes/stage3-preparation-asymptotics.md", "literature/stage3-preparation-sources.md",
    "results/stage3/preparation/diagnostics.json", "results/stage3/preparation/fixed-curve.pdf",
    "results/stage3/preparation/fixed-curve.png", "sources/kolchmeyer-liu-2411.08090v1.html",
    "sources/dlmf-5.11.html",
}
CG_SOURCE_ARCHIVES = {
    "sources/bros-et-al-0901.4223v1.html", "sources/bros-et-al-0901.4223v1.pdf",
    "sources/dlmf-14.24.html",
}
MOVING_RUNS = [f"results/stage3/moving-curve/{name}" for name in (
    "run-w24-chi40-R10.json", "run-w32-chi40-R10.json", "run-w24-chi64-R10.json")]
for moving_run in MOVING_RUNS:
    CG_ROOT_HASH_RECORDS[moving_run] = {
        "code/run_moving_curve.py", "code/moving_threshold_curve.py", "code/preparation_asymptotics.py",
        "manuscript/preparation.tex", "manuscript/logtime.tex", "sources/kolchmeyer-liu-2411.08090v1.html",
    }
CG_ROOT_HASH_RECORDS["results/stage3/moving-curve/diagnostics.json"] = set(MOVING_RUNS) | {
    "code/validate_moving_curve.py", "code/run_moving_curve.py", "code/moving_threshold_curve.py",
    "code/preparation_asymptotics.py",
}
CG_ROOT_HASH_RECORDS["results/stage3/moving-curve/evidence.json"] = set(MOVING_RUNS) | {
    "manuscript/movingcurve.tex", "notes/stage3-moving-threshold-curve.md",
    "literature/stage3-moving-curve-sources.md", "results/stage3/moving-curve/diagnostics.json",
    "results/stage3/moving-curve/moving-curve.pdf", "results/stage3/moving-curve/moving-curve.png",
    "manuscript/preparation.tex", "manuscript/logtime.tex", "manuscript/crosslimit.tex",
}
CG_SOURCE_RECORD = CG_AUDIT_DIR + "analytic-sources.json"
DENSITY_RUNS = [f"results/stage3/density-endpoints/{name}" for name in (
    "run-w24-chi40.json", "run-w32-chi40.json", "run-w24-chi64.json")]
for density_run in DENSITY_RUNS:
    CG_ROOT_HASH_RECORDS[density_run] = {
        "code/run_density_curve.py", "code/density_endpoints.py", "code/moving_threshold_curve.py",
        "code/preparation_asymptotics.py", "manuscript/preparation.tex", "manuscript/movingcurve.tex",
    }
CG_ROOT_HASH_RECORDS["results/stage3/density-endpoints/diagnostics.json"] = set(DENSITY_RUNS) | {
    "code/validate_density_endpoints.py", "code/run_density_curve.py", "code/density_endpoints.py",
    "code/moving_threshold_curve.py", "code/preparation_asymptotics.py",
}
CG_ROOT_HASH_RECORDS["results/stage3/density-endpoints/evidence.json"] = set(DENSITY_RUNS) | {
    "manuscript/densityendpoints.tex", "notes/stage3-density-endpoints.md",
    "literature/stage3-density-endpoint-sources.md", "results/stage3/density-endpoints/diagnostics.json",
    "results/stage3/density-endpoints/endpoint-curves.pdf", "results/stage3/density-endpoints/endpoint-curves.png",
    "manuscript/preparation.tex", "manuscript/movingcurve.tex", "manuscript/logtime.tex",
}
CG_VERTEX_DATA = [CG_AUDIT_DIR + "vertex-check.json", CG_AUDIT_DIR + "vertex-check-pilot.json"]
CG_ROOT_HASH_RECORDS["results/stage3/two-branch/diagnostics.json"] = set(DENSITY_RUNS + MOVING_RUNS) | {
    "code/two_branch_window.py", "code/logtime_limit.py", "code/crossing_assembly.py",
    "code/moving_threshold_curve.py", "code/preparation_asymptotics.py", "manuscript/twobranch.tex",
    "sources/kolchmeyer-liu-2411.08090v1.html",
}
CG_ROOT_HASH_RECORDS["results/stage3/two-branch/evidence.json"] = {
    "manuscript/twobranch.tex", "notes/stage3-two-branch-window.md", "literature/stage3-two-branch-sources.md",
    "manuscript/assembly.tex", "manuscript/logtime.tex", "manuscript/crosslimit.tex",
    "manuscript/preparation.tex", "manuscript/densityendpoints.tex",
    "results/stage3/two-branch/diagnostics.json", "results/stage3/two-branch/branch-window.pdf",
    "results/stage3/two-branch/branch-window.png", "code/two_branch_window.py",
    "sources/kolchmeyer-liu-2411.08090v1.html",
}
REQUIRED = [
    "notes/stage3-growing-box-windows.md", "literature/stage3-growing-box-sources.md", "code/growing_box_bounds.py",
    "notes/stage3-two-branch-window.md", "literature/stage3-two-branch-sources.md",
    "code/two_branch_window.py", "results/stage3/two-branch/branch-window.pdf", "results/stage3/two-branch/branch-window.png",
    "notes/stage3-density-endpoints.md", "literature/stage3-density-endpoint-sources.md",
    "code/density_endpoints.py", "code/run_density_curve.py", "code/validate_density_endpoints.py",
    "results/stage3/density-endpoints/endpoint-curves.pdf", "results/stage3/density-endpoints/endpoint-curves.png",
    "notes/stage3-moving-threshold-curve.md", "literature/stage3-moving-curve-sources.md",
    "code/moving_threshold_curve.py", "code/run_moving_curve.py", "code/validate_moving_curve.py",
    "results/stage3/moving-curve/moving-curve.pdf", "results/stage3/moving-curve/moving-curve.png",
    "notes/stage3-preparation-asymptotics.md", "literature/stage3-preparation-sources.md",
    "code/preparation_asymptotics.py", "results/stage3/preparation/fixed-curve.pdf",
    "results/stage3/preparation/fixed-curve.png",
    "notes/stage3-logarithmic-time.md", "literature/stage3-logtime-sources.md", "code/logtime_limit.py",
    "notes/stage3-assembled-fixed-time.md", "literature/stage3-assembly-sources.md",
    "code/crossing_assembly.py",
    "notes/stage3-hypergeometric-endpoint.md", "literature/stage3-endpoint-sources.md",
    "code/hyper_endpoint.py", "journal/2026-10-30.md",
    "notes/stage3-crossed-limit-reduction.md", "literature/stage3-crossed-large-mass-sources.md",
    "code/cross_limit_comparison.py",
    "notes/stage3-integrated-annihilation.md", "literature/stage3-annihilation-sources.md",
    "notes/stage3-direct-position-creation.md", "literature/stage3-position-action-sources.md",
    *MANUSCRIPT_FILES,
    "code/build_manuscript.ps1", "output/pdf/de-sitter-observer-foundations.pdf",
    "results/manuscript/main.log", "results/manuscript/main.blg", "results/manuscript/build-review.md",
    "notes/stage3-ambient-spectral-test-core.md", "notes/stage3-maximal-weak-frequency-space.md",
    "README.md", "PROJECT.md", "PLAN.md", "NEXT.md", "AGENTS.md", "requirements.txt",
    "literature/observers.md", "literature/models.md", "literature/novelty-and-next.md",
    "notes/stage1.md", "notes/math-independent.md", "notes/physics-audit.md",
    "journal/2026-10-29.md", "code/clock_model.py", "code/run_study.py",
    "code/validate_clock.py", "results/clock-costs.png", "results/compact-clock-timing.png",
    "sources/manifest.json", "literature/stage2-field-spectrum.md",
    "notes/stage2-wick-measures.md", "notes/stage2-physics-audit.md", "notes/stage2.md",
    "notes/stage2-localization-bridge.md", "notes/stage2-weyl-bridge.md",
    "code/field_model.py", "code/run_field_study.py", "code/validate_field.py",
    "code/validate_operator_bridge.py", "results/stage2/two-vs-four.png",
    "results/stage2/shared-clock-otoc.png", "literature/stage3-recoil-interface.md",
    "notes/stage3-baseline.md", "code/recoil_model.py",
    "code/run_recoil_study.py", "code/validate_recoil.py",
    "literature/stage3-quantized-observer-dictionary.md",
    "literature/stage3-observer-kernel-audit.md", "notes/stage3-packets.md",
    "code/quantized_observer.py", "code/run_quantized_observer.py",
    "notes/stage3-fourpoint-domains.md", "notes/stage3-semiclassical-clock-bridge.md",
    "literature/stage3-fourpoint-source-map.md", "notes/stage3-fourpoint-study.md",
    "code/observer_fourpoint.py", "code/run_observer_fourpoint.py",
    "code/semiclassical_clock.py", "code/run_semiclassical_clock.py",
    "code/observer_crossing.py", "code/validate_observer_crossing.py",
    "notes/stage3-clock-rate.md", "code/validate_clock_rate.py",
    "notes/stage3-closed-domain-next.md", "notes/stage3-clock-first-correction.md",
    "notes/stage3-crossing-regularization.md", "notes/stage3-crossing-integral-audit.md",
    "code/observer_crossing_regular.py", "code/validate_crossing_regular.py",
    "code/run_crossing_conditional.py", "code/run_crossing_refined.py",
    "code/run_crossing_cartesian.py", "code/validate_crossing_conditional.py",
    "notes/stage3-crossing-conditional-study.md", "literature/october-tqft-update.md",
    "notes/stage4-conformal-setup.md", "literature/stage4-microscopic-match.md",
    "code/conformal_shock.py", "code/run_conformal_shock.py",
    "code/validate_microscopic_comparison.py", "results/stage4/conformal-shock.png",
    "literature/stage4-time-corrections.md", "code/microscopic_time_corrections.py",
    "code/crossing_point_cache.py", "code/packet_quadrature.py", "code/run_crossing_packet.py",
    "code/validate_packet_setup.py", "code/packet_positive.py", "code/validate_packet_positive.py",
    "code/validate_packet_tail.py", "code/validate_crossing_packet.py",
    "notes/stage3-packet-quadrature-design.md", "notes/stage3-packet-tail-envelope.md",
    "literature/tqft-recoil-dictionary.md", "notes/stage3-cg-domain-investigation.md",
    "notes/stage3-vertex-source-check.md", "notes/stage3-vertex-analytic-audit.md",
    "code/validate_vertex_source.py", "code/validate_vertex_gamma_reference.py",
    "results/stage3/packet-crossing/point-cache.sqlite", PACKET_RAW_ARCHIVE,
    CG_AUDIT_DIR + "kolchmeyer-liu-2411.08090v1.pdf", CG_AUDIT_DIR + "v1-page76.png",
    CG_SOURCE_RECORD,
    "notes/stage3-packet-crossing-study.md", "notes/stage3-fp-multiplicity-bridge.md",
    CG_AUDIT_DIR + "v1-page26.png",
    "literature/stage3-two-particle-cg-map.md",
    "literature/stage3-group-averaging-update.md",
    "notes/stage3-ga-spectral-bridge.md",
    "notes/stage3-fp-two-particle-inner-product.md",
    "notes/stage3-creation-descent-audit.md",
    "notes/stage3-fp-orbit-core-candidate.md",
    "notes/stage3-radial-scattering-core.md",
    "notes/stage3-positive-frequency-orbit-core.md",
    "notes/stage3-product-seed-ga-density.md",
    "notes/stage3-finite-k-ga-extension.md",
    "notes/stage3-all-n-generated-space.md",
    "notes/stage3-regulated-local-descent.md",
    "code/validate_orbit_core_identities.py",
] + [rel for rel, _, _ in RECORDS] + HISTORICAL_RECORDS + sorted(REQUIRED_SOURCE_ARCHIVES) + list(CG_ROOT_HASH_RECORDS) + CG_VERTEX_DATA

CG_ROOT_HASH_RECORDS["results/stage3/growing-box/diagnostics.json"] = {
    "code/growing_box_bounds.py", "sources/dlmf-5.11.html", "sources/dlmf-15.8.html",
    "sources/kolchmeyer-liu-2411.08090v1.html", "manuscript/endpoint.tex", "manuscript/assembly.tex", "manuscript/twobranch.tex",
}
CG_ROOT_HASH_RECORDS["results/stage3/growing-box/evidence.json"] = {
    "manuscript/growingbox.tex", "notes/stage3-growing-box-windows.md", "literature/stage3-growing-box-sources.md",
    "manuscript/endpoint.tex", "manuscript/assembly.tex", "manuscript/twobranch.tex", "manuscript/crosslimit.tex",
    "manuscript/preparation.tex", "manuscript/densityendpoints.tex", "code/growing_box_bounds.py",
    "results/stage3/growing-box/diagnostics.json", "sources/dlmf-5.11.html", "sources/dlmf-15.8.html",
    "sources/kolchmeyer-liu-2411.08090v1.html",
}
REQUIRED.extend(["results/stage3/growing-box/diagnostics.json", "results/stage3/growing-box/evidence.json"])

CG_ROOT_HASH_RECORDS["results/stage4/wkb-overlap/diagnostics.json"] = {
    "code/wkb_overlap.py", "sources/chen-stanford-tang-yang-2607.14042v1.html",
    "sources/chandra-liska-post-2610.30247v1.html",
}
CG_ROOT_HASH_RECORDS["results/stage4/wkb-overlap/evidence.json"] = {
    "code/wkb_overlap.py", "code/build_wkb_note.ps1", "manuscript/finite-g-bridge/main.tex",
    "notes/stage4-wkb-overlap.md", "literature/stage4-wkb-overlap-sources.md",
    "results/stage4/wkb-overlap/diagnostics.json", "output/pdf/de-sitter-wkb-overlap.pdf",
    "results/stage4/wkb-overlap/main.log", "sources/chen-stanford-tang-yang-2607.14042v1.html",
    "sources/chandra-liska-post-2610.30247v1.html",
}
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS["results/stage4/wkb-overlap/evidence.json"]))
REQUIRED.extend(["results/stage4/wkb-overlap/evidence.json"])

CG_ROOT_HASH_RECORDS["results/stage4/gram-prefactor/diagnostics.json"] = {
    "code/gram_prefactor.py", "code/wkb_overlap.py",
    "sources/chen-stanford-tang-yang-2607.14042v1.html",
    "sources/chandra-liska-post-2610.30247v1.html",
}
CG_ROOT_HASH_RECORDS["results/stage4/gram-prefactor/evidence.json"] = {
    "code/gram_prefactor.py", "code/wkb_overlap.py", "code/build_gram_note.ps1",
    "manuscript/gram-prefactor/main.tex", "notes/stage4-gram-prefactor.md",
    "literature/stage4-gram-prefactor-sources.md",
    "results/stage4/gram-prefactor/diagnostics.json",
    "results/stage4/gram-prefactor/main.log",
    "results/stage4/gram-prefactor/build-review.md",
    "output/pdf/de-sitter-gram-prefactor.pdf",
    "sources/chen-stanford-tang-yang-2607.14042v1.html",
    "sources/chandra-liska-post-2610.30247v1.html",
}
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS["results/stage4/gram-prefactor/evidence.json"]))
REQUIRED.append("results/stage4/gram-prefactor/evidence.json")

CG_ROOT_HASH_RECORDS["results/stage4/recoil-vertex/diagnostics.json"] = {
    "code/recoil_vertex.py", "sources/chandra-liska-post-2610.30247v1.html",
    "sources/dlmf-5.11.html",
}
CG_ROOT_HASH_RECORDS["results/stage4/recoil-vertex/evidence.json"] = {
    "code/recoil_vertex.py", "code/build_vertex_note.ps1",
    "manuscript/recoil-vertex/main.tex", "notes/stage4-recoil-vertex.md",
    "literature/stage4-recoil-vertex-sources.md",
    "results/stage4/recoil-vertex/diagnostics.json", "results/stage4/recoil-vertex/main.log",
    "results/stage4/recoil-vertex/build-review.md", "output/pdf/de-sitter-recoil-vertex.pdf",
    "sources/chandra-liska-post-2610.30247v1.html", "sources/dlmf-5.11.html",
}
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS["results/stage4/recoil-vertex/evidence.json"]))
REQUIRED.append("results/stage4/recoil-vertex/evidence.json")

CG_ROOT_HASH_RECORDS["results/stage4/recoil-internal/diagnostics.json"] = {
    "code/recoil_internal.py", "sources/chandra-liska-post-2610.30247v1.html",
    "sources/dlmf-5.11.html",
}
CG_ROOT_HASH_RECORDS["results/stage4/recoil-internal/evidence.json"] = {
    "code/recoil_internal.py", "code/build_internal_note.ps1",
    "manuscript/recoil-internal/main.tex", "notes/stage4-recoil-internal.md",
    "literature/stage4-recoil-internal-sources.md",
    "results/stage4/recoil-internal/diagnostics.json", "results/stage4/recoil-internal/main.log",
    "results/stage4/recoil-internal/build-review.md",
    "results/stage4/recoil-internal/initial-build-failure.txt",
    "results/stage4/recoil-internal/missing-package-build-failure.txt",
    "output/pdf/de-sitter-recoil-internal.pdf", "manuscript/recoil-vertex/main.tex",
    "sources/chandra-liska-post-2610.30247v1.html", "sources/dlmf-5.11.html",
}
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS["results/stage4/recoil-internal/evidence.json"]))
REQUIRED.append("results/stage4/recoil-internal/evidence.json")

CG_ROOT_HASH_RECORDS["results/stage4/recoil-saddle/diagnostics.json"] = {
    "code/recoil_saddle.py", "code/recoil_internal.py", "code/gram_prefactor.py",
    "sources/chen-stanford-tang-yang-2607.14042v1.html",
    "sources/chandra-liska-post-2610.30247v1.html",
}
CG_ROOT_HASH_RECORDS["results/stage4/recoil-saddle/evidence.json"] = {
    "code/recoil_saddle.py", "code/build_saddle_note.ps1",
    "manuscript/recoil-saddle/main.tex", "notes/stage4-recoil-saddle.md",
    "literature/stage4-recoil-saddle-sources.md",
    "results/stage4/recoil-saddle/diagnostics.json", "results/stage4/recoil-saddle/main.log",
    "results/stage4/recoil-saddle/build-review.md",
    "results/stage4/recoil-saddle/initial-branch-diagnostics.json",
    "results/stage4/recoil-saddle/initial-overfull.log",
    "output/pdf/de-sitter-recoil-saddle.pdf", "code/recoil_internal.py",
    "code/gram_prefactor.py", "manuscript/recoil-vertex/main.tex",
    "sources/chen-stanford-tang-yang-2607.14042v1.html",
    "sources/chandra-liska-post-2610.30247v1.html",
}
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS["results/stage4/recoil-saddle/evidence.json"]))
REQUIRED.append("results/stage4/recoil-saddle/evidence.json")

CG_ROOT_HASH_RECORDS["results/stage4/recoil-local/diagnostics.json"] = {
    "code/recoil_local.py", "code/recoil_saddle.py", "code/recoil_internal.py", "code/gram_prefactor.py",
    "sources/chen-stanford-tang-yang-2607.14042v1.html", "sources/dlmf-5.11.html",
}
CG_ROOT_HASH_RECORDS["results/stage4/recoil-local/evidence.json"] = {
    "code/recoil_local.py", "code/build_local_note.ps1", "manuscript/recoil-local/main.tex",
    "notes/stage4-recoil-local.md", "literature/stage4-recoil-local-sources.md",
    "results/stage4/recoil-local/diagnostics.json", "results/stage4/recoil-local/main.log",
    "results/stage4/recoil-local/build-review.md", "results/stage4/recoil-local/initial-norm-grid.json",
    "results/stage4/recoil-local/initial-code.py.txt", "output/pdf/de-sitter-recoil-local.pdf",
    "code/recoil_saddle.py", "code/recoil_internal.py", "code/gram_prefactor.py",
    "manuscript/recoil-saddle/main.tex", "sources/chen-stanford-tang-yang-2607.14042v1.html",
    "sources/dlmf-5.11.html",
}
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS["results/stage4/recoil-local/evidence.json"]))
REQUIRED.append("results/stage4/recoil-local/evidence.json")


CG_ROOT_HASH_RECORDS["results/stage4/recoil-global/diagnostics.json"] = {
    "code/recoil_global_probe.py", "code/recoil_saddle.py", "code/recoil_local.py",
    "code/recoil_internal.py", "sources/dlmf-5.11.html",
}
CG_ROOT_HASH_RECORDS["results/stage4/recoil-global/evidence.json"] = {
    "code/recoil_global_probe.py", "code/build_global_note.ps1",
    "manuscript/recoil-global/main.tex", "notes/stage4-recoil-global.md",
    "literature/stage4-recoil-global-sources.md",
    "results/stage4/recoil-global/diagnostics.json",
    "results/stage4/recoil-global/main.log", "results/stage4/recoil-global/build-review.md",
    "output/pdf/de-sitter-recoil-global.pdf",
    "code/recoil_saddle.py", "code/recoil_local.py", "code/recoil_internal.py",
    "sources/dlmf-5.11.html",
}
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS["results/stage4/recoil-global/evidence.json"]))
REQUIRED.append("results/stage4/recoil-global/evidence.json")


CG_ROOT_HASH_RECORDS["results/stage4/recoil-bessel-limit/diagnostics.json"] = {
    "code/recoil_bessel_limit.py", "code/recoil_global_probe.py",
    "results/stage4/recoil-global/diagnostics.json", "sources/dlmf-5.11.html",
}
CG_ROOT_HASH_RECORDS["results/stage4/recoil-bessel-limit/evidence.json"] = {
    "code/recoil_bessel_limit.py", "code/build_bessel_limit_note.ps1",
    "manuscript/recoil-bessel-limit/main.tex", "notes/stage4-recoil-bessel-limit.md",
    "literature/stage4-recoil-bessel-limit-sources.md",
    "results/stage4/recoil-bessel-limit/diagnostics.json",
    "results/stage4/recoil-bessel-limit/main.log",
    "results/stage4/recoil-bessel-limit/build-review.md",
    "output/pdf/de-sitter-recoil-bessel-limit.pdf",
    "code/recoil_global_probe.py", "results/stage4/recoil-global/diagnostics.json",
    "notes/stage4-recoil-global.md", "sources/dlmf-5.11.html",
}
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS["results/stage4/recoil-bessel-limit/evidence.json"]))
REQUIRED.append("results/stage4/recoil-bessel-limit/evidence.json")


CG_ROOT_HASH_RECORDS["results/stage4/recoil-bessel-cancellation/diagnostics.json"] = {
    "code/recoil_bessel_cancellation.py", "code/recoil_bessel_limit.py",
    "results/stage4/recoil-bessel-limit/diagnostics.json",
    "results/stage4/recoil-global/diagnostics.json", "sources/dlmf-5.11.html",
}
CG_ROOT_HASH_RECORDS["results/stage4/recoil-bessel-cancellation/evidence.json"] = {
    "code/recoil_bessel_cancellation.py", "code/build_bessel_cancellation_note.ps1",
    "manuscript/recoil-bessel-cancellation/main.tex",
    "notes/stage4-recoil-bessel-cancellation.md",
    "literature/stage4-recoil-bessel-cancellation-sources.md",
    "results/stage4/recoil-bessel-cancellation/diagnostics.json",
    "results/stage4/recoil-bessel-cancellation/main.log",
    "results/stage4/recoil-bessel-cancellation/build-review.md",
    "output/pdf/de-sitter-recoil-bessel-cancellation.pdf",
    "code/recoil_bessel_limit.py",
    "results/stage4/recoil-bessel-limit/diagnostics.json",
    "results/stage4/recoil-global/diagnostics.json",
    "notes/stage4-recoil-bessel-limit.md", "sources/dlmf-5.11.html",
}
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS["results/stage4/recoil-bessel-cancellation/evidence.json"]))
REQUIRED.append("results/stage4/recoil-bessel-cancellation/evidence.json")

CG_ROOT_HASH_RECORDS["results/stage4/recoil-bessel-endpoints/diagnostics.json"] = {
    "code/recoil_bessel_endpoints.py", "code/recoil_bessel_limit.py",
    "notes/stage4-recoil-bessel-cancellation.md", "sources/dlmf-5.11.html",
}
CG_ROOT_HASH_RECORDS["results/stage4/recoil-bessel-endpoints/evidence.json"] = {
    "code/recoil_bessel_endpoints.py", "code/build_bessel_endpoints_note.ps1",
    "manuscript/recoil-bessel-endpoints/main.tex",
    "notes/stage4-recoil-bessel-endpoints.md",
    "literature/stage4-recoil-bessel-endpoints-sources.md",
    "results/stage4/recoil-bessel-endpoints/diagnostics.json",
    "results/stage4/recoil-bessel-endpoints/main.log",
    "results/stage4/recoil-bessel-endpoints/build-review.md",
    "output/pdf/de-sitter-recoil-bessel-endpoints.pdf",
    "code/recoil_bessel_limit.py", "notes/stage4-recoil-bessel-cancellation.md",
    "sources/dlmf-5.11.html",
}
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS["results/stage4/recoil-bessel-endpoints/evidence.json"]))
REQUIRED.append("results/stage4/recoil-bessel-endpoints/evidence.json")

CG_ROOT_HASH_RECORDS["results/stage4/recoil-bessel-joint/diagnostics.json"] = {
    "code/recoil_bessel_joint.py", "code/recoil_bessel_limit.py",
    "code/recoil_bessel_cancellation.py",
    "results/stage4/recoil-bessel-limit/diagnostics.json",
    "notes/stage4-recoil-bessel-endpoints.md", "sources/dlmf-5.11.html",
}
CG_ROOT_HASH_RECORDS["results/stage4/recoil-bessel-joint/evidence.json"] = {
    "code/recoil_bessel_joint.py", "code/build_bessel_joint_note.ps1",
    "manuscript/recoil-bessel-joint/main.tex",
    "notes/stage4-recoil-bessel-joint.md",
    "literature/stage4-recoil-bessel-joint-sources.md",
    "results/stage4/recoil-bessel-joint/diagnostics.json",
    "results/stage4/recoil-bessel-joint/initial-contour-grid.json",
    "results/stage4/recoil-bessel-joint/initial-contour-grid.py",
    "results/stage4/recoil-bessel-joint/main.log",
    "results/stage4/recoil-bessel-joint/build-review.md",
    "output/pdf/de-sitter-recoil-bessel-joint.pdf",
    "code/recoil_bessel_limit.py", "code/recoil_bessel_cancellation.py",
    "results/stage4/recoil-bessel-limit/diagnostics.json",
    "notes/stage4-recoil-bessel-endpoints.md", "sources/dlmf-5.11.html",
}
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS["results/stage4/recoil-bessel-joint/evidence.json"]))
REQUIRED.append("results/stage4/recoil-bessel-joint/evidence.json")

CG_ROOT_HASH_RECORDS["results/stage4/recoil-transverse/diagnostics.json"] = {
    "code/recoil_transverse.py", "code/recoil_bessel_joint.py", "code/recoil_bessel_limit.py",
    "code/recoil_internal.py", "sources/chandra-liska-post-2610.30247v1.html", "sources/dlmf-5.11.html",
}
CG_ROOT_HASH_RECORDS["results/stage4/recoil-transverse/consequences.json"] = {
    "code/transverse_consequences.py", "code/recoil_transverse.py", "code/recoil_bessel_joint.py",
    "results/stage4/recoil-transverse/diagnostics.json", "results/stage4/recoil-bessel-joint/diagnostics.json",
}
CG_ROOT_HASH_RECORDS["results/stage4/recoil-transverse/evidence.json"] = (
    CG_ROOT_HASH_RECORDS["results/stage4/recoil-transverse/diagnostics.json"]
    | CG_ROOT_HASH_RECORDS["results/stage4/recoil-transverse/consequences.json"] | {
    "code/build_transverse_note.ps1", "manuscript/recoil-transverse/main.tex",
    "notes/stage4-recoil-transverse.md", "literature/stage4-recoil-transverse-sources.md",
    "results/stage4/recoil-transverse/consequences.json", "results/stage4/recoil-transverse/main.log",
    "results/stage4/recoil-transverse/build-review.md", "output/pdf/de-sitter-recoil-transverse.pdf",
})
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS["results/stage4/recoil-transverse/evidence.json"]))
REQUIRED.append("results/stage4/recoil-transverse/evidence.json")

PAIRING_DIR = "results/stage4/recoil-pairing/"
CG_ROOT_HASH_RECORDS[PAIRING_DIR + "diagnostics.json"] = {
    "code/recoil_pairing.py", "code/recoil_internal.py", "code/recoil_vertex.py",
    "code/recoil_bessel_limit.py", "code/recoil_bessel_joint.py",
    "results/stage4/recoil-bessel-joint/diagnostics.json",
    "results/stage4/recoil-transverse/diagnostics.json",
    "results/stage4/recoil-transverse/consequences.json",
    "sources/chandra-liska-post-2610.30247v1.html",
    "sources/chen-stanford-tang-yang-2607.14042v1.html",
}
CG_ROOT_HASH_RECORDS[PAIRING_DIR + "evidence.json"] = (
    CG_ROOT_HASH_RECORDS[PAIRING_DIR + "diagnostics.json"] | {
    "code/build_pairing_note.ps1", "manuscript/recoil-pairing/main.tex",
    "notes/stage4-recoil-pairing.md", "literature/stage4-recoil-pairing-sources.md",
    "output/pdf/de-sitter-recoil-pairing.pdf",
    *[PAIRING_DIR + name for name in ["diagnostics.json", "initial-window-grid.json",
      "initial-window-grid.py", "main.log", "build-review.md"]],
})
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS[PAIRING_DIR + "evidence.json"]))
REQUIRED.extend([PAIRING_DIR + "evidence.json", "journal/2026-10-31.md"])

POLES_DIR = "results/stage4/recoil-poles/"
CG_ROOT_HASH_RECORDS[POLES_DIR + "diagnostics.json"] = {
    "code/recoil_poles.py", "code/validate_recoil_poles.py",
    "code/recoil_internal.py", "code/recoil_vertex.py",
    "sources/chandra-liska-post-2610.30247v1.html",
    "sources/chen-stanford-tang-yang-2607.14042v1.html",
}
CG_ROOT_HASH_RECORDS[POLES_DIR + "evidence.json"] = (
    CG_ROOT_HASH_RECORDS[POLES_DIR + "diagnostics.json"] | {
    "code/build_poles_note.ps1", "manuscript/recoil-poles/main.tex",
    "notes/stage4-recoil-poles.md", "literature/stage4-recoil-poles-sources.md",
    "output/pdf/de-sitter-recoil-poles.pdf",
    *[POLES_DIR + name for name in ["diagnostics.json", "main.log", "build-review.md"]],
})
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS[POLES_DIR + "evidence.json"]))
REQUIRED.append(POLES_DIR + "evidence.json")

SPECTRAL_DIR = "results/stage4/recoil-spectral/"
CG_ROOT_HASH_RECORDS[SPECTRAL_DIR + "diagnostics.json"] = {
    "code/recoil_spectral.py", "code/recoil_poles.py", "code/recoil_internal.py", "code/recoil_vertex.py",
    "sources/chandra-liska-post-2610.30247v1.html", "sources/chen-stanford-tang-yang-2607.14042v1.html",
}
CG_ROOT_HASH_RECORDS[SPECTRAL_DIR + "evidence.json"] = (
    CG_ROOT_HASH_RECORDS[SPECTRAL_DIR + "diagnostics.json"] | {
    "code/build_spectral_note.ps1", "manuscript/recoil-spectral/main.tex",
    "notes/stage4-recoil-spectral.md", "literature/stage4-recoil-spectral-sources.md",
    "output/pdf/de-sitter-recoil-spectral.pdf",
    *[SPECTRAL_DIR + n for n in ["diagnostics.json", "initial-diagnostics.json",
      "initial-recoil-spectral.py", "main.log", "build-review.md"]],
})
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS[SPECTRAL_DIR + "evidence.json"]))
REQUIRED.append(SPECTRAL_DIR + "evidence.json")

IDENTITY_DIR = "results/stage4/recoil-identity/"
CG_ROOT_HASH_RECORDS[IDENTITY_DIR + "diagnostics.json"] = {
    "code/recoil_identity.py", "sources/collier-eberhardt-zhang-2304.13650v3.html",
    "sources/chandra-liska-post-2610.30247v1.html", "sources/collier-et-al-2501.01486v1.html",
    "notes/stage4-recoil-spectral.md",
}
CG_ROOT_HASH_RECORDS[IDENTITY_DIR + "evidence.json"] = (
    CG_ROOT_HASH_RECORDS[IDENTITY_DIR + "diagnostics.json"] | {
    "code/build_identity_note.ps1", "manuscript/recoil-identity/main.tex",
    "notes/stage4-recoil-identity.md", "literature/stage4-recoil-identity-sources.md",
    "output/pdf/de-sitter-recoil-identity.pdf",
    *[IDENTITY_DIR + n for n in ["diagnostics.json", "main.log", "build-review.md"]],
})
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS[IDENTITY_DIR + "evidence.json"]))
REQUIRED.append(IDENTITY_DIR + "evidence.json")

RACAH_DIR = "results/stage4/global-racah/"
RACAH_SOURCES = {
    "sources/derkachov-spiridonov-1711.07073v2.html",
    "sources/liu-et-al-1808.00612v3.html",
    "sources/karateev-kravchuk-simmons-duffin-1809.05111v1.html",
}
REQUIRED_SOURCE_ARCHIVES.update(RACAH_SOURCES)
CG_ROOT_HASH_RECORDS[RACAH_DIR + "diagnostics.json"] = RACAH_SOURCES | {
    "code/global_racah.py", "code/validate_global_racah.py",
    "code/recoil_poles.py", "code/recoil_internal.py", "code/recoil_vertex.py",
}
CG_ROOT_HASH_RECORDS[RACAH_DIR + "evidence.json"] = (
    CG_ROOT_HASH_RECORDS[RACAH_DIR + "diagnostics.json"] | {
    "code/build_global_racah_note.ps1", "manuscript/global-racah/main.tex",
    "notes/stage4-global-racah.md", "literature/stage4-global-racah-sources.md",
    "output/pdf/de-sitter-global-racah.pdf",
    *[RACAH_DIR + n for n in ["diagnostics.json", "main.log", "build-review.md",
                              "initial-driver.py.txt", "initial-exception.txt"]],
})
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS[RACAH_DIR + "evidence.json"]))
REQUIRED.append(RACAH_DIR + "evidence.json")

DS3_DIR = "results/stage4/ds3-local-vertex/"
DS3_SOURCES = {
    "sources/lagogiannis-maloney-wang-1106.2846v1.html",
    "sources/marolf-morrison-0810.5163v1.html", "sources/dlmf-16.4.html",
}
REQUIRED_SOURCE_ARCHIVES.update(DS3_SOURCES)
CG_ROOT_HASH_RECORDS[DS3_DIR + "diagnostics.json"] = DS3_SOURCES | {
    "code/ds3_local_vertex.py", "code/validate_ds3_local_vertex.py",
    "sources/derkachov-spiridonov-1711.07073v2.html",
    "sources/bros-et-al-0901.4223v1.html", "sources/bros-et-al-0901.4223v1.pdf",
    "sources/kolchmeyer-liu-2411.08090v1.html",
}
CG_ROOT_HASH_RECORDS[DS3_DIR + "clp-comparison.json"] = {
    "code/compare_ds3_local_clp.py", "code/ds3_local_vertex.py", "code/recoil_poles.py",
    "notes/stage4-recoil-poles.md", "notes/stage4-recoil-spectral.md",
}
CG_ROOT_HASH_RECORDS[DS3_DIR + "evidence.json"] = (
    CG_ROOT_HASH_RECORDS[DS3_DIR + "diagnostics.json"] |
    CG_ROOT_HASH_RECORDS[DS3_DIR + "clp-comparison.json"] | {
    "code/build_ds3_local_vertex_note.ps1", "manuscript/ds3-local-vertex/main.tex",
    "notes/stage4-ds3-local-vertex.md", "literature/stage4-ds3-local-vertex-sources.md",
    "output/pdf/de-sitter-ds3-local-vertex.pdf",
    *[DS3_DIR + n for n in ["diagnostics.json", "clp-comparison.json", "main.log",
                           "initial-tex.log", "build-review.md"]],
})
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS[DS3_DIR + "evidence.json"]))
REQUIRED.append(DS3_DIR + "evidence.json")

MASS_DIR = "results/stage4/ds3-mass-spectrum/"
MASS_INPUTS = {
    "code/ds3_mass_spectrum.py", "code/validate_ds3_mass_spectrum.py", "code/ds3_local_vertex.py",
    "sources/lagogiannis-maloney-wang-1106.2846v1.html",
    "sources/marolf-morrison-0810.5163v1.html", "sources/kolchmeyer-liu-2411.08090v1.html",
    "notes/stage4-ds3-local-vertex.md",
}
CG_ROOT_HASH_RECORDS[MASS_DIR + "initial-diagnostics.json"] = MASS_INPUTS
CG_ROOT_HASH_RECORDS[MASS_DIR + "diagnostics.json"] = MASS_INPUTS | {
    "code/complete_ds3_mass_diagnostics.py", MASS_DIR + "initial-diagnostics.json",
}
CG_ROOT_HASH_RECORDS[MASS_DIR + "evidence.json"] = (
    CG_ROOT_HASH_RECORDS[MASS_DIR + "diagnostics.json"] | {
    "code/build_ds3_mass_spectrum_note.ps1", "manuscript/ds3-mass-spectrum/main.tex",
    "notes/stage4-ds3-mass-spectrum.md", "literature/stage4-ds3-mass-spectrum-sources.md",
    "output/pdf/de-sitter-ds3-mass-spectrum.pdf",
    *[MASS_DIR + n for n in ["diagnostics.json", "main.log", "build-review.md"]],
})
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS[MASS_DIR + "evidence.json"]))
REQUIRED.append(MASS_DIR + "evidence.json")

GROUP_DIR = "results/stage4/ds3-group-pairing/"
GROUP_SOURCE = "sources/yang-zhang-zheng-2505.08647v1.html"
REQUIRED_SOURCE_ARCHIVES.add(GROUP_SOURCE)
CG_ROOT_HASH_RECORDS[GROUP_DIR + "diagnostics.json"] = {
    "code/ds3_group_pairing.py", "code/validate_ds3_group_pairing.py",
    "code/ds3_local_vertex.py", "code/ds3_mass_spectrum.py",
    "code/recoil_spectral.py", "code/recoil_poles.py",
    "sources/marolf-yu-2610.24073v1.html", "sources/marolf-morrison-0810.5163v1.html",
    "sources/kolchmeyer-liu-2411.08090v1.html", "sources/derkachov-spiridonov-1711.07073v2.html",
}
CG_ROOT_HASH_RECORDS[GROUP_DIR + "offdiagonal.json"] = {
    "code/validate_ds3_group_offdiagonal.py", "code/ds3_group_pairing.py",
    "code/ds3_local_vertex.py",
}
CG_ROOT_HASH_RECORDS[GROUP_DIR + "evidence.json"] = (
    CG_ROOT_HASH_RECORDS[GROUP_DIR + "diagnostics.json"]
    | CG_ROOT_HASH_RECORDS[GROUP_DIR + "offdiagonal.json"] | {
        GROUP_SOURCE, "code/audit_ds3_group_pairing.py", "code/build_ds3_group_pairing_note.ps1",
        "manuscript/ds3-group-pairing/main.tex", "notes/stage4-ds3-group-pairing.md",
        "literature/stage4-ds3-group-pairing-sources.md", "output/pdf/de-sitter-ds3-group-pairing.pdf",
        *[GROUP_DIR+n for n in ["diagnostics.json","offdiagonal.json","main.log","build-review.md"]],
})
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS[GROUP_DIR + "evidence.json"]))
REQUIRED.append(GROUP_DIR + "evidence.json")

TWO_DIR = "results/stage4/ds3-two-matter/"
CG_ROOT_HASH_RECORDS[TWO_DIR + "diagnostics.json"] = {
    "code/ds3_two_matter.py", "code/validate_ds3_two_matter.py",
    "code/ds3_group_pairing.py", "code/ds3_local_vertex.py", "code/ds3_mass_spectrum.py",
    "sources/derkachov-spiridonov-1711.07073v2.html", "sources/kolchmeyer-liu-2411.08090v1.html",
    "sources/marolf-morrison-0810.5163v1.html", "sources/marolf-yu-2610.24073v1.html",
}
CG_ROOT_HASH_RECORDS[TWO_DIR + "angular.json"] = {
    "code/validate_ds3_two_matter_angular.py", "code/ds3_two_matter.py",
    "code/validate_ds3_group_offdiagonal.py", "code/ds3_local_vertex.py", "code/ds3_group_pairing.py",
}
CG_ROOT_HASH_RECORDS[TWO_DIR + "evidence.json"] = (
    CG_ROOT_HASH_RECORDS[TWO_DIR+"diagnostics.json"] | CG_ROOT_HASH_RECORDS[TWO_DIR+"angular.json"] | {
        "code/audit_ds3_two_matter.py", "code/build_ds3_two_matter_note.ps1",
        "manuscript/ds3-two-matter/main.tex", "notes/stage4-ds3-two-matter.md",
        "literature/stage4-ds3-two-matter-sources.md", "output/pdf/de-sitter-ds3-two-matter.pdf",
        *[TWO_DIR+n for n in ["diagnostics.json","angular.json","main.log","initial-tex.log","build-review.md"]],
})
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS[TWO_DIR+"evidence.json"]))
REQUIRED.append(TWO_DIR+"evidence.json")

CROSS_DICT_DIR = "results/stage4/ds3-crossing-dictionary/"
CG_ROOT_HASH_RECORDS[CROSS_DICT_DIR+"evidence.json"] = {
    "code/ds3_crossing_dictionary.py", "code/ds3_crossing_tail.py",
    "code/validate_ds3_crossing_dictionary.py", "code/validate_ds3_crossing_geometry.py",
    "code/run_ds3_crossing_euclidean.py", "code/audit_ds3_crossing_dictionary.py",
    "code/build_ds3_crossing_dictionary.ps1", "code/global_racah.py",
    "code/ds3_local_vertex.py", "code/ds3_two_matter.py", "code/ds3_group_pairing.py",
    "code/ds3_mass_spectrum.py", "manuscript/ds3-crossing-dictionary/main.tex",
    "notes/stage4-ds3-crossing-dictionary.md",
    "literature/stage4-ds3-crossing-dictionary-sources.md",
    "output/pdf/de-sitter-ds3-crossing-dictionary.pdf",
    "sources/liu-et-al-1808.00612v3.html",
    "sources/karateev-kravchuk-simmons-duffin-1809.05111v1.html",
    "sources/derkachov-spiridonov-1711.07073v2.html",
    *[CROSS_DICT_DIR+n for n in ["diagnostics.json", "geometry.json", "euclidean.json",
        "euclidean-refined.json", "executed-kernel.py.txt", "source-correction.json",
        "build-review.md", "main.log"]],
}
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS[CROSS_DICT_DIR+"evidence.json"]))
REQUIRED.append(CROSS_DICT_DIR+"evidence.json")

PACKET4_DIR = "results/stage4/ds3-crossed-packet/"
PACKET4_LABELS = ["pilot-n4-L8", "pilot-precision-n4-L8", "refined-n6-L10",
                  "precision-n6-L12", "refined-n8-L12"]
CG_ROOT_HASH_RECORDS[PACKET4_DIR+"evidence.json"] = {
    "code/ds3_racah_fast.py", "code/ds3_crossed_packet.py", "code/ds3_packet_tail.py",
    "code/ds3_cutoff_witness.py", "code/run_ds3_crossed_packet.py",
    "code/validate_ds3_racah_fast.py", "code/validate_ds3_packet_bounds.py",
    "code/replay_ds3_crossed_packet.py", "code/audit_ds3_crossed_packet.py",
    "code/build_ds3_crossed_packet.ps1", "code/ds3_crossing_dictionary.py",
    "code/global_racah.py", "code/ds3_local_vertex.py", "code/ds3_group_pairing.py",
    "code/ds3_mass_spectrum.py", "code/ds3_two_matter.py",
    "manuscript/ds3-crossed-packet/main.tex", "manuscript/ds3-crossed-packet/numerics.tex",
    "notes/stage4-ds3-crossed-packet.md", "literature/stage4-ds3-crossed-packet-sources.md",
    "output/pdf/de-sitter-ds3-crossed-packet.pdf",
    "sources/liu-et-al-1808.00612v3.html", "sources/derkachov-spiridonov-1711.07073v2.html",
    "sources/karateev-kravchuk-simmons-duffin-1809.05111v1.html",
    *[PACKET4_DIR+label+ext for label in PACKET4_LABELS for ext in [".json",".npz"]],
    *[PACKET4_DIR+n for n in ["summary.json","bounds-validation.json","accelerator-validation.json",
                             "build-review.md","main.log"]],
}
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS[PACKET4_DIR+"evidence.json"]))
REQUIRED.append(PACKET4_DIR+"evidence.json")

BOSE_DIR = "results/stage4/ds3-bose/"
BOSE_LABELS = ["pilot-n4-L8", "fixed-n4-L12", "refined-n6-L12"]
CG_ROOT_HASH_RECORDS[BOSE_DIR+"evidence.json"] = {
    "code/ds3_bose_packet.py", "code/validate_ds3_bose.py", "code/ds3_bose_precision.py",
    "code/check_ds3_bose_nodes.py", "code/replay_ds3_bose.py", "code/audit_ds3_bose.py",
    "code/build_ds3_bose.ps1", "code/ds3_racah_fast.py", "code/ds3_crossing_dictionary.py",
    "code/ds3_packet_tail.py", "code/ds3_crossed_packet.py", "code/replay_ds3_crossed_packet.py",
    "code/global_racah.py", "manuscript/ds3-bose/main.tex", "manuscript/ds3-bose/numerics.tex",
    "notes/stage4-ds3-bose.md", "literature/stage4-ds3-bose-sources.md",
    "output/pdf/de-sitter-ds3-bose.pdf", "sources/kolchmeyer-liu-2411.08090v1.html",
    "sources/liu-et-al-1808.00612v3.html",
    *[BOSE_DIR+label+ext for label in BOSE_LABELS for ext in [".json", ".npz"]],
    *[BOSE_DIR+n for n in ["summary.json", "validation.json", "precision-repair.json",
                          "weighted-node-checks.json", "build-review.md", "main.log", "initial-tex.log"]],
}
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS[BOSE_DIR+"evidence.json"]))
REQUIRED.append(BOSE_DIR+"evidence.json")

# Post-Bose stages retain their original failed outputs and historical PDFs.
for stage, main_code in [
    ("ds3-collision", "ds3_controlled_omega"),
    ("ds3-full-tail", "ds3_full_tail_bound"),
    ("ds3-shadow-duality", "ds3_shadow_duality"),
    ("ds3-peak-limit", "ds3_peak_normalization"),
    ("ds3-transverse-phase", "ds3_transverse_phase"),
    ("ds3-uniform-local", "ds3_uniform_local"),
]:
    evidence_path = "results/stage4/"+stage+"/evidence.json"
    expected = {"code/"+main_code+".py", "manuscript/"+stage+"/main.tex",
                "output/pdf/de-sitter-"+stage+".pdf", "notes/stage4-"+stage+".md",
                "literature/stage4-"+stage+"-sources.md", "results/stage4/"+stage+"/study.json"}
    CG_ROOT_HASH_RECORDS[evidence_path] = expected
    REQUIRED.extend(sorted(expected))
    REQUIRED.append(evidence_path)
REQUIRED.extend(["code/audit_ds3_recent_stages.py", "results/stage4/ds3-collision/completion.md",
                 "results/stage4/ds3-collision/controlled-validation.json",
                 "results/stage4/ds3-transverse-phase/refinement.json"])

PHYSICAL_LOCAL_DIR = "results/stage4/ds3-physical-local/"
CG_ROOT_HASH_RECORDS[PHYSICAL_LOCAL_DIR+"evidence.json"] = {
    "code/ds3_physical_local.py", "code/check_ds3_physical_local.py",
    "code/audit_ds3_physical_local.py", "code/build_ds3_physical_local.ps1",
    "manuscript/ds3-physical-local/main.tex", "output/pdf/de-sitter-ds3-physical-local.pdf",
    "notes/stage4-ds3-physical-local.md", "literature/stage4-ds3-physical-local-sources.md",
    PHYSICAL_LOCAL_DIR+"diagnostics.json", PHYSICAL_LOCAL_DIR+"audit.json",
    PHYSICAL_LOCAL_DIR+"build-review.md", PHYSICAL_LOCAL_DIR+"main.log",
}
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS[PHYSICAL_LOCAL_DIR+"evidence.json"]))
REQUIRED.append(PHYSICAL_LOCAL_DIR+"evidence.json")

COMPACT_GLOBAL_DIR = "results/stage4/ds3-compact-global/"
CG_ROOT_HASH_RECORDS[COMPACT_GLOBAL_DIR+"evidence.json"] = {
    "code/ds3_compact_global.py", "code/run_ds3_compact_global.py",
    "code/diagnose_ds3_compact_collision.py", "code/audit_ds3_compact_global.py",
    "code/audit_ds3_compact_global_refined.py", "code/build_ds3_compact_global.ps1",
    "manuscript/ds3-compact-global/main.tex", "output/pdf/de-sitter-ds3-compact-global.pdf",
    "notes/stage4-ds3-compact-global.md", "literature/stage4-ds3-compact-global-sources.md",
    COMPACT_GLOBAL_DIR+"study.json", COMPACT_GLOBAL_DIR+"collision-diagnosis.json",
    COMPACT_GLOBAL_DIR+"audit.json", COMPACT_GLOBAL_DIR+"audit-refined.json",
    COMPACT_GLOBAL_DIR+"runs.jsonl", COMPACT_GLOBAL_DIR+"build-review.md",
    COMPACT_GLOBAL_DIR+"main.log",
}
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS[COMPACT_GLOBAL_DIR+"evidence.json"]))
REQUIRED.append(COMPACT_GLOBAL_DIR+"evidence.json")


FOCK_VACUUM_DIR = "results/stage4/ds3-fock-vacuum/"
CG_ROOT_HASH_RECORDS[FOCK_VACUUM_DIR+"evidence.json"] = {
    "code/"+name for name in ["ds3_fock_vacuum.py", "check_ds3_fock_vacuum.py",
    "recover_ds3_fock_vacuum.py", "audit_ds3_fock_vacuum.py", "audit_ds3_fock_vacuum_refined.py",
    "ds3_gaussian_field.py", "check_ds3_gaussian_field.py", "audit_ds3_gaussian_field.py",
    "record_ds3_fock_vacuum.py", "build_ds3_fock_vacuum.ps1"]
} | {"manuscript/ds3-fock-vacuum/main.tex", "output/pdf/de-sitter-ds3-fock-vacuum.pdf",
     "notes/stage4-ds3-fock-vacuum.md", "literature/stage4-ds3-fock-vacuum-sources.md"} | {
    FOCK_VACUUM_DIR+name for name in ["vacuum-runs.jsonl", "diagnostics-recovered.json",
    "audit.json", "audit-refined.json", "gaussian-runs.jsonl", "gaussian.json",
    "gaussian-audit.json", "preservation.json", "build-review.md", "main.log"]
}
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS[FOCK_VACUUM_DIR+"evidence.json"]))
REQUIRED.append(FOCK_VACUUM_DIR+"evidence.json")

GAUSSIAN_EIKONAL_DIR = "results/stage4/gaussian-eikonal/"
CG_ROOT_HASH_RECORDS[GAUSSIAN_EIKONAL_DIR+"evidence.json"] = {
    "code/"+name for name in ["gaussian_eikonal.py", "run_gaussian_eikonal.py", "check_gaussian_eikonal.py",
    "refine_gaussian_eikonal.py", "record_gaussian_eikonal.py", "build_gaussian_eikonal.ps1"]
} | {"manuscript/gaussian-eikonal/main.tex", "output/pdf/de-sitter-gaussian-eikonal.pdf",
     "notes/stage4-gaussian-eikonal.md", "literature/stage4-gaussian-eikonal-sources.md"} | {
    GAUSSIAN_EIKONAL_DIR+name for name in ["study.json", "runs.jsonl", "diagnostics.json",
    "refinement.json", "preservation.json", "build-review.md", "main.log"]}
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS[GAUSSIAN_EIKONAL_DIR+"evidence.json"]))
REQUIRED.append(GAUSSIAN_EIKONAL_DIR+"evidence.json")
for file, expected in [("study.json", {"code/gaussian_eikonal.py", "code/run_gaussian_eikonal.py", "code/conformal_shock.py", "sources/chen-stanford-tang-yang-2607.14042v1.html"}),
                       ("diagnostics.json", {"code/gaussian_eikonal.py", "code/check_gaussian_eikonal.py", GAUSSIAN_EIKONAL_DIR+"study.json"}),
                       ("refinement.json", {"code/refine_gaussian_eikonal.py", "code/gaussian_eikonal.py", GAUSSIAN_EIKONAL_DIR+"study.json", GAUSSIAN_EIKONAL_DIR+"diagnostics.json"})]:
    CG_ROOT_HASH_RECORDS[GAUSSIAN_EIKONAL_DIR+file] = expected

PRINCIPAL_EIKONAL_DIR = "results/stage4/principal-eikonal/"
CG_ROOT_HASH_RECORDS[PRINCIPAL_EIKONAL_DIR+"evidence.json"] = {
    "code/"+name for name in ["principal_eikonal.py", "run_principal_eikonal.py", "check_principal_eikonal.py",
    "refine_principal_eikonal.py", "late_principal_eikonal.py", "audit_principal_eikonal.py",
    "record_principal_eikonal.py", "build_principal_eikonal.ps1"]
} | {"manuscript/principal-eikonal/main.tex", "output/pdf/de-sitter-principal-eikonal.pdf",
     "notes/stage4-principal-eikonal.md", "literature/stage4-principal-eikonal-sources.md"} | {
    PRINCIPAL_EIKONAL_DIR+name for name in ["study.json", "runs.jsonl", "moments.json", "diagnostics.json",
    "refinement.json", "refinement-runs.jsonl", "late.json", "late-runs.jsonl", "audit.json",
    "preservation.json", "source-record.json", "build-review.md", "main.log", "build-initial-overflow.log",
    "preservation-initial-failure.txt"]}
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS[PRINCIPAL_EIKONAL_DIR+"evidence.json"]))
REQUIRED.append(PRINCIPAL_EIKONAL_DIR+"evidence.json")
for file, expected in [
    ("study.json", {"code/principal_eikonal.py", "code/run_principal_eikonal.py", "sources/chen-stanford-tang-yang-2607.14042v1.html"}),
    ("diagnostics.json", {"code/principal_eikonal.py", "code/check_principal_eikonal.py", PRINCIPAL_EIKONAL_DIR+"moments.json"}),
    ("refinement.json", {"code/principal_eikonal.py", "code/check_principal_eikonal.py", "code/refine_principal_eikonal.py", PRINCIPAL_EIKONAL_DIR+"study.json"}),
    ("late.json", {"code/principal_eikonal.py", "code/check_principal_eikonal.py", "code/late_principal_eikonal.py", PRINCIPAL_EIKONAL_DIR+"moments.json"}),
    ("audit.json", {"code/audit_principal_eikonal.py"})]:
    CG_ROOT_HASH_RECORDS[PRINCIPAL_EIKONAL_DIR+file] = expected

MASS_FILTER_DIR = 'results/stage4/mass-filter/'
CG_ROOT_HASH_RECORDS[MASS_FILTER_DIR+"evidence.json"] = {'results/stage4/mass-filter/source-record.json', 'results/stage4/mass-filter/audit.json', 'manuscript/mass-filter/main.tex', 'results/stage4/mass-filter/main.log', 'results/stage4/mass-filter/packet.json', 'code/record_mass_filter.py', 'results/stage4/mass-filter/study.json', 'literature/stage4-mass-filter-sources.md', 'code/run_mass_filter.py', 'code/audit_mass_filter.py', 'results/stage4/mass-filter/preservation.json', 'results/stage4/mass-filter/build-review.md', 'code/packet_mass_filter.py', 'results/stage4/mass-filter/runs.jsonl', 'code/build_mass_filter.ps1', 'notes/stage4-mass-filter.md', 'code/mass_filter.py', 'code/check_mass_filter.py', 'output/pdf/de-sitter-mass-filter.pdf', 'results/stage4/mass-filter/packet-runs.jsonl', 'results/stage4/mass-filter/diagnostics.json'}
REQUIRED.extend(sorted(CG_ROOT_HASH_RECORDS[MASS_FILTER_DIR+"evidence.json"]))
REQUIRED.append(MASS_FILTER_DIR+"evidence.json")
CG_ROOT_HASH_RECORDS[MASS_FILTER_DIR+'study.json'] = {'code/mass_filter.py', 'code/principal_eikonal.py', 'code/run_mass_filter.py'}
CG_ROOT_HASH_RECORDS[MASS_FILTER_DIR+'diagnostics.json'] = {'code/check_mass_filter.py', 'code/mass_filter.py', 'code/principal_eikonal.py'}
CG_ROOT_HASH_RECORDS[MASS_FILTER_DIR+'packet.json'] = {'code/packet_mass_filter.py', 'code/principal_eikonal.py'}
CG_ROOT_HASH_RECORDS[MASS_FILTER_DIR+'audit.json'] = {'code/audit_mass_filter.py'}

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--snapshot", action="store_true", help="Requires a new --label; never overwrites")
    parser.add_argument("--label", help="Plain filename stem, e.g. 2026-10-29-operator-bridge-finite-g")
    args = parser.parse_args()
    if args.snapshot and not args.label:
        parser.error("--snapshot requires an explicit --label")
    if args.label and not re.fullmatch(r"[a-zA-Z0-9][a-zA-Z0-9_-]{0,100}", args.label):
        parser.error("label must contain only letters, digits, hyphens or underscores")
    target = ROOT/"snapshots"/(str(args.label)+".zip")
    if args.snapshot and (target.exists() or target.with_suffix(".sha256").exists()):
        parser.error("snapshot or checksum already exists; choose a new label")

    errors, warnings = [], []
    for name in REQUIRED:
        if not (ROOT/name).is_file():
            errors.append(f"missing: {name}")
    checked_links = 0
    for path in ROOT.rglob("*.md"):
        if ".git" in path.parts:
            continue
        markdown = path.read_text(encoding="utf-8-sig")
        # Exclude mathematical [A](N+1) products and examples before parsing links.
        markdown = re.sub(r"\x60{3}.*?\x60{3}", "", markdown, flags=re.S)
        markdown = re.sub(r"\\\[.*?\\\]", "", markdown, flags=re.S)
        markdown = re.sub(r"\\\(.*?\\\)", "", markdown, flags=re.S)
        for link in re.findall(r"(?<!\\)\[[^\]\n]*\]\(([^)\n]+)\)", markdown):
            if link.startswith(("http:", "https:", "#", "mailto:")):
                continue
            link = link.strip("<>").split("#", 1)[0]
            if link and not (path.parent/link).exists():
                errors.append(f"broken local link: {path.relative_to(ROOT)} -> {link}")
            checked_links += 1

    def read_json(path):
        try:
            return json.loads(path.read_text(encoding="utf-8-sig"))
        except (OSError, ValueError) as exc:
            errors.append(f"unreadable JSON: {path.relative_to(ROOT)}: {exc}")
            return {}

    scientific, completed_studies, checked_hashes = {}, {}, 0
    for rel, label, is_validation in RECORDS:
        path = ROOT/rel
        data = read_json(path)
        packet_schema = label.startswith("stage3-packets")
        crossing_schema = label == "stage3-crossing"
        rate_schema = label in ROOT_HASH_RECORDS
        regular_schema = label == "stage3-crossing-regular"
        if crossing_schema:
            hashes = data.get("provenance", {})
        elif packet_schema:
            hashes = data.get("metadata", {}).get("code_sha256", {})
        elif rate_schema:
            hashes = data.get("metadata", {}).get("hashes", {})
        else:
            hashes = data.get("code_sha256", {})
        if not isinstance(hashes, dict) or not hashes:
            errors.append(f"missing code/provenance hashes: {rel}")
            hashes = {}
        if crossing_schema:
            actual_paths = {name.replace("\\", "/") for name in hashes}
            for missing in sorted(CROSSING_PROVENANCE - actual_paths):
                errors.append(f"missing crossing provenance: {missing}")
        if rate_schema:
            actual_paths = {name.replace("\\", "/") for name in hashes}
            for missing in sorted(EXPECTED_RATE_HASHES - actual_paths):
                errors.append(f"missing rate provenance: {label} / {missing}")
        for filename, expected in hashes.items():
            filename = filename.replace("\\", "/")
            source = ROOT/filename if packet_schema or crossing_schema or rate_schema else ROOT/"code"/filename
            if not source.is_file() or digest(source) != expected:
                errors.append(f"stale result: {label} / {filename}")
            checked_hashes += 1
        inputs = data.get("input_sha256", {})
        if not isinstance(inputs, dict):
            errors.append(f"invalid input hashes: {rel}")
            inputs = {}
        if (regular_schema or label in COMPLETED_STUDIES or label in EXPECTED_INPUT_HASHES) and not inputs:
            errors.append(f"missing input hashes: {rel}")
        input_paths = {name.replace("\\", "/") for name in inputs}
        for missing in sorted(EXPECTED_INPUT_HASHES.get(label, set()) - input_paths):
            errors.append(f"missing expected input hash: {label} / {missing}")
        for filename, expected in inputs.items():
            source = ROOT/filename.replace("\\", "/")
            if not source.is_file() or digest(source) != expected:
                errors.append(f"stale input: {label} / {filename}")
            checked_hashes += 1
        if label in COMPLETED_STUDIES:
            rows = data.get("rows")
            finished = (data.get("status") == "completed"
                        and isinstance(rows, list) and len(rows) > 0
                        and all(isinstance(row, dict) and "failure" not in row for row in rows)
                        and data.get("exception") is None)
            if not finished:
                errors.append(f"study not completed: {label}")
            if label == "packet-crossing-study" and data.get("frozen_point_archive") != PACKET_RAW_ARCHIVE:
                errors.append("packet study frozen archive does not match the required exported points")
            completed_studies[label] = {"status": data.get("status"),
                                       "rows": len(rows) if isinstance(rows, list) else None}
        if is_validation:
            if crossing_schema:
                summary = data.get("summary", {})
                count, failed = summary.get("checks"), summary.get("failed")
                details, exceptions = data.get("checks"), data.get("exceptions")
                exception_count = summary.get("evaluation_exceptions")
                points = data.get("points")
                passed = (data.get("schema") == "observer-crossing-validation-v1"
                          and type(count) is int and count > 0
                          and failed == 0 and summary.get("passed") == count
                          and summary.get("all_passed") is True
                          and isinstance(details, list) and len(details) == count
                          and all(isinstance(row, dict) and row.get("passed") is True
                                  for row in details)
                          and type(exception_count) is int and exception_count == 0
                          and isinstance(exceptions, list) and len(exceptions) == 0
                          and isinstance(points, list) and len(points) > 0
                          and all(isinstance(row, dict) and "failure" not in row for row in points))
                status = "passed" if passed else "failed"
                if not (type(exception_count) is int and exception_count == 0
                        and isinstance(exceptions, list) and len(exceptions) == 0):
                    errors.append("crossing evaluation exceptions are nonempty or inconsistent")
            elif rate_schema:
                summary = data.get("summary", {})
                count, failed = summary.get("total"), summary.get("failed")
                details = data.get("checks")
                passed = (type(count) is int and count > 0 and type(failed) is int and failed == 0
                          and summary.get("passed") == count and isinstance(details, list)
                          and len(details) == count
                          and all(isinstance(row, dict) and row.get("passed") is True
                                  for row in details))
                status = "passed" if passed else "failed"
            elif packet_schema:
                count, failed = data.get("total"), data.get("failed")
                details = data.get("checks")
                passed = (isinstance(count, int) and count > 0 and failed == 0
                          and data.get("passed") == count and isinstance(details, list)
                          and len(details) == count
                          and all(isinstance(row, dict) and row.get("passed") is True
                                  for row in details))
                status = "passed" if passed else "failed"
            else:
                count, failed, status = data.get("checks"), data.get("failed"), data.get("status")
                details = data.get("check_details")
            if regular_schema:
                exception_count = data.get("evaluation_exceptions")
                exceptions = data.get("exceptions")
                if not (type(exception_count) is int and exception_count == 0
                        and isinstance(exceptions, list) and len(exceptions) == 0):
                    errors.append("regular crossing evaluation exceptions are nonempty or inconsistent")
                    status = "failed"
            if (regular_schema or label == "conditional-crossing" or label in PACKET_VALIDATIONS) and details is None:
                errors.append(f"missing validation details: {label}")
            if status != "passed" or type(failed) is not int or failed != 0 or type(count) is not int or count <= 0:
                errors.append(f"scientific validation not passed: {label}")
            if details is not None:
                if (not isinstance(details, list) or len(details) != count
                        or any(not isinstance(row, dict) or row.get("passed") is not True
                               for row in details)):
                    errors.append(f"validation summary/detail mismatch: {label}")
            if label in EXPECTED_VALIDATION_COUNTS and count != EXPECTED_VALIDATION_COUNTS[label]:
                errors.append(f"unexpected validation count: {label} / {count}")
            scientific[label] = {"checks": count, "failed": failed, "status": status}
            if crossing_schema or regular_schema:
                scientific[label]["evaluation_exceptions"] = exception_count
                scientific[label]["exception_records"] = len(exceptions) if isinstance(exceptions, list) else None

    historical = []
    for rel in HISTORICAL_RECORDS:
        data = read_json(ROOT/rel)
        historical.append({"file": rel, "recorded_status": data.get("status"),
                           "recorded_failed": data.get("failed"),
                           "hash_policy": "retained historical bytes; old code hashes not compared to current code"})

    supporting_records = []
    for rel, required_hashes in CG_ROOT_HASH_RECORDS.items():
        data = read_json(ROOT/rel)
        hashes = data.get("hashes", {})
        if not isinstance(hashes, dict) or not hashes:
            errors.append(f"missing supporting provenance hashes: {rel}")
            hashes = {}
        actual_paths = {name.replace("\\", "/") for name in hashes}
        for missing in sorted(required_hashes - actual_paths):
            errors.append(f"missing supporting provenance: {rel} / {missing}")
        for filename, expected in hashes.items():
            source = ROOT/filename.replace("\\", "/")
            if not source.is_file() or digest(source) != expected:
                errors.append(f"stale supporting input: {rel} / {filename}")
            checked_hashes += 1
        supporting_records.append({"file": rel, "kind": "evidence/root-relative provenance",
                                   "recorded_hashes": len(hashes)})

    comparison_rel = "results/stage3/cross-limit/comparison.json"
    comparison = read_json(ROOT/comparison_rel)
    comparison_checks = comparison.get("checks")
    if (comparison.get("status") != "passed"
            or not isinstance(comparison_checks, list) or len(comparison_checks) != 19
            or any(not isinstance(row, dict) or row.get("passed") is not True
                   for row in comparison_checks)
            or len(comparison.get("rows", [])) != 6):
        errors.append("invalid projected comparison-model diagnostics")
    supporting_records.append({"file": comparison_rel,
                               "kind": "projected comparison model only; not exact 6j convergence or interval certification",
                               "checks": len(comparison_checks) if isinstance(comparison_checks, list) else None})

    endpoint_rel = "results/stage3/endpoint/diagnostics.json"
    endpoint = read_json(ROOT/endpoint_rel)
    endpoint_checks = endpoint.get("connection_and_precision_checks", [])
    if (endpoint.get("status") != "passed" or endpoint.get("checks") != 15
            or len(endpoint_checks) != 15
            or any(row.get("passed") is not True for row in endpoint_checks)
            or len(endpoint.get("single_hypergeometric_rows", [])) != 9
            or len(endpoint.get("four_branch_diagnostics", [])) != 3):
        errors.append("invalid finite-point endpoint diagnostics")
    supporting_records.append({"file": endpoint_rel,
                               "kind": "finite-point endpoint diagnostics, not distribution-limit certification",
                               "checks": len(endpoint_checks)})

    assembly_rel = "results/stage3/assembly/diagnostics.json"
    assembly = read_json(ROOT/assembly_rel)
    assembly_checks = assembly.get("checks", [])
    if (assembly.get("status") != "passed" or len(assembly_checks) != 21
            or any(row.get("passed") is not True for row in assembly_checks)
            or len(assembly.get("exact_kernel_rows", [])) != 8
            or len(assembly.get("leading_distribution_tests", [])) != 4):
        errors.append("invalid assembled-kernel finite-point diagnostics")
    supporting_records.append({"file": assembly_rel,
                               "kind": "finite-point algebra and asymptotic diagnostics; leading-kernel test only, not exact distribution quadrature",
                               "checks": len(assembly_checks)})

    logtime_rel = "results/stage3/logtime/diagnostics.json"
    logtime = read_json(ROOT/logtime_rel)
    logtime_checks = logtime.get("checks", [])
    if (logtime.get("status") != "passed" or len(logtime_checks) != 25
            or any(row.get("passed") is not True for row in logtime_checks)
            or len(logtime.get("coefficient_rows", [])) != 3
            or len(logtime.get("leading_distribution_rows", [])) != 12
            or len(logtime.get("plateau_rows", [])) != 4):
        errors.append("invalid logarithmic-time coefficient/PV diagnostics")
    supporting_records.append({"file": logtime_rel,
                               "kind": "coefficient and compact leading-expression PV diagnostics; not exact 6j or physical packet quadrature",
                               "checks": len(logtime_checks)})

    prep_rel = "results/stage3/preparation/diagnostics.json"
    prep = read_json(ROOT/prep_rel)
    prep_checks = prep.get("checks", [])
    if (prep.get("status") != "passed" or len(prep_checks) != 37
            or any(row.get("passed") is not True for row in prep_checks)
            or len(prep.get("factorization_rows", [])) != 3
            or len(prep.get("corner_rows", [])) != 6
            or any(len(prep.get("fixed_runs", {}).get(str(n), {}).get("rows", [])) != 8 for n in (24,40,64))):
        errors.append("invalid preparation curve/corner diagnostics")
    supporting_records.append({"file": prep_rel,
                               "kind": "empirical F limiting curve and Q local corner; not full Q curve, finite-mass or interval certification",
                               "checks": len(prep_checks)})

    moving_rel = "results/stage3/moving-curve/diagnostics.json"
    moving = read_json(ROOT/moving_rel)
    moving_checks = moving.get("checks", [])
    if (moving.get("status") != "passed" or len(moving_checks) != 122
            or any(row.get("passed") is not True for row in moving_checks)
            or len(moving.get("geometry_rows", [])) != 7
            or len(moving.get("curve_refinement_rows", [])) != 25):
        errors.append("invalid complete moving-threshold curve diagnostics")
    for rel in MOVING_RUNS:
        run = read_json(ROOT/rel)
        spectrum = run.get("spectrum", {})
        if (len(run.get("rows", [])) != 25
                or len(run.get("dense_curve", {}).get("T", [])) != 257
                or not (len(spectrum.get("chi", [])) == len(spectrum.get("weights", []))
                        == len(spectrum.get("JQ", [])) > 1000)):
            errors.append(f"incomplete moving-threshold raw run: {rel}")
    supporting_records.append({"file": moving_rel,
                               "kind": "actual complex limiting Q curves; empirical transverse/Fourier refinements, no finite-M or total interval certification",
                               "checks": len(moving_checks)})

    wkb_rel = "results/stage4/wkb-overlap/diagnostics.json"
    wkb = read_json(ROOT/wkb_rel)
    wkb_checks = wkb.get("checks", [])
    if (wkb.get("status") != "passed" or len(wkb_checks) != 42
            or any(row.get("passed") is not True for row in wkb_checks)
            or len(wkb.get("samples", [])) != 6
            or len(wkb.get("critical_scaling", [])) != 9
            or len(wkb.get("prefactors", [])) != 6):
        errors.append("invalid reduced WKB overlap diagnostics")
    supporting_records.append({"file": wkb_rel,
                               "kind": "reduced recoil and action-only integral checks; not finite-G amplitude certification",
                               "checks": len(wkb_checks)})

    gram_rel = "results/stage4/gram-prefactor/diagnostics.json"
    gram = read_json(ROOT/gram_rel)
    gram_checks = gram.get("checks", [])
    if (gram.get("status") != "passed" or len(gram_checks) != 83
            or any(row.get("passed") is not True for row in gram_checks)
            or any(len(gram.get(key, [])) != count for key, count in {
                "general_geometry": 5, "fixed_angle": 9, "pure_recoil": 6,
                "weighted_integrals": 6, "transverse_hessian": 9,
            }.items())):
        errors.append("invalid Gram and transverse Hessian diagnostics")
    supporting_records.append({"file": gram_rel,
                               "kind": "exact Gram algebra and classical Hessian checks; conditional local prefactor only, not full finite-G amplitude",
                               "checks": len(gram_checks)})

    vertex_rel = "results/stage4/recoil-vertex/diagnostics.json"
    vertex = read_json(ROOT/vertex_rel)
    vertex_checks = vertex.get("checks", [])
    if (vertex.get("status") != "passed" or len(vertex_checks) != 42
            or any(row.get("passed") is not True for row in vertex_checks)
            or any(len(vertex.get(key, [])) != count for key, count in {
                "gamma_rows": 4, "shift_rows": 6, "lattice_rows": 4,
                "box_rows": 9, "observer_rows": 3,
            }.items())):
        errors.append("invalid combined recoil vertex diagnostics")
    supporting_records.append({"file": vertex_rel,
                               "kind": "balanced Barnes moments, limiting vertices and local bounds; not full 6j or physical four-point normalization",
                               "checks": len(vertex_checks)})

    internal_rel = "results/stage4/recoil-internal/diagnostics.json"
    internal = read_json(ROOT/internal_rel)
    internal_checks = internal.get("checks", [])
    if (internal.get("status") != "passed" or len(internal_checks) != 34
            or any(row.get("passed") is not True for row in internal_checks)
            or any(len(internal.get(key, [])) != count for key, count in {
                "point_rows": 9, "integral_rows": 3, "tail_rows": 8,
            }.items())
            or not any("residue" in row and "cancellation_ratio" in row
                       for row in internal.get("integral_rows", []))):
        errors.append("invalid fixed-mass internal recoil diagnostics")
    supporting_records.append({"file": internal_rel,
                               "kind": "exact finite-kappa shifts and fixed-mass Barnes limit diagnostics; not physical external-energy integration",
                               "checks": len(internal_checks)})

    saddle_rel = "results/stage4/recoil-saddle/diagnostics.json"
    saddle = read_json(ROOT/saddle_rel)
    saddle_checks = saddle.get("checks", [])
    if (saddle.get("status") != "passed" or len(saddle_checks) != 57
            or any(row.get("passed") is not True for row in saddle_checks)
            or len(saddle.get("saddle_rows", [])) != 5
            or len(saddle.get("scaling_rows", [])) != 9):
        errors.append("invalid reflected Barnes saddle diagnostics")
    supporting_records.append({"file": saddle_rel,
                               "kind": "saddle algebra and local numerical comparisons; not full contour or physical OTOC certification",
                               "checks": len(saddle_checks)})

    local_rel = "results/stage4/recoil-local/diagnostics.json"
    local = read_json(ROOT/local_rel)
    local_checks = local.get("checks", [])
    if (local.get("status") != "passed" or len(local_checks) != 50
            or any(row.get("passed") is not True for row in local_checks)
            or len(local.get("angle_rows", [])) != 3
            or len(local.get("local_rows", [])) != 3):
        errors.append("invalid local Gamma and general-action diagnostics")
    initial_local = read_json(ROOT/"results/stage4/recoil-local/initial-norm-grid.json")
    if digest(ROOT/"results/stage4/recoil-local/initial-code.py.txt") != initial_local.get("hashes", {}).get("code/recoil_local.py"):
        errors.append("initial local failure does not match archived initial code")
    supporting_records.append({"file": local_rel,
                               "kind": "general angle algebra and actual localized Gamma integrals; not the full Barnes cycle",
                               "checks": len(local_checks)})

    global_rel = "results/stage4/recoil-global/diagnostics.json"
    global_data = read_json(ROOT/global_rel)
    global_checks = global_data.get("checks", [])
    global_rows = global_data.get("full_rows", [])
    if (global_data.get("status") != "passed" or len(global_checks) != 34
            or any(row.get("passed") is not True for row in global_checks)
            or len(global_rows) != 7
            or len(global_data.get("bessel_comparisons", [])) != 4
            or len(global_data.get("negative_identities", [])) != 3
            or any("full_over_local" not in row or "segment_cancellation" not in row
                   for row in global_rows)
            or any("phase_corrected_relative_difference" not in row
                   for row in global_data.get("bessel_comparisons", []))):
        errors.append("invalid full-contour and auxiliary Bessel diagnostics")
    supporting_records.append({"file": global_rel,
                               "kind": "finite-parameter complete contour diagnostics and auxiliary Bessel identities; not a uniform full Barnes limit",
                               "checks": len(global_checks)})

    bessel_rel = "results/stage4/recoil-bessel-limit/diagnostics.json"
    bessel = read_json(ROOT/bessel_rel)
    bessel_checks = bessel.get("checks", [])
    if (bessel.get("status") != "passed" or len(bessel_checks) != 69
            or any(row.get("passed") is not True for row in bessel_checks)
            or len(bessel.get("tail_bound_rows", [])) != 36
            or len(bessel.get("compact_rows", [])) != 6
            or len(bessel.get("integral_rows", [])) != 6
            or not bessel.get("original_contour_comparison")
            or any("scaled_nu2_absolute_error" not in row or "cancellation" not in row
                   for row in bessel.get("integral_rows", []))):
        errors.append("invalid fixed-K full Barnes-to-Bessel diagnostics")
    supporting_records.append({"file": bessel_rel,
                               "kind": "recurrence/tail samples and fixed-K full contour correction checks accompanying a self-reviewed proof; not joint growing-K or interval certification",
                               "checks": len(bessel_checks)})

    cancellation_rel = "results/stage4/recoil-bessel-cancellation/diagnostics.json"
    cancellation = read_json(ROOT/cancellation_rel)
    cancellation_checks = cancellation.get("checks", [])
    if (cancellation.get("status") != "passed" or len(cancellation_checks) != 65
            or any(row.get("passed") is not True for row in cancellation_checks)
            or len(cancellation.get("zero_rows", [])) != 6
            or len(cancellation.get("zero_integral_rows", [])) != 6
            or len(cancellation.get("compact_rows", [])) != 6
            or len(cancellation.get("fixed_second_rows", [])) != 6
            or len(cancellation.get("joint_coefficient_comparisons", [])) != 4
            or len(cancellation.get("joint_tail_rows", [])) != 18
            or any("relative_prediction_error" not in row or "cancellation" not in row
                   for row in cancellation.get("zero_integral_rows", []))):
        errors.append("invalid cancellation/zero-lifting/joint-tail diagnostics")
    supporting_records.append({"file": cancellation_rel,
                               "kind": "fixed-parameter coefficient and zero-lifting proof diagnostics plus joint positive-tail samples; not a complete joint asymptotic or physical critical point",
                               "checks": len(cancellation_checks)})

    endpoints_rel = "results/stage4/recoil-bessel-endpoints/diagnostics.json"
    endpoints = read_json(ROOT/endpoints_rel)
    endpoint_checks = endpoints.get("checks", [])
    if (endpoints.get("status") != "passed" or len(endpoint_checks) != 300
            or len({row.get("name") for row in endpoint_checks}) != 300
            or any(row.get("passed") is not True for row in endpoint_checks)
            or len(endpoints.get("vertical_rows", [])) != 120
            or len(endpoints.get("base_rows", [])) != 12
            or len(endpoints.get("positive_base_samples", [])) != 36
            or len(endpoints.get("local_multiplier_samples", [])) != 36
            or len(endpoints.get("central_rows", [])) != 5):
        errors.append("invalid joint endpoint / central absolute-norm diagnostic schema")
    try:
        with localcontext() as ctx:
            ctx.prec = 100
            for row in endpoint_checks:
                err, tol = Decimal(row["error"]), Decimal(row["tolerance"])
                if not err.is_finite() or err < 0 or err > tol:
                    errors.append("endpoint check numeric mismatch: " + row["name"])
            for row in endpoints["vertical_rows"]:
                actual = Decimal(row["log_ratio"])
                if (abs(actual-Decimal(row["product_log_ratio"])) > Decimal("1e-55")
                        or actual-Decimal(row["log_upper"]) > Decimal("1e-55")):
                    errors.append("endpoint vertical raw-row identity/bound mismatch")
            for row in endpoints["central_rows"]:
                if (row["dps"] != 75 or row["order"] != 32
                        or Decimal(row["U_absolute_over_E"]) <= 0
                        or Decimal(row["Q_absolute_over_E"]) <= 0
                        or len(row["quadrature_changes"]) != 2
                        or any(Decimal(v) > Decimal("1e-16") for v in row["quadrature_changes"])):
                    errors.append("endpoint central absolute quadrature mismatch")
            if (len(endpoints["central_rows"][-1]["precision_changes"]) != 2
                    or any(Decimal(v) > Decimal("1e-55")
                           for v in endpoints["central_rows"][-1]["precision_changes"])):
                errors.append("endpoint central precision mismatch")
    except (KeyError, TypeError, ValueError, ArithmeticError, IndexError) as exc:
        errors.append("malformed endpoint numeric records: " + str(exc))
    supporting_records.append({"file": endpoints_rel,
                               "kind": "joint endpoint identity/bound diagnostics and central absolute norms; not full signed-amplitude or interval certification",
                               "checks": len(endpoint_checks)})

    joint_rel = "results/stage4/recoil-bessel-joint/diagnostics.json"
    joint = read_json(ROOT/joint_rel)
    joint_checks = joint.get("checks", [])
    if (joint.get("status") != "passed" or len(joint_checks) != 64
            or len({row.get("name") for row in joint_checks}) != 64
            or any(row.get("passed") is not True for row in joint_checks)
            or len(joint.get("finite_contours", [])) != 3
            or len(joint.get("positive_moments", [])) != 3
            or len(joint.get("vertical_moments", [])) != 3
            or len(joint.get("legacy_rows", [])) != 6
            or len(joint.get("joint_rows", [])) != 12):
        errors.append("invalid phase-preserving joint diagnostic schema")
    try:
        with localcontext() as ctx:
            ctx.prec = 100
            def decimal_complex(s):
                vals = re.findall(r'[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?', s.replace(" ", ""))
                if len(vals) != 2:
                    raise ValueError("invalid complex record")
                return tuple(map(Decimal, vals))
            for row in joint_checks:
                err, tol = Decimal(row["error"]), Decimal(row["tolerance"])
                if not err.is_finite() or err < 0 or err > tol:
                    errors.append("joint check numeric mismatch: " + row["name"])
            for row in joint["finite_contours"]:
                z, main, boundary = map(decimal_complex,
                    (row["original"], row["branch_main"], row["boundary_correction"]))
                residual = sum((z[j]-main[j]-boundary[j])**2 for j in [0, 1]).sqrt()
                if residual/Decimal(row["envelope"]) > Decimal("1e-35"):
                    errors.append("joint finite contour signed identity mismatch")
            for row in joint["joint_rows"]:
                z, pred = map(decimal_complex, (row["value"], row["leading"]))
                actual = sum((z[j]-pred[j])**2 for j in [0, 1]).sqrt()/Decimal(row["envelope"])
                if abs(actual-Decimal(row["envelope_error"])) > Decimal("1e-60"):
                    errors.append("joint raw complex / envelope error mismatch")
                if Decimal(row["quadrature_envelope_change"]) > Decimal("1e-18"):
                    errors.append("joint quadrature record mismatch")
            for row in joint["legacy_rows"]:
                if Decimal(row["legacy_envelope_change"]) > Decimal("1e-13"):
                    errors.append("joint frozen full-contour comparison mismatch")
        initial = read_json(ROOT/"results/stage4/recoil-bessel-joint/initial-contour-grid.json")
        if (initial.get("status") != "failed" or len(initial.get("checks", [])) != 55
                or sum(row.get("passed") is not True for row in initial["checks"]) != 2
                or initial["hashes"]["code/recoil_bessel_joint.py"] != digest(
                    ROOT/"results/stage4/recoil-bessel-joint/initial-contour-grid.py")):
            errors.append("joint initial failed-grid provenance mismatch")
    except (KeyError, TypeError, ValueError, ArithmeticError, IndexError) as exc:
        errors.append("malformed phase-preserving joint records: " + str(exc))
    supporting_records.append({"file": joint_rel,
                               "kind": "exact contour-split and moment checks supporting a self-reviewed joint leading estimate; no numerical threshold or finite-G certification",
                               "checks": len(joint_checks), "retained_initial_failed_checks": 2})

    transverse_rel = "results/stage4/recoil-transverse/diagnostics.json"
    transverse = read_json(ROOT/transverse_rel)
    consequences = read_json(ROOT/"results/stage4/recoil-transverse/consequences.json")
    try:
        for rec, count in [(transverse, 67), (consequences, 16)]:
            checks = rec["checks"]
            if (rec["status"] != "passed" or len(checks) != count
                    or len({r["name"] for r in checks}) != count):
                errors.append("invalid transverse diagnostic schema")
            for row in checks:
                err, tol = Decimal(row["error"]), Decimal(row["tolerance"])
                if row["passed"] is not True or not err.is_finite() or err < 0 or err > tol:
                    errors.append("transverse check numeric mismatch: " + row["name"])
        if (len(transverse["finite_contours"]) != 2 or len(transverse["joint_rows"]) != 12
                or len(transverse["transverse_rows"]) != 9 or len(transverse["tail_rows"]) != 16
                or len(consequences["multiplier_rows"]) != 16 or len(consequences["growth_rows"]) != 3):
            errors.append("incomplete transverse raw records")
        with localcontext() as ctx:
            ctx.prec = 110
            for row in transverse["finite_contours"]:
                z, main, boundary = map(decimal_complex, (row["original"], row["main"], row["correction"]))
                residual = sum((z[j]-main[j]-boundary[j])**2 for j in [0,1]).sqrt()
                if residual/Decimal(row["envelope"]) > Decimal("1e-32"):
                    errors.append("transverse signed contour identity mismatch")
            for row in transverse["joint_rows"]:
                z, pred = map(decimal_complex, (row["value"], row["leading"]))
                actual = sum((z[j]-pred[j])**2 for j in [0,1]).sqrt()/Decimal(row["envelope"])
                K,nu,k,alpha = map(Decimal,(row["K"],row["nu"],row["k"],row["alpha"]))
                if (abs(actual-Decimal(row["envelope_error"])) > Decimal("1e-60")
                        or abs(Decimal(row["rho"])-K**3/nu**2) > Decimal("1e-60")
                        or not 0 < k <= K or Decimal(row["theta"]) >= 1
                        or Decimal(row["quadrature_change"]) > Decimal("1e-17")):
                    errors.append("transverse raw amplitude/scale mismatch")
            for row in consequences["multiplier_rows"]:
                z=decimal_complex(row["log_residual"])
                actual=sum(x*x for x in z).sqrt()/Decimal(row["rho"])
                if abs(actual-Decimal(row["residual_over_rho"])) > Decimal("1e-65"):
                    errors.append("transverse multiplier residual mismatch")
            for row in consequences["growth_rows"]:
                nr=next(r for r in transverse["joint_rows"] if r["family"]=="fractional"
                        and Decimal(r["K"])==Decimal(row["K"]) and Decimal(r["alpha"])==2*Decimal(row["K"]))
                sr=next(r for r in joint["joint_rows"] if Decimal(r["K"])==Decimal(row["K"])
                        and Decimal(r["alpha"])==2*Decimal(row["K"]))
                zn,zs=map(decimal_complex,(nr["value"],sr["value"]))
                ratio=(sum(x*x for x in zn)/sum(x*x for x in zs)).sqrt()*Decimal(row["C_modulus_ratio"])
                if abs(ratio/Decimal(row["raw_internal_modulus_ratio"])-1)>Decimal("1e-65"):
                    errors.append("transverse growth does not match saved complex amplitudes")
    except (KeyError, TypeError, ValueError, ArithmeticError, IndexError, StopIteration) as exc:
        errors.append("malformed transverse records: " + str(exc))
    supporting_records.append({"file": transverse_rel, "checks": 83,
        "kind": "strict-gap transverse internal-graph estimate and branch-growth diagnostics; no physical outer-integral or interval certification"})

    pairing = read_json(ROOT/(PAIRING_DIR + "diagnostics.json"))
    try:
        checks = pairing["checks"]
        if (pairing["status"] != "passed" or len(checks) != 37
                or len({r["name"] for r in checks}) != 37):
            errors.append("invalid pairing diagnostic schema")
        for row in checks:
            if (row["passed"] is not True or not Decimal(row["error"]).is_finite()
                    or not 0 <= Decimal(row["error"]) <= Decimal(row["tolerance"])):
                errors.append("pairing numeric check mismatch: " + row["name"])
        for key, count in [("vertex_rows", 3), ("bessel_rows", 4), ("symmetric_rows", 12),
                           ("transverse_rows", 3), ("window_integrals", 6)]:
            if len(pairing[key]) != count:
                errors.append("incomplete pairing data: " + key)
        initial = read_json(ROOT/(PAIRING_DIR + "initial-window-grid.json"))
        if (initial["status"] != "failed" or len(initial["checks"]) != 37
                or sum(r["passed"] is False for r in initial["checks"]) != 4
                or initial["hashes"]["code/recoil_pairing.py"] != digest(ROOT/(PAIRING_DIR + "initial-window-grid.py"))):
            errors.append("invalid preserved pairing failure or original code hash")
        checked_hashes += 1
        with localcontext() as ctx:
            ctx.prec = 110
            for row in pairing["bessel_rows"]:
                actual = sum(Decimal(row[k]) for k in ["K_square", "I_square", "interference"])
                if abs(actual/Decimal(row["exact_square"])-1) > Decimal("1e-70"):
                    errors.append("paired Bessel interference arithmetic mismatch")
            for row in pairing["symmetric_rows"]:
                raw = next(r for r in joint["joint_rows"] if r["K"] == row["K"] and r["alpha"] == row["alpha"])
                ref = next(r for r in joint["joint_rows"] if r["K"] == row["K"] and Decimal(r["alpha"]) == 0)
                factor = (sum(x*x for x in decimal_complex(raw["value"])) /
                          sum(x*x for x in decimal_complex(ref["value"])))
                z, vertex = map(decimal_complex, (row["complete_graph_ratio"], row["vertex_ratio"]))
                if max(abs(z[j]-vertex[j]*factor) for j in [0,1])/max(abs(x) for x in z) > Decimal("1e-65"):
                    errors.append("paired symmetric graph differs from saved complex amplitudes")
            for row in pairing["transverse_rows"]:
                old = next(r for r in consequences["growth_rows"] if r["K"] == row["K"])
                factor = Decimal(old["raw_internal_modulus_ratio"])**2
                z, vertex = map(decimal_complex, (row["complete_transverse_ratio"], row["analytic_vertex_ratio"]))
                if max(abs(z[j]-vertex[j]*factor) for j in [0,1])/max(abs(x) for x in z) > Decimal("1e-65"):
                    errors.append("paired transverse graph arithmetic mismatch")
            for row in pairing["window_integrals"]:
                actual = Decimal(row["scaled_integral"])/Decimal(row["endpoint_scaled_prediction"])
                if (abs(actual-Decimal(row["endpoint_ratio"])) > Decimal("1e-65")
                        or Decimal(row["relative_refinement"]) > Decimal("1e-20")):
                    errors.append("paired auxiliary-window mismatch")
    except (KeyError, TypeError, ValueError, ArithmeticError, IndexError, StopIteration) as exc:
        errors.append("malformed pairing records: " + str(exc))
    supporting_records.append({"file": PAIRING_DIR + "diagnostics.json", "checks": 37,
        "kind": "Schwarz pairing and topology graph; reused finite-branch data and auxiliary-window quadrature, not physical OTOC certification"})

    poles = read_json(ROOT/(POLES_DIR + "diagnostics.json"))
    try:
        checks = poles["checks"]
        if (poles["status"] != "passed" or len(checks) != 74
                or len({r["name"] for r in checks}) != 74):
            errors.append("invalid recoil-pole diagnostic schema")
        for row in checks:
            if (row["passed"] is not True or not Decimal(row["error"]).is_finite()
                    or not 0 <= Decimal(row["error"]) <= Decimal(row["tolerance"])):
                errors.append("recoil-pole numeric mismatch: " + row["name"])
        for key,count in [("identity_rows",3),("pinch_rows",6),("direct_c_rows",3),
                          ("shift_rows",6),("direct_c_poles",2)]:
            if len(poles[key]) != count:
                errors.append("incomplete recoil-pole raw data: " + key)
        def complex_relative_difference(z, w):
            return (sum((z[j]-w[j])**2 for j in [0,1])/sum(x*x for x in w)).sqrt()
        with localcontext() as ctx:
            ctx.prec = 95
            for row in poles["identity_rows"]:
                for first,second in [("left_H","right_H"),("direct_graph","regular_graph")]:
                    if complex_relative_difference(decimal_complex(row[first]),decimal_complex(row[second])) > Decimal("1e-27"):
                        errors.append("pole representation mismatch")
            for row in poles["pinch_rows"]:
                if row["multiplicity"] != {"simple":1,"double":2,"quadruple":4}[row["family"]] or len(row["approaches"]) != 3:
                    errors.append("pole multiplicity/approach schema mismatch")
                z0=decimal_complex(row["graph0"])
                if complex_relative_difference(z0,decimal_complex(row["graph0_high"])) > Decimal("1e-28"):
                    errors.append("pole value precision mismatch")
                for approach in row["approaches"]:
                    b=decimal_complex(approach["B"])
                    scaled=tuple(x*Decimal(approach["epsilon"])**row["multiplicity"] for x in b)
                    if complex_relative_difference(scaled,decimal_complex(approach["scaled_pole_coefficient"])) > Decimal("1e-35"):
                        errors.append("scaled pole coefficient arithmetic mismatch")
                    for actual,reference,key in [(scaled,decimal_complex(row["leading_B"]),"relative_pole_limit"),
                            (decimal_complex(approach["graph"]),z0,"relative_graph_limit")]:
                        if abs(complex_relative_difference(actual,reference)-Decimal(approach[key])) > Decimal("1e-33"):
                            errors.append("pole approach residual arithmetic mismatch")
            for row in poles["direct_c_rows"]:
                c,cs=map(decimal_complex,(row["chiral_C"],row["sharp_C"]))
                pair=(c[0]*cs[0]-c[1]*cs[1],c[0]*cs[1]+c[1]*cs[0])
                paired_target = (decimal_complex(row["paired_C"]) if "j" in row["paired_C"]
                                 else (Decimal(row["paired_C"]), Decimal(0)))
                if complex_relative_difference(pair,paired_target) > Decimal("1e-32"):
                    errors.append("direct C Schwarz product mismatch")
            for row in poles["direct_c_poles"]:
                for approach in row["approaches"]:
                    scaled=tuple(x*Decimal(approach["epsilon"]) for x in decimal_complex(approach["value"]))
                    if abs(complex_relative_difference(scaled,decimal_complex(row["residue"]))-Decimal(approach["relative"])) > Decimal("1e-32"):
                        errors.append("direct C residue arithmetic mismatch")
    except (KeyError, TypeError, ValueError, ArithmeticError, IndexError) as exc:
        errors.append("malformed recoil-pole records: " + str(exc))
    supporting_records.append({"file": POLES_DIR + "diagnostics.json", "checks": 74,
        "kind": "regularized recoil graph and bare direct-C0 factors; includes merged pinches, not physical sewing or finite-G certification"})

    spectral = read_json(ROOT/(SPECTRAL_DIR + "diagnostics.json"))
    try:
        checks = spectral["checks"]
        if (spectral["status"] != "passed" or len(checks) != 71
                or len({r["name"] for r in checks}) != 71):
            errors.append("invalid spectral diagnostic schema")
        for row in checks:
            if (row["passed"] is not True or not Decimal(row["error"]).is_finite()
                    or not 0 <= Decimal(row["error"]) <= Decimal(row["tolerance"])):
                errors.append("spectral numeric mismatch: " + row["name"])
        for key,count in [("shape_rows",9),("protocol_rows",6),("fourier_rows",5),
                          ("freedom_rows",9),("trace_rows",6),("time_rows",10)]:
            if len(spectral[key]) != count:
                errors.append("incomplete spectral records: " + key)
        initial = read_json(ROOT/(SPECTRAL_DIR + "initial-diagnostics.json"))
        if (initial["status"] != "passed" or len(initial["checks"]) != 59
                or initial["hashes"]["code/recoil_spectral.py"] != digest(ROOT/(SPECTRAL_DIR + "initial-recoil-spectral.py"))):
            errors.append("spectral initial-record provenance mismatch")
        checked_hashes += 1
        with localcontext() as ctx:
            ctx.prec = 90
            for row in spectral["shape_rows"]:
                if abs(Decimal(row["assembled"])/Decimal(row["compact"])-1) > Decimal("1e-45"):
                    errors.append("assembled two-point spectrum mismatch")
            for row in spectral["protocol_rows"]:
                numerator=Decimal(row["density_initial"])*Decimal(row["forward"])
                denominator=Decimal(row["density_final"])*Decimal(row["reverse_other_state"])
                if abs(numerator/denominator-1) > Decimal("1e-45"):
                    errors.append("reversed-transition density mismatch")
            for row in spectral["freedom_rows"]:
                actual=4*Decimal(row["M"])**2*(Decimal(row["ratio"])-1)/Decimal(row["omega"])**2
                if abs(actual-Decimal(row["scaled_correction"])) > Decimal("1e-45"):
                    errors.append("normalization-family arithmetic mismatch")
            for row in spectral["time_rows"]:
                val,low,coef=map(Decimal,(row["value"],row["low"],row["tail_coefficient"]))
                trunc=Decimal(row["baseline_tail_bound"])+Decimal(row["series_time_bound"])
                if (abs(val*Decimal(row["t"])**2/coef-Decimal(row["tail_ratio"])) > Decimal("1e-45")
                        or abs(abs(low/val-1)-Decimal(row["relative_refinement"])) > Decimal("1e-42")
                        or trunc/abs(val) > Decimal("1e-24")
                        or abs(trunc/abs(val)-Decimal(row["relative_analytic_truncation_bound"])) > Decimal("1e-65")):
                    errors.append("trace time/tail arithmetic mismatch")
                if row["extended_value"] is not None and abs(Decimal(row["extended_value"])-val) > Decimal(row["baseline_tail_bound"])*Decimal("1.01"):
                    errors.append("trace cutoff extension mismatch")
    except (KeyError, TypeError, ValueError, ArithmeticError, IndexError) as exc:
        errors.append("malformed spectral records: " + str(exc))
    supporting_records.append({"file": SPECTRAL_DIR + "diagnostics.json", "checks": 71,
        "kind": "two-point calibration and declared n=1 spectral trace model; insertion freedom retained, no physical four-point or interval certification"})

    identity = read_json(ROOT/(IDENTITY_DIR + "diagnostics.json"))
    try:
        def identity_number(value):
            if 'j' not in value:
                return Decimal(value), Decimal(0)
            return decimal_complex(value)
        checks = identity["checks"]
        if (identity["status"] != "passed" or len(checks) != 106
                or len({r["name"] for r in checks}) != 106):
            errors.append("invalid identity-endpoint diagnostic schema")
        with localcontext() as ctx:
            ctx.prec = 90
            for row in checks:
                error, tol = map(Decimal, (row["absolute_error"], row["tolerance"]))
                z, w = map(identity_number, (row["observed"], row["expected"]))
                actual = sum((z[j]-w[j])**2 for j in [0,1]).sqrt()
                if (row["passed"] is not True or not error.is_finite()
                        or not 0 <= error <= tol or actual > tol):
                    errors.append("identity-endpoint numeric mismatch: " + row["name"])
            for row in identity["rows"]:
                if row["kind"] == "pole_crossing":
                    residual = (Decimal(row["real_axis"])+Decimal(row["residue_correction"])
                                -Decimal(row["analytically_continued"]))
                    if abs(residual) > Decimal("1e-30"):
                        errors.append("identity-endpoint saved crossing mismatch")
                else:
                    e = Decimal(row["epsilon"])
                    if abs(Decimal(row["variance"])-(2*e-e*e)/3) > Decimal("1e-30"):
                        errors.append("identity-endpoint saved moment mismatch")
            for row in identity["weights"]:
                hp, hm = map(identity_number, (row["h_plus"], row["h_minus"]))
                if abs(hp[0]+hm[0]-1)+abs(hp[1]+hm[1]) > Decimal("1e-30"):
                    errors.append("identity-endpoint paired weights mismatch")
        for key,count in [("rows",13),("finite_mass",12),("weights",9),("scaling",9)]:
            if len(identity[key]) != count:
                errors.append("incomplete identity-endpoint records: " + key)
    except (KeyError, TypeError, ValueError, ArithmeticError) as exc:
        errors.append("malformed identity-endpoint records: " + str(exc))
    supporting_records.append({"file": IDENTITY_DIR + "diagnostics.json", "checks": 106,
        "kind": "reduced spectral identity limit and label algebra, including repeated vacuum checks; not physical sewing or interval certification"})

    racah = read_json(ROOT/(RACAH_DIR + "diagnostics.json"))
    try:
        def racah_number(value):
            value = str(value)
            return decimal_complex(value) if 'j' in value else (Decimal(value), Decimal(0))
        def racah_mul(z,w):
            return z[0]*w[0]-z[1]*w[1], z[0]*w[1]+z[1]*w[0]
        def racah_div(z,w):
            den = w[0]**2+w[1]**2
            return (z[0]*w[0]+z[1]*w[1])/den, (z[1]*w[0]-z[0]*w[1])/den
        def racah_norm(z):
            return (z[0]**2+z[1]**2).sqrt()
        def racah_relative(z,w):
            return racah_norm((z[0]-w[0],z[1]-w[1]))/max(racah_norm(w),Decimal("1e-80"))
        def racah_rectangle(values):
            z = list(map(racah_number,values))
            return racah_div(racah_mul(z[0],z[3]),racah_mul(z[1],z[2]))
        checks = racah["checks"]
        if (racah["status"] != "passed" or len(checks) != 44
                or len({r["name"] for r in checks}) != 44):
            errors.append("invalid global Racah diagnostic schema")
        if (racah["external"] != [".4",".7",".93",".6"]
                or racah["s_values"] != [".65",".85"]
                or racah["t_values"] != [".57",".81"]):
            errors.append("global Racah candidate rectangle changed")
        with localcontext() as ctx:
            ctx.prec = 90
            for row in checks:
                err,tol = map(Decimal,(row["relative_error"],row["tolerance"]))
                actual = racah_relative(racah_number(row["observed"]),racah_number(row["expected"]))
                if (row["passed"] is not True or not err.is_finite()
                        or not 0 <= err < tol or actual >= tol):
                    errors.append("global Racah check mismatch: "+row["name"])
            if [d["dps"] for d in racah["datasets"]] != [40,60]:
                errors.append("global Racah precision datasets incomplete")
            for dataset in racah["datasets"]:
                rows = dataset["rows"]
                if len(rows) != 4:
                    raise ValueError("rectangle needs four rows")
                for key in ["J","H","Z"]:
                    values = [r["sixj"]["value"] if key=="J" else r[key] for r in rows]
                    actual = racah_rectangle(values)
                    if racah_relative(actual,racah_number(dataset["ratios"][key])) > Decimal("1e-35"):
                        errors.append("global Racah saved rectangle mismatch: "+key)
                for row in rows:
                    if racah_relative(racah_number(row["H"]),racah_number(row["H_left"])) > Decimal("1e-30"):
                        errors.append("global Racah recoil representation mismatch")
                hr = racah_number(dataset["ratios"]["H"])
                if racah_relative(racah_number(dataset["ratios"]["Z"]),
                                  (hr[0]**2+hr[1]**2,Decimal(0))) > Decimal("1e-35"):
                    errors.append("global Racah Z/H rectangle mismatch")
            high = racah["datasets"][1]
            for key in ["H","Z"]:
                quotient = racah_div(racah_number(high["ratios"][key]),racah_number(high["ratios"]["J"]))
                actual = racah_norm((quotient[0]-1,quotient[1]))
                if abs(actual-Decimal(racah["gaps"][key])) > Decimal("1e-50"):
                    errors.append("global Racah raw gap mismatch: "+key)
                if not Decimal(".04") < actual < Decimal(".06"):
                    errors.append("global Racah exclusion gap not reproduced")
            if racah_relative(racah_number(racah["direct"]["value"]),
                              racah_number(high["rows"][0]["sixj"]["value"])) > Decimal("1e-28"):
                errors.append("global Racah direct Euler mismatch")
            for row in racah["shadow_rows"]:
                actual = racah_rectangle(row["values"])
                if (racah_relative(actual,racah_number(row["ratio"])) > Decimal("1e-50")
                        or racah_relative(actual,racah_number(high["ratios"]["J"])) > Decimal("1e-45")):
                    errors.append("global Racah raw shadow rectangle mismatch")
        if len(racah["shadow_rows"]) != 3 or len(racah["shadow_norm"]) != 3:
            errors.append("global Racah shadow records incomplete")
    except (KeyError,TypeError,ValueError,ArithmeticError) as exc:
        errors.append("malformed global Racah records: "+str(exc))
    supporting_records.append({"file":RACAH_DIR+"diagnostics.json","checks":44,
        "kind":"specified scalar dictionary comparison; raw rectangles and gaps independently reconstructed with Decimal; not an interval certificate or full physical recoupling"})

    ds3 = read_json(ROOT/(DS3_DIR + "diagnostics.json"))
    ds3_clp = read_json(ROOT/(DS3_DIR + "clp-comparison.json"))
    try:
        checks = ds3["checks"]
        if (ds3["status"] != "passed" or len(checks) != 249
                or len({r["name"] for r in checks}) != 249):
            raise ValueError("expected 249 unique passed diagnostics")
        for key, count in [("modes",36),("vertices",24),("boundary",29),("ward",27),
                           ("direct",4),("shadow",8),("precision",9)]:
            if len(ds3[key]) != count:
                raise ValueError("incomplete raw rows: " + key)
        with localcontext() as ctx:
            ctx.prec = 90
            for row in checks:
                z,w = map(racah_number,(row["observed"],row["expected"]))
                scale = Decimal(str(row["scale"]))
                if not scale.is_finite() or scale <= 0:
                    raise ValueError("nonpositive diagnostic scale")
                actual = racah_norm((z[0]-w[0],z[1]-w[1]))/scale
                err,tol = map(Decimal,(row["scaled_error"],row["tolerance"]))
                if (row["passed"] is not True or not err.is_finite()
                        or not 0 <= err < tol or actual >= tol
                        or abs(actual-err)*scale > Decimal("1e-58")):
                    errors.append("dS3 raw diagnostic mismatch: "+row["name"])
            for row in ds3["vertices"]:
                c = racah_number(row["coefficient"])
                norm = c[0]**2+c[1]**2
                if (racah_relative(racah_number(row["value"]),racah_number(row["closed"])) > Decimal("1e-45")
                        or abs(norm/Decimal(row["norm"])-1) > Decimal("1e-50")
                        or abs(norm/Decimal(row["closed_norm"])-1) > Decimal("1e-45")):
                    errors.append("dS3 raw spherical coefficient mismatch")
            for row in ds3["boundary"]:
                if "bulk" in row:
                    predicted = racah_mul(racah_number(row["coefficient"]),racah_number(row["boundary"]))
                    if racah_relative(racah_number(row["bulk"]),predicted) > Decimal("1e-42"):
                        errors.append("dS3 raw bulk/boundary mismatch")
            for row in ds3["ward"]:
                terms = list(map(racah_number,row["terms"]))
                total = tuple(sum(t[i] for t in terms) for i in [0,1])
                scale = sum(racah_norm(t) for t in terms)
                saved = racah_number(row["residual"])
                if (racah_norm(total)/scale >= Decimal("1e-40")
                        or racah_norm((total[0]-saved[0],total[1]-saved[1])) > Decimal("1e-58")
                        or abs(scale/Decimal(row["scale"])-1) > Decimal("1e-55")):
                    errors.append("dS3 raw boost Ward mismatch")
            for row in ds3["direct"]:
                beta,refined,direct = map(racah_number,(row["beta"],row["refined"],row["direct"]))
                tail = Decimal(row["analytic_tail_bound"])
                actual = racah_relative(refined,beta)
                if (row["T"] != 96 or row["step"] != [3,2] or tail < 0
                        or not tail.is_finite() or actual >= Decimal("1e-34")
                        or racah_relative(direct,refined) >= Decimal("1e-35")
                        or abs(actual-Decimal(row["relative_difference"])) > Decimal("1e-55")
                        or abs(tail/racah_norm(beta)-Decimal(row["tail_relative"])) > Decimal("1e-55")):
                    errors.append("dS3 raw quadrature/tail mismatch")
            rows = ds3_clp["rows"]
            if (ds3_clp["status"] != "passed" or len(rows) != 24
                    or len({(tuple(r["nus"]),tuple(r["conjugates"])) for r in rows}) != 24):
                raise ValueError("expected 24 distinct CLP rows")
            pi = Decimal("3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628")
            for row in rows:
                c = racah_number(row["coefficient"])
                norm = c[0]**2+c[1]**2
                ratio = norm/Decimal(row["direct_C_pair"])
                nus = list(map(Decimal,row["nus"]))
                cc = row["conjugates"]
                predicted = 1/(64*pi**4)
                for nu in nus:
                    e = (pi*nu).exp()
                    predicted *= nu/e/((e-1/e)/2)
                if sum(cc) in [1,2]:
                    r = cc.index(True) if sum(cc)==1 else cc.index(False)
                    e = (pi*nus[r]).exp()
                    predicted *= ((e+1/e)/2)**2
                actual = abs(ratio/predicted-1)
                if (row["passed"] is not True or actual >= Decimal(row["tolerance"])
                        or abs(norm/Decimal(row["norm"])-1) > Decimal("1e-55")
                        or abs(ratio/Decimal(row["ratio"])-1) > Decimal("1e-55")
                        or abs(predicted/Decimal(row["expected"])-1) > Decimal("1e-55")
                        or abs(actual-Decimal(row["relative_error"])) > Decimal("1e-55")):
                    errors.append("dS3 raw CLP leg factor mismatch")
    except (KeyError,TypeError,ValueError,ArithmeticError) as exc:
        errors.append("malformed dS3 local vertex records: "+str(exc))
    supporting_records.append({"file":DS3_DIR+"diagnostics.json","checks":249,
        "kind":"fixed-mass KG local vertex; saved errors and raw identities reconstructed with Decimal; not observer spectral normalization or interval quadrature"})
    supporting_records.append({"file":DS3_DIR+"clp-comparison.json","checks":24,
        "kind":"separable local leg factors reconstructed with Decimal; not physical two-point normalization"})

    mass = read_json(ROOT/(MASS_DIR + "diagnostics.json"))
    mass_initial = read_json(ROOT/(MASS_DIR + "initial-diagnostics.json"))
    try:
        initial_failed = [r for r in mass_initial["checks"] if not r["passed"]]
        if (mass_initial["status"] != "failed" or len(mass_initial["checks"]) != 146
                or [r["name"] for r in initial_failed] != ["packet_norm_20","packet_antipodal_20"]
                or mass["historical_failed_targets"] != initial_failed):
            raise ValueError("short-box failure history changed")
        checks = mass["checks"]
        if (mass["status"] != "passed" or len(checks) != 150
                or len({r["name"] for r in checks}) != 150):
            raise ValueError("expected 150 unique diagnostics")
        for key,count in [("gram",4),("modes",24),("bounds",10),("parseval",4),("packets",3),("legs",4)]:
            if len(mass[key]) != count:
                raise ValueError("incomplete mass-spectrum raw rows: "+key)
        with localcontext() as ctx:
            ctx.prec = 90
            pi = Decimal("3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628")
            for row in checks:
                z,w=map(racah_number,(row["observed"],row["expected"]))
                scale=Decimal(row["scale"])
                actual=racah_norm((z[0]-w[0],z[1]-w[1]))/scale
                err,tol=map(Decimal,(row["scaled_error"],row["tolerance"]))
                if (scale <= 0 or row["passed"] is not True or not 0<=err<tol
                        or actual>=tol or abs(actual-err)>Decimal("1e-42")):
                    errors.append("mass-spectrum check mismatch: "+row["name"])
            for row in mass["gram"]:
                nu=Decimal(row["nu"]);e=(pi*nu).exp();q=1/e
                a=pi/nu*(e+q)/(e-q);b=-2*pi/nu/(e-q)
                if (abs(Decimal(row["a"])/a-1)>Decimal("1e-38")
                        or abs(Decimal(row["b"])/b-1)>Decimal("1e-38")
                        or abs(Decimal(row["determinant"])/(a*a-b*b)-1)>Decimal("1e-38")
                        or abs(Decimal(row["mu"])*a-1)>Decimal("1e-38")
                        or abs(Decimal(row["overlap"])+2/(e+q))>Decimal("1e-38")):
                    errors.append("mass-spectrum raw Gram mismatch")
                expected=[[1/(1+q*q),-q/(1+q*q)],[-q/(1+q*q),q*q/(1+q*q)]]
                if any(abs(Decimal(row["projector"][i][j])-expected[i][j])>Decimal("1e-38")
                       for i in [0,1] for j in [0,1]):
                    errors.append("mass-spectrum raw projector mismatch")
            for row in mass["modes"]:
                nu=Decimal(row["nu"]);q=(-pi*nu).exp()
                N=(2*nu*(1-q*q)).sqrt()
                factor=(2*pi).sqrt()/N
                em,ep=map(racah_number,(row["e_minus"],row["e_plus"]))
                bd=tuple(factor*(em[i]-q*ep[i]) for i in [0,1])
                anti=tuple(factor*(ep[i]-q*em[i]) for i in [0,1])
                if (racah_relative(racah_number(row["bd"]),bd)>Decimal("1e-38")
                        or racah_relative(racah_number(row["antipode"]),anti)>Decimal("1e-38")):
                    errors.append("mass-spectrum raw channel mismatch")
            for row in mass["bounds"]:
                if (abs(Decimal(row["norm"])-1)>Decimal("1e-38")
                        or racah_norm(racah_number(row["ode_residual"]))>Decimal("1e-38")
                        or racah_norm(racah_number(row["overlap"]))>Decimal("1e-30")):
                    errors.append("mass-spectrum raw bound-state mismatch")
            for row in mass["parseval"]:
                discrete=sum(racah_norm(racah_number(v))**2 for v in row["bound_overlaps"])
                continuous=Decimal(row["continuum"]);total=continuous+discrete
                if (continuous<0 or abs(discrete-Decimal(row["bound_weight"]))>Decimal("1e-38")
                        or abs(total-Decimal(row["total"]))>Decimal("1e-38")
                        or abs(total-Decimal(row["expected"]))>Decimal("1e-32")
                        or racah_norm(tuple(a-b for a,b in zip(racah_number(row["transform_direct"]),
                            racah_number(row["transform_beta"]))))>Decimal("1e-28")):
                    errors.append("mass-spectrum raw Parseval mismatch")
            energy=Decimal(str(mass["packet_tail"]["estimates"][-1]["derivative_energy"]))
            for row in mass["packets"]:
                r=row["refined"];T=Decimal(r["T"]);tail=energy/T**8
                if (tail<=0 or abs(tail/Decimal(str(row["tail_estimate"]))-1)>Decimal("1e-14")
                        or abs(Decimal(str(r["norm"]))-1)>tail+Decimal("1e-12")):
                    errors.append("mass-spectrum raw packet-tail mismatch")
                antipodal=racah_number(r["antipodal_overlap"])
                expected=Decimal(str(r["expected_antipodal"]))
                if (racah_norm((antipodal[0]-expected,antipodal[1]))>tail+Decimal("1e-12")
                        or racah_norm(racah_number(r["orthogonal_overlap"]))>tail+Decimal("1e-12")):
                    errors.append("mass-spectrum raw packet overlap mismatch")
                if T==80 and abs(antipodal[0])<Decimal(".05"):
                    errors.append("mass-spectrum antipodal nonorthogonality lost")
            for row in mass["legs"]:
                ratio=racah_norm(racah_number(row["coefficient"]))**2/racah_norm(racah_number(row["old_coefficient"]))**2
                expected=Decimal(1)
                for s in row["nus"][:2]:
                    nu=Decimal(s);e=(2*pi*nu).exp();expected*=nu/pi*(e-1)/(e+1)
                if (abs(ratio/expected-1)>Decimal("1e-38")
                        or abs(ratio/Decimal(row["norm_ratio"])-1)>Decimal("1e-38")):
                    errors.append("mass-spectrum observer-leg mismatch")
    except (KeyError,TypeError,ValueError,ArithmeticError) as exc:
        errors.append("malformed mass-spectrum records: "+str(exc))
    supporting_records.append({"file":MASS_DIR+"diagnostics.json","checks":150,
        "kind":"spacetime L2/BD projector and packet checks; raw Gram, Parseval and overlaps reconstructed with Decimal; two initial short-box failures retained; not physical group averaging or interval certification"})

    try:
        from audit_ds3_group_pairing import audit as audit_group_pairing
        errors.extend(audit_group_pairing(ROOT, racah_number, racah_norm))
    except (KeyError,TypeError,ValueError,ArithmeticError) as exc:
        errors.append("malformed group-pairing records: "+str(exc))
    supporting_records.append({"file": GROUP_DIR+"diagnostics.json", "checks":81,
        "kind":"G=0 group norms, linear gluing and physical two-point; saved errors and raw algebra independently replayed with Decimal; not full Fock or finite gravity"})
    supporting_records.append({"file": GROUP_DIR+"offdiagonal.json", "checks":16,
        "kind":"independent non-spherical group integrals; saved rank-one products and tails replayed; not interval quadrature"})

    try:
        from audit_ds3_two_matter import audit as audit_two_matter
        errors.extend(audit_two_matter(ROOT, racah_number, racah_norm))
    except (KeyError,TypeError,ValueError,ArithmeticError) as exc:
        errors.append("malformed two-matter records: "+str(exc))
    supporting_records.append({"file": TWO_DIR+"diagnostics.json", "checks":56,
        "kind":"two-species scalar observer channel; DS/group norms, shifted uncrossed proper-time limit and residuals replayed with Decimal; not crossed spin completeness or interval certification"})
    supporting_records.append({"file": TWO_DIR+"angular.json", "checks":10,
        "kind":"non-spherical four-body group integrals; l1 polynomial and tails independently replayed"})

    try:
        from audit_ds3_crossing_dictionary import audit as audit_crossing_dictionary
        errors.extend(audit_crossing_dictionary(ROOT))
    except (KeyError,TypeError,ValueError,ArithmeticError) as exc:
        errors.append("malformed crossing dictionary records: "+str(exc))
    supporting_records.append({"file":CROSS_DICT_DIR+"diagnostics.json","checks":33,
        "kind":"physical G=0 scalar recoupling dictionary; independent Euclidean complex overlap, phase/reflection norms and analytic puncture tail; Decimal replay and source AST check, not interval quadrature or full scalar unitarity"})
    supporting_records.append({"file":CROSS_DICT_DIR+"geometry.json","checks":10,
        "kind":"six-chart coverage and analytic tail primitive; does not certify the hypergeometric quadrature remainder"})

    try:
        from audit_ds3_crossed_packet import audit as audit_crossed_packet
        errors.extend(audit_crossed_packet(ROOT, racah_number, racah_norm))
    except (KeyError,TypeError,ValueError,ArithmeticError) as exc:
        errors.append("malformed crossed packet records: "+str(exc))
    supporting_records.append({"file":PACKET4_DIR+"summary.json", "kernel_points":3442832,
        "kind":"actual regulated compact initial packet, five complete complex 6j grids and independent E/W/rho replay; analytic mass-tail and finite-K nonzero witness; numerical magnitude not total interval certified"})
    supporting_records.append({"file":PACKET4_DIR+"bounds-validation.json", "checks":24,
        "kind":"vertex, missing-interval commutator, positive norms, tail primitive, and independently group-averaged witness"})

    try:
        from audit_ds3_bose import audit as audit_bose
        errors.extend(audit_bose(ROOT))
    except (KeyError, TypeError, ValueError, ArithmeticError) as exc:
        errors.append("malformed Bose records: "+str(exc))
    supporting_records.append({"file":BOSE_DIR+"summary.json", "kernel_points":1216480,
        "kind":"Bose direct/exchange forms and short-time model witness; three raw grids independently replayed; initial 160% pointwise error retained and selected high-precision fallback checked; no total accuracy certification"})

    growing_rel = "results/stage3/growing-box/diagnostics.json"
    try:
        from audit_ds3_recent_stages import audit as audit_recent
        recent_errors, recent_records = audit_recent(ROOT)
        errors.extend(recent_errors)
        supporting_records.extend(recent_records)
    except (KeyError, TypeError, ValueError, ArithmeticError) as exc:
        errors.append("malformed post-Bose stage records: "+str(exc))

    growing = read_json(ROOT/growing_rel)
    try:
        from audit_ds3_uniform_local import audit as audit_uniform_local
        uniform_record = audit_uniform_local(ROOT)
        errors.extend("uniform-local: "+e for e in uniform_record["errors"])
        if uniform_record["checks"] != 446:
            errors.append("uniform-local: unexpected check count")
        supporting_records.append({"stage":"ds3-uniform-local", "checks":uniform_record["checks"],
            "status":uniform_record["status"], "kind":"uniform merged tails and local dr test integrals; no full physical mass integral"})
    except (KeyError, TypeError, ValueError, ArithmeticError) as exc:
        errors.append("malformed uniform-local stage records: "+str(exc))
    try:
        from audit_ds3_physical_local import audit as audit_physical_local
        physical_record = audit_physical_local(ROOT)
        errors.extend("physical-local: "+e for e in physical_record["errors"])
        supporting_records.append({"stage":"ds3-physical-local", "checks":physical_record["checks"],
            "status":physical_record["status"], "kind":"independent Decimal physical-weight replay on existing U nodes; local cutoff only"})
    except (KeyError, TypeError, ValueError, ArithmeticError) as exc:
        errors.append("malformed physical-local stage records: "+str(exc))
    try:
        from audit_ds3_compact_global import audit as audit_compact_global
        compact_replay = audit_compact_global(ROOT)
        compact_refined = read_json(ROOT/COMPACT_GLOBAL_DIR/"audit-refined.json")
        if (compact_replay["errors"] != ["central coefficient induction"]
                or compact_replay["checks"] != 442
                or compact_refined["status"] != "passed"
                or compact_refined["checks"] != 444 or compact_refined["errors"]):
            errors.append("compact-global: inconsistent original/refined audit")
        supporting_records.append({"stage":"ds3-compact-global", "checks":444,
            "status":compact_refined["status"],
            "kind":"compact kernel and Hilbert-tail replay; original structural-equality failure retained with exact-polynomial repair; no full physical quadrature"})
    except (KeyError, TypeError, ValueError, ArithmeticError) as exc:
        errors.append("malformed compact-global stage records: "+str(exc))
    try:
        from audit_ds3_fock_vacuum_refined import refined as replay_fock_vacuum
        from audit_ds3_gaussian_field import audit as replay_gaussian_field
        for stage, replay, count in [("ds3-fock-vacuum", replay_fock_vacuum, 154),
                                      ("ds3-gaussian-field", replay_gaussian_field, 115)]:
            record = replay(ROOT)
            errors.extend(stage+": "+e for e in record["errors"])
            if record["status"] != "passed" or record["checks"] != count:
                errors.append(stage+": unexpected replay status/count")
            supporting_records.append({"stage":stage, "checks":record["checks"], "status":record["status"],
                "kind":"regulated vacuum-return arithmetic and independent quadratures; no full fourfield or all-spin numerical certification"})
    except (KeyError, TypeError, ValueError, ArithmeticError, AssertionError) as exc:
        errors.append("malformed fock-vacuum stage records: "+str(exc))
    try:
        for name, count in [("diagnostics.json",52),("refinement.json",11)]:
            result = read_json(ROOT/GAUSSIAN_EIKONAL_DIR/name)
            if result["status"] != "passed" or len(result["checks"]) != count or not all(x["passed"] for x in result["checks"]):
                errors.append("gaussian-eikonal: invalid saved diagnostics "+name)
            for item in result["checks"]:
                if "error" in item and item["error"] > item["tolerance"]:
                    errors.append("gaussian-eikonal: criterion mismatch")
        study = read_json(ROOT/GAUSSIAN_EIKONAL_DIR/"study.json")
        refinement = read_json(ROOT/GAUSSIAN_EIKONAL_DIR/"refinement.json")
        if len(study["rows"]) != 8 or len(refinement["rows"]) != 8:
            errors.append("gaussian-eikonal: incomplete parameter rows")
        for old,new in zip(study["rows"],refinement["rows"]):
            if (old["gamma"],old["sigma"],old["u"]) != (new["gamma"],new["sigma"],new["u"]):
                errors.append("gaussian-eikonal: parameter mismatch")
            diff = abs(complex(*old["fine"]["value"])-complex(*new["value"]["value"]))
            if abs(diff-new["difference"]) > 1e-15 or diff >= 2e-7:
                errors.append("gaussian-eikonal: refinement replay")
        supporting_records.append({"stage":"gaussian-eikonal", "diagnostic_checks":52, "refinement_checks":11,
            "kind":"original momentum and four-coordinate controls plus finer quadrature; empirical errors, not a full gravity or principal-Fock identification"})
    except (KeyError, TypeError, ValueError, ArithmeticError) as exc:
        errors.append("malformed gaussian-eikonal records: "+str(exc))
    try:
        from audit_principal_eikonal import audit as replay_principal_eikonal
        replay = replay_principal_eikonal(ROOT)
        errors.extend("principal-eikonal: "+e for e in replay["errors"])
        if replay["status"] != "passed" or replay["checks"] != 727:
            errors.append("principal-eikonal: incomplete replay")
        caches = read_json(ROOT/PRINCIPAL_EIKONAL_DIR/"source-record.json")["sources"]
        if len(caches) != 3:
            errors.append("principal-eikonal: source count")
        for cache in caches:
            if digest(ROOT/cache["file"]) != cache["sha256"]:
                errors.append("principal-eikonal: source mismatch "+cache["file"])
        supporting_records.append({"stage":"principal-eikonal", "checks":727, "full_scrambled_runs":144,
            "kind":"raw arithmetic replay, recorded Bessel/rotation controls; empirical quadrature, late u8 sign unresolved"})
    except (KeyError, TypeError, ValueError, ArithmeticError) as exc:
        errors.append("malformed principal-eikonal records: "+str(exc))
    try:
        from audit_mass_filter import audit as replay_mass_filter
        replay = replay_mass_filter(ROOT)
        errors.extend('mass-filter: '+e for e in replay['errors'])
        if replay['status'] != 'passed' or replay['checks'] != 444:
            errors.append('mass-filter: incomplete replay')
        sources = read_json(ROOT/MASS_FILTER_DIR/'source-record.json')
        if len(sources['sources']) != 2 or len(sources['prior_evidence']) != 2:
            errors.append('mass-filter: source count')
        for cache in sources['sources']:
            if digest(ROOT/cache['file']) != cache['sha256']:
                errors.append('mass-filter: source mismatch '+cache['file'])
        for file, expected in sources['prior_evidence'].items():
            if digest(ROOT/file) != expected:
                errors.append('mass-filter: prior evidence mismatch '+file)
        supporting_records.append(dict(stage='mass-filter', checks=444, full_scrambled_runs=80,
            kind='Exact filter transform, conditional and assumed mass-average comparisons; not an established finite-kappa physical trace'))
    except (KeyError, TypeError, ValueError, ArithmeticError) as exc:
        errors.append('malformed mass-filter records: '+str(exc))
    growing_checks = growing.get("checks", [])
    if (growing.get("status") != "passed" or len(growing_checks) != 290
            or any(row.get("passed") is not True for row in growing_checks)
            or len(growing.get("gamma_ratio_rows", [])) != 16
            or len(growing.get("Euler_derivative_rows", [])) != 216
            or len(growing.get("energy_rows", [])) != 3
            or len(growing.get("cutoff_tail_rows", [])) != 27):
        errors.append("invalid growing-box proof-ingredient diagnostics")
    supporting_records.append({"file": growing_rel,
                               "kind": "finite-point Gamma/Euler/energy/tail ingredient checks; not uniform constants or exact finite-M packet certification",
                               "checks": len(growing_checks)})

    two_branch_rel = "results/stage3/two-branch/diagnostics.json"
    two_branch = read_json(ROOT/two_branch_rel)
    two_branch_checks = two_branch.get("checks", [])
    if (two_branch.get("status") != "passed" or len(two_branch_checks) != 300
            or any(row.get("passed") is not True for row in two_branch_checks)
            or len(two_branch.get("identities", [])) != 6
            or any(len(two_branch.get("runs", {}).get(name, [])) != 24
                   for name in ("w24-chi40", "w32-chi40", "w24-chi64"))):
        errors.append("invalid two-branch leading-packet diagnostics")
    supporting_records.append({"file": two_branch_rel,
                               "kind": "leading two-branch physical packet identities and empirical small-real refinements; not exact finite-M certification",
                               "checks": len(two_branch_checks)})

    density_rel = "results/stage3/density-endpoints/diagnostics.json"
    density = read_json(ROOT/density_rel)
    density_checks = density.get("checks", [])
    if (density.get("status") != "passed" or len(density_checks) != 384
            or any(row.get("passed") is not True for row in density_checks)
            or len(density.get("geometry_rows", [])) != 21
            or len(density.get("corner_rows", [])) != 20
            or len(density.get("gap_corner_rows", [])) != 5
            or len(density.get("refinement_rows", [])) != 75):
        errors.append("invalid density-endpoint diagnostics")
    for rel in DENSITY_RUNS:
        run = read_json(ROOT/rel)
        node_count = len(run.get("spectrum_nodes", {}).get("chi", []))
        if node_count <= 1000 or len(run.get("spectrum_nodes", {}).get("weights", [])) != node_count:
            errors.append(f"incomplete density-endpoint nodes: {rel}")
        for name in ("gamma1", "gamma2", "gap1"):
            model = run.get("models", {}).get(name, {})
            if (len(model.get("rows", [])) != 25 or len(model.get("JQ", [])) != node_count
                    or len(model.get("dense_curve", {}).get("T", [])) != 257):
                errors.append(f"incomplete density-endpoint model: {rel}/{name}")
    supporting_records.append({"file": density_rel,
                               "kind": "Gamma/gapped packet endpoint and full-curve diagnostics; empirical small-real precision limits retained",
                               "checks": len(density_checks)})

    orbit_identity = read_json(ROOT/ORBIT_IDENTITY_RECORD)
    identity_rows = orbit_identity.get("details")
    if (orbit_identity.get("status") != "passed" or orbit_identity.get("checks") != 51
            or orbit_identity.get("failed") != 0 or orbit_identity.get("working_dps") != 70
            or not isinstance(identity_rows, list) or len(identity_rows) != 51
            or any(not isinstance(row, dict) or row.get("passed") is not True for row in identity_rows)):
        errors.append("invalid orbit-core finite-point identity diagnostics")
    supporting_records.append({"file": ORBIT_IDENTITY_RECORD,
                               "kind": "finite-point special-function diagnostics; not interval or proof certification",
                               "checks": orbit_identity.get("checks"), "failed": orbit_identity.get("failed")})

    for rel in CG_VERTEX_DATA:
        data = read_json(ROOT/rel)
        expected = data.get("code_sha256")
        source = ROOT/"code/validate_vertex_source.py"
        if not isinstance(expected, str) or not source.is_file() or digest(source) != expected:
            errors.append(f"stale vertex source code: {rel}")
        checked_hashes += 1
        runs = data.get("runs")
        if (not isinstance(runs, list) or not runs
                or any(not isinstance(row, dict) or "exception" in row or "failure" in row for row in runs)):
            errors.append(f"invalid vertex data runs: {rel}")
        supporting_records.append({"file": rel, "kind": "raw numerical data; no passed-check claim",
                                   "runs": len(runs) if isinstance(runs, list) else None})

    analytic_sources = read_json(ROOT/CG_SOURCE_RECORD)
    analytic_records = analytic_sources.get("sources", [])
    if not isinstance(analytic_records, list) or not analytic_records:
        errors.append("missing analytic source archives")
        analytic_records = []
    analytic_paths = []
    for row in analytic_records:
        if not isinstance(row, dict) or not isinstance(row.get("file"), str):
            errors.append("invalid analytic source archive record")
            continue
        rel = row["file"].replace("\\", "/")
        analytic_paths.append(rel)
        source = ROOT/rel
        if (not source.is_file() or digest(source) != row.get("sha256")
                or source.stat().st_size != row.get("bytes")):
            errors.append(f"analytic source mismatch: {rel}")
        checked_hashes += 1
    for missing in sorted(CG_SOURCE_ARCHIVES - set(analytic_paths)):
        errors.append(f"missing analytic source provenance: {missing}")
    if len(analytic_paths) != len(set(analytic_paths)):
        errors.append("duplicate analytic source paths")
    supporting_records.append({"file": CG_SOURCE_RECORD, "kind": "versioned source archives",
                               "source_archives": len(analytic_paths)})

    sources = read_json(ROOT/"sources/manifest.json")
    source_records = sources.get("sources", [])
    if not source_records:
        errors.append("no source archives in manifest")
    archived_paths = {source["file"].replace("\\", "/") for source in source_records}
    for missing in sorted(REQUIRED_SOURCE_ARCHIVES - archived_paths):
        errors.append(f"source missing from manifest: {missing}")
    if len(archived_paths) != len(source_records):
        errors.append("duplicate source paths in manifest")
    for source in source_records:
        path = ROOT/source["file"].replace("\\", "/")
        if not path.is_file() or digest(path) != source["sha256"]:
            errors.append(f"source mismatch: {source['name']}")
    for item in sources.get("unavailable", []):
        warnings.append(f"Source remains a link, not an archive: {item['name']}")
    output = {
        "status": "passed" if not errors else "failed", "errors": errors, "warnings": warnings,
        "verified_at_utc": datetime.now(timezone.utc).isoformat(),
        "required_files": len(REQUIRED), "local_links_checked": checked_links,
        "code_and_input_hashes_checked": checked_hashes, "scientific_checks": scientific,
        "completed_studies": completed_studies,
        "historical_records_retained": historical,
        "supporting_data_and_provenance": supporting_records,
        "source_archives_checked": len(source_records),
        "scope": "File and link checks, source cache hashes and recorded code/input hashes. A fresh inventory is built on each run, not compared against a prior inventory; unbound CSV/image/result content is not independently certified. Snapshot mode verifies archived bytes against this inventory. Historical initial/failed outputs retain their old hashes. Scientific checks retain their own assumptions and tolerances; this audit does not certify a full holographic dual.",
    }
    (ROOT/"results/project-audit.json").write_text(json.dumps(output, indent=2)+"\n", encoding="utf-8")
    if errors:
        print(json.dumps(output, indent=2))
        raise SystemExit(1)

    files = sorted(p for p in ROOT.rglob("*") if p.is_file()
                   and not set(p.relative_to(ROOT).parts) & {"__pycache__", "snapshots", ".git", "tmp"}
                   and p.name != "project-files.json")
    inventory = [{"file": p.relative_to(ROOT).as_posix(), "bytes": p.stat().st_size,
                  "sha256": digest(p)} for p in files]
    manifest = ROOT/"results/project-files.json"
    manifest.write_text(json.dumps(inventory, indent=2)+"\n", encoding="utf-8")
    if args.snapshot:
        target.parent.mkdir(exist_ok=True)
        expected = {row["file"]: row["sha256"] for row in inventory}
        expected[manifest.relative_to(ROOT).as_posix()] = digest(manifest)
        with zipfile.ZipFile(target, "x", zipfile.ZIP_DEFLATED) as archive:
            for rel, sha in expected.items():
                content = (ROOT/rel).read_bytes()
                if hashlib.sha256(content).hexdigest() != sha:
                    raise AssertionError(f"file changed during freeze: {rel}")
                archive.writestr(rel, content)
        with zipfile.ZipFile(target) as archive:
            for rel, sha in expected.items():
                if hashlib.sha256(archive.read(rel)).hexdigest() != sha:
                    raise AssertionError(f"archive mismatch: {rel}")
        with target.with_suffix(".sha256").open("x", encoding="ascii") as checksum:
            checksum.write(digest(target)+"  "+target.name+"\n")
        print(f"Snapshot verified: {target} ({target.stat().st_size} bytes)")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
