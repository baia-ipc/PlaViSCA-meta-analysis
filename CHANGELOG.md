# Changelog

## Unreleased

### Added

- `scripts/integration_assessment.R` and the `integration-assessment` pixi
  task (`cd scripts && pixi run integration-assessment`). This is a
  read-only, deterministic post-integration assessment of the
  `plavisca-v1-dec01-dec02-2026-09-21` release object. It replaces the ad hoc
  computation behind the manuscript's integration metrics (Table 3, Fig 2).
  It writes the study/stage/cluster composition tables, global kNN-mixing and
  silhouette metrics, per-stage cross-study neighbour sharing (observed vs
  expected under complete mixing), and Fig 2 to
  `results/integration_assessment/`. Outputs are recorded in
  `build_manifest.tsv`. The script fails loudly if the input hash, the
  package library, or the reproduced values diverge from the published
  ones. See `results/README.md`.
