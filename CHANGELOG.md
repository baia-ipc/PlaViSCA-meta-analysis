# Changelog

## Unreleased

### Changed

- `audit/phase3_preflight/final_release_manifest.tsv` and its generator no
  longer record absolute paths inside the removed
  `pre_process_data_gametocytes/` worktree. The study objects and the
  Ruberto2022_1 author object are now repo-relative (`silva2022.rds`,
  `data/Hep59.1.2.seu_20aug2025.rds`, …). The released outputs point at
  `../deployments/data/releases/plavisca-v1-gametocytes-2026-10-08/`
  (override with `PLAVISCA_RELEASE_DIR`/`PLAVISCA_CANDIDATE_ROOT`). All 13
  hashes and sizes are unchanged; AT23 and AT01–AT31 pass (31/31).
- The 2026-09-21 release and the phase-2 atlas are archived under
  `../archive/` (see `../archive/README.md`). `integration_assessment.R` and
  `validate_dec01_dec02_final.R` now say which archive to extract when their
  default input is missing.

### Notes

- `build_manifest.tsv` has nine rows timestamped 2026-10-08 12:14:29 whose
  artifact paths point to a temporary scratch directory
  (`/tmp/claude-1026/…/assess_published2/`) from an ad hoc check of the
  published release. The directory no longer exists. The rows are kept
  because the log is append-only.

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
