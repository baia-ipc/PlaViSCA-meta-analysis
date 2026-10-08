# plavisca-v1-gametocytes-2026-10-08

Release directory: `deployments/data/releases/plavisca-v1-gametocytes-2026-10-08/`
(atlas sha256 `e1a41126…`; `SHA256SUMS` alongside). Previous release
`plavisca-v1-dec01-dec02-2026-09-21` is unchanged.
Authoritative manifest: `audit/phase3_preflight/final_release_manifest.tsv`.
Acceptance: `audit/phase3_preflight/acceptance_matrix_AT01_AT31.tsv`, 31 PASS / 0 FAIL.

## Changes

1. **Gametocyte precedence (singleR.R, DEC07 resolved).** In the previous
   build, Part 5a filled `harmonized_life_cycle_stage` with the IDC bin for
   every eligible blood cell, and Part 8 only wrote the cluster call where
   the stage was still NA (`overwrite_target <- male_female_call &
   is.na(harmonized)`). As a result no gametocyte label could reach the
   field. Part 9 imported the Sa2020 labels into `source_sex_annotation` only.
   The order is now: source-defined label > cluster-level gametocyte call >
   IDC bin (`annotation_precedence.tsv`). `adjudicated_display_stage` equals
   `harmonized_life_cycle_stage`.
2. **rRNA removal (integration.R).** 49 rRNA genes, all carrying counts only in
   Mancio-Silva2022 (206,690 UMIs), are removed before normalisation. silva.R
   never filtered rRNA, and the feature union re-added those rows to every
   study. Atlas features: 6,860 → 6,811.
   The 19 protein-coding rRNA-processing genes are kept.
3. **Integration features (integration.R).** Measured feature sets are the
   quantification panels: the STARsolo panel, Mancio-Silva genes, and
   Ruberto2022_1 `author_present`. The integration features drop from 5,203
   to 4,709. 210 of the previous 2,000 HVGs were structurally zero in
   Ruberto2022_1; none are now.

## Verification

- Before the changes, the regenerated integration output plus the unchanged
  singleR.R reproduced the 2026-09-21 atlas exactly
  (`integration_reproduction.log`).
- The cell-ID set equals the previous release (105,949 cells;
  `final_vs_2026-09-21_membership.tsv`).
- 23,370 cells changed harmonized stage, all from an IDC bin to a gametocyte
  label (`final_vs_2026-09-21_stage_*.tsv`). The precedence-only intermediate
  (old embedding) changed 22,910 cells (`precedence_only_*`,
  `assessment_precedence_only/`).
- Per-resolution gametocyte cluster tables: `resolution_release_2026-09-21/`
  (old atlas) and `resolution_new_atlas/`. Production resolution is still
  0.05. On the new atlas, every resolution from 0.02 to 0.1 agrees with 0.05
  on ≥ 99.5% of cells.
- App (PlaViSca branch `gametocyte-harmonized-stage`): 75/75 tests, startup
  contract 33/33, Shiny app object and writexl smoke export pass.
