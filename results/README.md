# results/

Versioned analysis outputs derived **read-only** from released pipeline
artifacts. Nothing in this directory feeds back into the pipeline, and no
script that writes here modifies a release object.

## integration_assessment/

Post-integration assessment of the PlaViSca v1 release object
`deployments/data/releases/plavisca-v1-dec01-dec02-2026-09-21/pv_all_studies.rds`
(SHA256 `32068878fc92952e9be539befcc8ff2f42a6818716a912a55f5c45930cc830b9`,
105,949 cells). These are the numbers reported in the manuscript (Methods
"Assessment of integration", Results "Stage-resolved integration across
studies", Table 3, Fig 2).

Produced by `scripts/integration_assessment.R`.

### Regenerating

```bash
cd scripts
pixi run integration-assessment
```

The task runs in the locked pixi environment and clears `R_LIBS_USER`, so a
personal R library cannot shadow the locked package versions. The script
fails loudly (`pipeline_fail()`) under any of these conditions:

- the input's SHA256 differs from the release object above (override the
  path with `PLAVISCA_RELEASE_RDS`, but the hash check still applies);
- a package is loaded from outside the pixi environment library;
- any reproduced value deviates by more than 0.005 from the published value,
  or the study × stage totals differ.

Each run appends one row per output file to `build_manifest.tsv`. The row
records the output SHA256, git commit, R/Seurat versions and seed, and its
notes give the input path and SHA256 plus the Matrix/cluster/ggplot2/
patchwork/ggrastr/ragg versions.

### Files

| File | Contents |
|---|---|
| `study_by_stage.tsv` | `table(study_label, harmonized_life_cycle_stage, useNA = "ifany")` with margins |
| `study_by_broad_stage.tsv` | same for `parasite_broad_stage` |
| `cluster_by_study.tsv` | `RNA_snn_res.0.05` cluster × `study_label` |
| `cluster_by_stage.tsv` | `RNA_snn_res.0.05` cluster × `harmonized_life_cycle_stage` |
| `integration_metrics.tsv` | the four global metrics below, unintegrated vs integrated |
| `per_stage_cross_study_sharing.tsv` | per-stage cross-study neighbour sharing, observed vs expected |
| `fig2_integration.pdf` / `.png` | Fig 2 (A–D: UMAPs by study/stage before and after integration; E: per-stage observed/expected sharing) |

### Metrics

All metrics use Euclidean distances (`dist()`) on PCs 1–30 of
`pca_unintegrated` (log-normalised PCA before Harmony) and `pca_integrated`
(Harmony-corrected PCA). Neighbourhoods are the k = 30 nearest neighbours of
each cell, excluding the cell itself.

**Global metrics** use one subsample of 8,000 cells (`set.seed(1); sample(ncol, 8000)`).
The same cells are used for both embeddings.

| Metric | Definition | Better integration if |
|---|---|---|
| `knn_other_study` | mean over cells of the fraction of its 30 neighbours from a different `study_label` | higher (studies mix) |
| `knn_other_stage` | same with `harmonized_life_cycle_stage`; cells or neighbours with NA stage are dropped | lower (stages stay separate) |
| `silhouette_study` | mean silhouette width (`cluster::silhouette`) grouping cells by study | lower (study is not a cluster) |
| `silhouette_stage` | mean silhouette width grouping cells by stage (NA-stage cells dropped) | higher (stage is a cluster) |

**Per-stage cross-study sharing** is computed for each harmonized stage
present in more than one study. Stages are processed in the order of
`names(which(rowSums(table(stage, study) > 0) > 1))`, continuing the same
seeded RNG stream. For each stage:

1. Draw up to 3,000 of that stage's cells (`sample(sel, 3000)` when the stage
   has more).
2. Compute the mean fraction of each cell's 30 nearest *same-stage*
   neighbours that come from another study (`observed_*`).
3. Compute the value expected under complete mixing:
   `expected_complete_mixing = 1 - sum(p_i^2)`. Here `p_i` is study *i*'s
   share of *all* cells of that stage, not just the subsample.
4. Report `ratio_* = observed / expected` (1 = as mixed as study composition
   allows).

Interpret this as a within-stage comparison: it asks whether same-stage cells
from different studies become neighbours, independently of how many cells
each study contributes to the stage.
