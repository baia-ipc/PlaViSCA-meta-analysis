# PlaViSca preprocessing pipeline

Builds the data of [PlaViSca](https://github.com/baia-ipc/PlaViSca), an atlas of
six public *Plasmodium vivax* single-cell RNA-seq studies, from the count
matrices of each study. It produces:

- the data files of the PlaViSca app;
- the integrated atlas (Seurat object);
- the tables, values and figures of the manuscript.

| Study | PMID | Sample types | Counts |
|---|---|---|---|
| Sa2020 | [32365102](https://pubmed.ncbi.nlm.nih.gov/32365102/) | blood stages | STARsolo |
| Ruberto2022_1 | [36093191](https://pubmed.ncbi.nlm.nih.gov/36093191/) | liver stages | authors' matrix |
| Ruberto2022_2 | [35926062](https://pubmed.ncbi.nlm.nih.gov/35926062/) | salivary gland sporozoites | STARsolo |
| Hazzard2022 | [36525464](https://pubmed.ncbi.nlm.nih.gov/36525464/) | blood stages, salivary gland sporozoites | STARsolo |
| Hazzard2024 | [39223117](https://pubmed.ncbi.nlm.nih.gov/39223117/) | blood stages | STARsolo |
| Mancio-Silva2022 | [35443155](https://pubmed.ncbi.nlm.nih.gov/35443155/) | liver stages | authors' matrix |

## Stages

Each stage is one script in `R/`. Each one reads the outputs of the stages
before it.

| Stage | What it does | Main outputs |
|---|---|---|
| `00_reference` | gene list, names and rRNA genes from the PlasmoDB-68 GFF | `genes.tsv` |
| `10_import_counts` | count matrices; which barcodes are cells | one object per study |
| `20_cell_metadata` | study, run and author annotations of each cell | one object per study |
| `30_integrate` | merge, normalise, Harmony, UMAP, t-SNE, clusters | `atlas_integrated.rds` |
| `35_blood_embedding` | embedding and clusters of the blood-stage cells alone | `blood.rds` |
| `40_annotate_stages` | life-cycle stage of each cell | `atlas.rds` |
| `50_export_app_data` | data files of the app | `results/app_data/` |
| `60_assess_integration` | integration metrics and figure | `results/integration/` |
| `70_tables_and_values` | tables, cited values, parameters, software versions, figures | `results/tables_and_values/` |

[docs/pipeline.md](docs/pipeline.md) describes each stage in detail. What has
been verified about the input data is in [docs/validation.md](docs/validation.md).

## Layout

```
R/                 stage scripts (NN_name.R), lib/ (shared code), tools/
starsolo/          from reads to STARsolo count matrices (own environment)
config/            params.yaml (analysis parameters), local.env.example
metadata/          curated tables: studies, runs, stages, markers, columns,
                   expected cell counts, input files and checksums
tests/             tests of the metadata and of the shared code
docs/              pipeline description, validation summary
run.sh             runs the stages
pixi.toml/.lock    software environment
```

Data are never stored in this directory. They live in a **workspace** whose
location is set in `config/local.env`:

```
<workspace>/inputs/          input data (read only)
  starsolo_counts/<PMID>/<run>_solo_out/   STARsolo output per run
  author_data/               authors' objects and tables
  reference/                 GFF, reference time course, orthologs
  study_count_objects/       per-study count objects (optional entry point)
<workspace>/runs/<run>/      one directory per run (PLAVISCA_RUN, default "main")
  outs/NN_stage/             outputs of each stage
  results/                   app data, tables, figures
  logs/                      one log per stage
  outputs.tsv                every file written, with checksum and code version
```

`metadata/input_files.tsv` lists the input files: where each comes from, its
size and its checksum.

## Running

Requirements: [pixi](https://pixi.sh) on Linux. The environment
(R 4.3, Seurat 5, Harmony, SingleR) is created from `pixi.lock` on first use.

```bash
cp config/local.env.example config/local.env   # then set PLAVISCA_WORKSPACE
./run.sh check-inputs                          # verify the input checksums
./run.sh all                                   # run every stage
./run.sh 40 50                                 # or run single stages
pixi run test                                  # tests (no data needed)
```

**Starting from the reads.** `starsolo/` makes the STARsolo count
matrices, the input of stage 10, from the reads of the runs in
`metadata/runs.tsv`. It reproduces the counts exactly; see
[starsolo/README.md](starsolo/README.md).

**Starting from count objects.** Importing the STARsolo matrices (stage 10)
is the slowest step. `R/tools/import_count_objects.R` can be used instead. It
starts from per-study count objects (`inputs/study_count_objects/`, checksums
in `metadata/study_count_objects.sha256`) and writes the same stage 10 outputs.

```bash
pixi run Rscript --vanilla R/tools/import_count_objects.R
pixi run Rscript --vanilla R/10_import_counts.R Ruberto2022_1 Mancio-Silva2022
./run.sh 20 30 40 50 60 70
```

## Checks

The stages stop with an explanation when something is not as expected. For
example:
- an input file is missing;
- a cell cannot be matched to its run or to the authors' data;
- a cell gets a stage that is not possible for its tissue;
- cell numbers differ from `metadata/expected_counts.tsv`.

To accept new cell numbers on purpose, update that file in the same commit
as the change that causes them.

## Contact

Sopheap Oeng (osopheap@pasteur-kh.org), Giorgio Gonnella
(ggonnella@pasteur-kh.org), BAIA, Institut Pasteur du Cambodge.
