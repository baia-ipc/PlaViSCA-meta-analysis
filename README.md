# PlaViSca — Pre-processing Pipeline

[![R](https://img.shields.io/badge/R-%3E%3D%204.3-blue?logo=r&logoColor=white)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-v5-brightgreen)](https://satijalab.org/seurat/)
[![SingleR](https://img.shields.io/badge/SingleR-annotation-orange)](https://bioconductor.org/packages/SingleR/)
[![License](https://img.shields.io/badge/license-CC--BY%204.0-lightgrey)](LICENSE)

This repository contains the full preprocessing pipeline used to convert raw STARsolo count matrices from six *Plasmodium vivax* single-cell RNA-seq studies into the analysis-ready RDS files consumed by the [PlaViSca](https://github.com/Sopheap15/PlaViSca) interactive Shiny application.

---

## Overview

The pipeline covers all steps from per-study quality control and Seurat object construction through cross-study integration, automated cell-type annotation, and final data flattening:

```
Raw counts (STARsolo)
        │
        ▼
Per-study QC & Seurat objects     ← scripts/[study]_pv_analysis_script.R
        │
        ▼
Cross-study merge & integration   ← scripts/integration.R
        │
        ▼
Cell-type annotation (SingleR)    ← scripts/singleR.R
        │
        ▼
Flatten to RDS for PlaViSca app   ← scripts/flatten_data.R
```

Supporting utilities:

| Script | Purpose |
|---|---|
| `scripts/gff_gene_name_extraction.R` | Parse PlasmoDB GFF to produce `ref/gff_data.rds` |
| `scripts/data_source_manipulation.R` | Build `data/data_source.csv` publication table |
| `scripts/hazzard2024_merge_all.R` | Merge per-batch Hazzard 2024 Seurat objects |
| `scripts/run_all_hazzard2024.sh` | Shell wrapper for batch Hazzard 2024 scripts |
| `scripts/integration_assessment.R` | Read-only integration metrics + Fig 2 for the release object (`pixi run integration-assessment`; see [results/README.md](results/README.md)) |

---

## Integrated Studies

| Study label | PMID | Journal | Cells | Sample type |
|---|---|---|---|---|
| Sa2020 | [32365102](https://pubmed.ncbi.nlm.nih.gov/32365102/) | PLoS Biol | 9 766 | Blood (mammalian host) |
| Ruberto2022_1 | [36093191](https://pubmed.ncbi.nlm.nih.gov/36093191/) | Front. Cell. Infect. Microbiol. | 1 438 | Hepatocytes (mammalian host) |
| Ruberto2022_2 | [35926062](https://pubmed.ncbi.nlm.nih.gov/35926062/) | PLoS Negl Trop Dis | 9 947 | Salivary gland (vector host) |
| Hazzard2022 | [36525464](https://pubmed.ncbi.nlm.nih.gov/36525464/) | PLoS Negl Trop Dis | 3 294 | Blood & salivary gland |
| Hazzard2024 | [39223117](https://pubmed.ncbi.nlm.nih.gov/39223117/) | Nat Commun | 80 024 | Blood (mammalian host) |
| Mancio-Silva2022 | [35443155](https://pubmed.ncbi.nlm.nih.gov/35443155/) | Cell Host & Microbe | 1 480 | Hepatocytes (human) |

---

## Repository Structure

```
.
├── scripts/
│   ├── ref/
│   │   └── PlasmoDB-68_PvivaxP01.gff      # P. vivax P01 gene annotation (GFF3)
│   ├── sa2020_pv_analysis_script.R         # Sa et al. 2020 — QC & Seurat object
│   ├── ruberto2022_1_pv_analysis_script.R  # Ruberto 2022 (hepatocytes) — QC & Seurat object
│   ├── ruberto2022_2_pv_analysis_script.R  # Ruberto 2022 (sporozoites) — QC & Seurat object
│   ├── hazzard2022_pv_analysis_script.R    # Hazzard 2022 — QC & Seurat object
│   ├── hazzard2024_pv_combined_*.R         # Hazzard 2024 — per-batch QC & Seurat objects
│   ├── hazzard2024_merge_all.R             # Merge all Hazzard 2024 batches
│   ├── run_all_hazzard2024.sh              # Shell wrapper to run Hazzard 2024 scripts
│   ├── silva.R                             # Mancio-Silva 2022 — Seurat object from processed data
│   ├── integration.R                       # Merge all studies, normalize, PCA, UMAP, t-SNE, Harmony
│   ├── singleR.R                           # SingleR annotation: IDC stage, parasite stages, gametocytes
│   ├── flatten_data.R                      # Extract metadata + embeddings + expression → RDS
│   ├── gff_gene_name_extraction.R          # Parse GFF3 → ref/gff_data.rds
│   ├── data_source_manipulation.R          # Build data_source.csv publication table
│   ├── pixi.toml                           # Pixi environment definition (R + Bioconductor packages)
│   └── pixi.lock                           # Locked dependency versions
│
├── counts/                                 # STARsolo GeneFull output directories (one per SRR)
│   ├── 32365102/                           # Sa2020 — SRR110082{69-78}
│   ├── 35926062/                           # Ruberto2022_2
│   ├── 36093191/                           # Ruberto2022_1
│   ├── 36525464/                           # Hazzard2022
│   └── 39223117/                           # Hazzard2024 / Mancio-Silva2022
│
├── data/                                   # Reference and supplementary input files
│   ├── GenesByOrthologs_Summary.xlsx       # PVX → PVP01 ortholog mapping (for SingleR)
│   ├── GenesByOrthologs_Summary.xlsx       # PVX → PVP01 gene ID conversion
│   ├── pbio.3000711.s034.xlsx              # Sa2020 gametocyte cell metadata
│   ├── Proccessed_Data.txt                 # Hazzard2024 processed cell barcodes
│   ├── PvData_CHM_Final.RDS               # Mancio-Silva2022 raw Seurat object
│   ├── TableS2_Cell_Metadata.xls           # Mancio-Silva2022 cell metadata
│   ├── Zhu_SciReps_2016.xls               # IDC timecourse reference (for SingleR)
│   └── 6_PvSPZ.BS.combined.rds            # Additional reference data
│
├── [study].rds                             # Per-study Seurat objects (intermediate outputs)
│   ├── hazzard2022.rds
│   ├── hazzard2024.rds
│   ├── ruberto2022_1.rds
│   ├── ruberto2022_2.rds
│   └── sa2020.rds
│
└── [study]_knee_plots.pdf                  # Knee/barcode rank plots (empty-droplet QC)
```

---

## Prerequisites

### Raw count data

Count matrices must be generated upstream using STARsolo. Two pipelines are provided:

- **Download FASTQ**: [download_SRA](https://github.com/Sopheap15/download_SRA) — downloads SRA data from NCBI.
- **Alignment & feature counting**: [sRNA_nextflow_starsolo](https://github.com/Sopheap15/sRNA_nextflow_starsolo) — aligns reads with STARsolo using the *P. vivax* P01 reference genome.

Place STARsolo outputs under `counts/<PMID>/<SRR_ID>_solo_out/Solo.out/GeneFull/`.

### Reference annotation

Download the *P. vivax* P01 GFF3 from [PlasmoDB](https://plasmodb.org/common/downloads/Current_Release/PvivaxP01/) (file: `PlasmoDB-68_PvivaxP01.gff`) and place it in `scripts/ref/`.

### R environment

Dependencies are managed with [Pixi](https://pixi.sh) (conda-based). Install Pixi, then from the `scripts/` directory:

```bash
pixi install
pixi shell
```

Key R packages (R ≥ 4.3):

| Package | Source | Purpose |
|---|---|---|
| Seurat (v5) | CRAN | scRNA-seq analysis |
| DropletUtils | Bioconductor | Empty droplet detection |
| SingleR | Bioconductor | Cell-type annotation |
| SingleCellExperiment | Bioconductor | SCE data structure |
| scater | Bioconductor | QC utilities |
| rtracklayer | Bioconductor | GFF3 parsing |
| tidyverse | CRAN | Data manipulation |
| janitor | CRAN | Column name cleaning |
| scCustomize | CRAN | Seurat plotting helpers |

Without Pixi, install manually:

```r
install.packages(c("Seurat", "tidyverse", "janitor", "scCustomize"))
BiocManager::install(c("DropletUtils", "SingleR", "SingleCellExperiment",
                       "scater", "rtracklayer"))
```

---

## Step-by-step Usage

### 1. Extract GFF gene annotations

```r
source("scripts/gff_gene_name_extraction.R")
# Output: ref/gff_data.rds
```

### 2. Per-study QC and Seurat object construction

Run each study script from the `pre_process_data/` root. Each script:

1. Loads STARsolo count matrices via `read10xCounts()`
2. Removes rRNA genes
3. Generates a knee/barcode-rank plot (PDF) for empty-droplet inspection
4. Calls `emptyDrops()` (where applicable) to filter empty droplets
5. Creates a `SeuratObject` and attaches study metadata
6. Saves a per-study `.rds` file

```bash
# From pre_process_data/
Rscript scripts/sa2020_pv_analysis_script.R
Rscript scripts/ruberto2022_1_pv_analysis_script.R
Rscript scripts/ruberto2022_2_pv_analysis_script.R
Rscript scripts/hazzard2022_pv_analysis_script.R

# Hazzard 2024 (multiple batches — run in parallel or sequentially)
bash scripts/run_all_hazzard2024.sh
Rscript scripts/hazzard2024_merge_all.R

# Mancio-Silva 2022 (uses pre-processed data object)
Rscript scripts/silva.R
```

### 3. Cross-study integration

`scripts/integration.R` merges all per-study Seurat objects and runs the full integration workflow:

- `NormalizeData()` + `FindVariableFeatures()` (2 000 HVGs) + `ScaleData()` + `RunPCA()`
- Unintegrated UMAP (3D, 30 dims) stored as `umap_unintegrated`
- Harmony integration (`IntegrateLayers()`, grouped by `study_label`) → `pca_integrated`
- Integrated UMAP (3D, 30 dims) stored as `umap_integrated`
- `FindNeighbors()` + `FindClusters()` (resolutions 0.01–0.10)
- t-SNE (3D, integrated PCA) stored as `tsne_integrated`
- Output: `pv_all_studies.rds`

```r
source("scripts/integration.R")
```

### 4. Cell-type annotation (SingleR)

`scripts/singleR.R` annotates cells using the Zhu et al. 2016 IDC timecourse as reference:

- Assigns **hours post invasion (HPI)** per blood-stage cell
- Derives **parasite stages**: Ring, Trophozoite, Schizont (blood), Merozoite, Hypnozoite, Schizont (liver), Sporozoite, Female gametocyte, Male gametocyte
- Derives **development phase**: Blood stages, Liver stages, Sporozoite
- Integrates confirmed gametocyte labels from Sa2020 supplementary data
- Overwrites `pv_all_studies.rds` with annotated metadata

```r
source("scripts/singleR.R")
```

### 5. Flatten to app-ready RDS files

`scripts/flatten_data.R` extracts the data layers from the final Seurat object and produces three flat data frames:

| Output file | Contents |
|---|---|
| `data/normalize_df.rds` | Metadata + embeddings + log-normalized expression (cells × genes) |
| `data/raw_df.rds` | Metadata + embeddings + raw UMI counts (cells × genes) |
| `data/scale_df.rds` | Metadata + embeddings + scaled expression (cells × genes) |
| `data/cleaned_dataset.rds` | List: `mr_data` (metadata+embeddings), `pca_df` (per-study PCA stdev), `hiv_data` (HVF info) |

```r
source("scripts/flatten_data.R")
# Copy outputs to PlaViSca/data/ for use by the app
```

### 6. Build the data source table (optional)

```r
source("scripts/data_source_manipulation.R")
# Output: data/data_source.csv
```

---

## Output Files for PlaViSca

Copy the following files to the `data/` directory of the PlaViSca Shiny app:

```
data/normalize_df.rds      → PlaViSca/data/normalize_df.rds
data/raw_df.rds            → PlaViSca/data/raw_df.rds
data/cleaned_dataset.rds   → PlaViSca/data/cleaned_dataset.rds
data/data_source.csv       → PlaViSca/data/data_source.csv
ref/gff_data.rds           → PlaViSca/ref/PvivaxP01_gff_data.rds
```

---

## Metadata Columns

The processed `mr_data` data frame (inside `cleaned_dataset.rds`) includes the following metadata columns:

| Column | Description |
|---|---|
| `study_label` | Study identifier (e.g., `Sa2020`, `Hazzard2022`) |
| `study_pmid` | PubMed ID |
| `run_id` | SRA run accession (SRR…) |
| `host_species` | Host organism (e.g., *Saimiri boliviensis*, *Homo sapiens*) |
| `sample_type` | Tissue of origin (blood, hepatocyte, salivary gland) |
| `strain` | Parasite strain |
| `treatment` | Drug treatment or `No_Treatment` |
| `hour_post_invasion` | Predicted IDC hour (SingleR, blood-stage cells only) |
| `parasite_stages` | Specific stage (9 categories, see below) |
| `development_phase` | Broad phase: Blood stages / Liver stages / Sporozoite |
| `sc_technology` | Sequencing technology (e.g., `10x_Chromium_V3`, `Seq_Well`) |
| `sequencer` | Sequencing instrument |
| `umap_i_1/2/3` | Integrated UMAP coordinates |
| `umap_u_1/2/3` | Unintegrated UMAP coordinates |
| `pca_i_1/2/3` | Integrated PCA (Harmony) coordinates |
| `pca_u_1/2/3` | Unintegrated PCA coordinates |
| `t_sne_1/2/3` | t-SNE coordinates (integrated) |

**Parasite stages (9 categories)**:

| Stage | Development phase |
|---|---|
| Ring stage | Blood stages |
| Trophozoite | Blood stages |
| Schizont (blood) | Blood stages |
| Merozoite | Blood stages |
| Female gametocyte | Blood stages |
| Male gametocyte | Blood stages |
| Hypnozoite | Liver stages |
| Schizont (liver) | Liver stages |
| Sporozoite | Sporozoite |

---

## Contributing

1. Open an issue to discuss major changes.
2. Fork the repository and create a feature branch.
3. Submit a pull request with a clear description and reproducible examples.

---

## Related Repositories

- [PlaViSca](https://github.com/Sopheap15/PlaViSca) — Interactive Shiny application
- [download_SRA](https://github.com/Sopheap15/download_SRA) — FASTQ download pipeline
- [sRNA_nextflow_starsolo](https://github.com/Sopheap15/sRNA_nextflow_starsolo) — Alignment pipeline

---

## Contact

For questions or support, open an issue or contact the maintainers:

- Sopheap Oeng — osopheap@pasteur-kh.org
- Giorgio Gonnella — ggonnella@pasteur-kh.org

---

© 2025 BAIA Team, Institut Pasteur du Cambodge. All rights reserved.
