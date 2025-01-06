# Inflammatome Analysis Codebase

This repository contains the R code, data, and scripts for reproducing the analyses described in the paper **"[Your Paper Title Here]"**. The repository is organized to allow seamless execution of all analyses, from data preprocessing to figure generation.

---

## Repository Structure

- `data/`: Processed data and outputs from each analysis step.
  - `raw/`: Raw count matrices and metadata used in differential gene expression (DGE) analyses.
  - Outputs from each script are stored here with prefixes matching the number of their generating script (e.g., `01_microarray_DGE` for outputs from `01_microarray_DGE.R`).

- `scripts/`: Main R scripts for performing analyses.
  - `00_runall.R`: Master script to run all other scripts in sequence.
  - `01_microarray_DGE.R`: Differential gene expression analysis for microarray datasets.
  - `02_RNAseq_DGE.R`: Differential gene expression analysis for RNA-seq datasets.
  - `03_contrast_benchmark.R`: Benchmarking contrasts against inflammatory markers.
  - `04_rank_aggregation.R`: Rank aggregation to define the inflammatome and inflammation signature.
  - `05_inflammation_score_UseCase.R`: Use case demonstrating the calculation of the inflammation score.
  - `06_GSEA_UseCase.R`: Gene set enrichment analysis (GSEA) on selected contrasts.
  - `07_network_annotations.R`: Network annotation and analysis for inflammation-related genes.
  - `08_paper_figures.R`: Scripts to generate figures included in the paper.
  - `99_project_functions.R`: Custom functions used across multiple scripts.
  - `pre-analysis/`: Contains auxiliary scripts used before main analyses:
    - `gene_annotation.R`: Processes and annotates gene information.
    - `textmined_markers.R`: Retrieves and processes text-mined inflammatory markers.

- `figures/`: Figures generated during the analyses. Outputs are stored with prefixes matching the generating script.

---

## How to Use This Repository

### Prerequisites
Ensure the following are installed:
- R (version >= [version])
- Required R packages (listed in `00_runall.R`)

### Running the Analysis
1. Clone this repository to your local machine:
   ```bash
   git clone [repository URL]

