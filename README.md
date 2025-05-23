# Mining Inflammation in Human Transcriptomes Codebase

This repository contains the R code, data, and scripts for reproducing the analyses described in the paper **"Mining inflammation in human transcriptomes: a consensus signature across diseases and tissues"**. The repository is organized to allow seamless execution of all analyses, from data preprocessing to figure generation.

---

## Repository Structure

- `data/`: Processed data and outputs from each analysis step.
  - `raw/`: Raw count matrices and metadata used in differential gene expression (DGE) analyses. Download links and expected paths for data hosted in GEO and SRA are provided in `raw_data_downloads` for reproducibility. Alternatively, the analysis can be run starting with script 03_.
  - Outputs from each script are stored here with prefixes matching the number of their generating script (e.g., `01_` for outputs from `01_microarray_DGE.R`).

- `scripts/`: Main R scripts for performing analyses.
  - `00_runall.R`: Master script to run all other scripts in sequence.
  - `01_microarray_DGE.R`: Differential gene expression analysis for microarray datasets.
  - `02_RNAseq_DGE.R`: Differential gene expression analysis for RNA-seq datasets.
  - `03_contrast_benchmark.R`: Benchmarking contrasts against a gold standard set of inflammatory markers.
  - `04_rank_aggregation.R`: Rank aggregation to define the inflammatome and inflammation signature.
  - `05_inflammation_score_UseCase.R`: Use case demonstrating the correlation of the inflammation signature-based score to disease severity scores in a ulcerative colitis proteomic dataset.
  - `06_GSEA_UseCase.R`: Gene set enrichment analysis (GSEA) with inflammation signature on selected contrasts and volcano plot showcasing the inflammatome on DE results.
  - `07_network_annotations.R`: Exploring overlaps with a mouse inflammatome and the GOBP inflammatory response set and generating network annotations for use in cytoscape.
  - `08_paper_figures.R`: Scripts to generate figures included in the paper.
  - `99_project_functions.R`: Custom functions.
  - `pre-analysis/`: Contains auxiliary scripts used before main analyses. The output of these scripts is frozen for reproducibility:
    - `gene_annotation.R`: Retrieving list of protein-coding genes to which results are mapped throughout the project for consistency. Output is `data/Final_Annotation_List_Biomart.tsv`
    - `textmined_markers.R`: Retrieves text-mined inflammatory markers that will be used as gold standard in combination with the MSigDB Hallmarks inflammatory response gene set. Output is `data/Inflammation_markers_Top100.tsv`

- `figures/`: Figures generated during the analyses. Outputs are stored with prefixes matching the generating script.

---

## How to Use This Repository

### Prerequisites
Ensure the following are installed:
- **R** (version ≥ 4.4.2 recommended)  
- **RStudio** (optional but recommended)  
- **Git** (to clone the repository)  
- The `renv` R package:
  
  ```r
  install.packages("renv")

### Running the Analysis
1. Clone this repository to your local machine:
    
   ```bash
   git clone https://github.com/isadpc/Inflammatome.git
   cd Inflammatome
2. Open the repository folder in R or RStudio.
3. Restore the package environment using `renv`:
   
   ```r
   renv::restore()
4. Run the full analysis pipeline by sourcing the master script:

  ```r
  source("scripts/00_runall.R")


