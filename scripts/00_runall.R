# Set WD -----------------------------------------------------------------------
setwd("..")

# Create directory to store figures --------------------------------------------
figures="figures/"
ifelse(!dir.exists(figures), dir.create(figures), FALSE)

# Run all scripts --------------------------------------------------------------
source(file = 'scripts/01_microarray_DGE.R')
source(file = 'scripts/02_RNAseq_DGE.R')
source(file = 'scripts/03_contrast_benchmark.R')
source(file = 'scripts/04_rank_aggregation.R')
source(file = 'scripts/05_inflammation_score_UseCase.R')
source(file = 'scripts/06_GSEA_UseCase.R')
source(file = 'scripts/07_overlap_networkAnnotations.R')
source(file = 'scripts/08_paper_figures.R')
