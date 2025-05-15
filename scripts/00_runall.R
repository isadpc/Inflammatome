## may need to install
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("illuminaHumanv4.db", "hgug4112a.db", "hgu133plus2.db", 
# "AnnotationDbi", "hta20transcriptcluster.db", "HsAgilentDesign026652.db", "hgu219.db"))


figures="figures/"
resultdir="results/"
ifelse(!dir.exists(figures), dir.create(figures), FALSE)
ifelse(!dir.exists(resultdir), dir.create(resultdir), FALSE)


# Run all scripts --------------------------------------------------------------
source(file = 'scripts/01_microarray_DGE.R')
source(file = 'scripts/02_RNAseq_DGE.R')
source(file = 'scripts/03_contrast_benchmark.R')
source(file = 'scripts/04_rank_aggregation.R')
source(file = 'scripts/05_inflammation_score_UseCase.R')
source(file = 'scripts/06_GSEA_UseCase.R')
source(file = 'scripts/07_overlap_networkAnnotations.R')
source(file = 'scripts/08_paper_figures.R')