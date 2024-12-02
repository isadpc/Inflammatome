## may need to install
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("illuminaHumanv4.db", "hgug4112a.db", "hgu133plus2.db", 
# "AnnotationDbi", "hta20transcriptcluster.db", "HsAgilentDesign026652.db", "hgu219.db"))

# Run all scripts --------------------------------------------------------------
source(file = 'scripts/01_microarray_DGE.R')
source(file = 'scripts/02_RNAseq_DGE.R')
#source(file = 'scripts/04_benchmark_plot.R')