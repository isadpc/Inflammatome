# Clear workspace --------------------------------------------------------------
rm(list = ls())

# Load packages ----------------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(sjmisc)
library(janitor)
library(GEOquery)
library(stringr)

# Define working directory & load external functions ---------------------------
wdir <- paste0(dirname(dirname(rstudioapi::getSourceEditorContext()$path)), "/")
setwd(wdir)

source(file = "scripts/99_project_functions.R")


# Data set GSE109142 (transcriptomics), UC, recal mucosal biopsies (score) -----


#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109142

# As case-control study (unpaired, UC vs CTL)
## Count matrix generated from raw data by us

meta_raw.GSE109142 <- read_tsv(paste0("data/raw/", "SraRunTable.GSE109142.txt"))

meta.GSE109142 <- meta_raw.GSE109142 %>% 
  mutate(Condition = case_when(
    diagnosis == "Ulcerative Colitis" ~ "UC",
    TRUE ~ as.character(diagnosis)
  )) %>% 
  dplyr::rename("Sample_ID" = "Run") 

count.GSE109142 <- read_tsv(paste0("data/raw/", "cts.raw.GSE109142.txt")) 

#Analyze (diff expression: logFold changes)
# Raw counts fed into DESeq2

coldata_prel <- meta.GSE109142 %>% 
  column_to_rownames(var = "Sample_ID")
coldata_prel$Condition <- factor(coldata_prel$Condition)
coldata <- coldata_prel 

table(coldata$Condition) #Get overview of sample numbers per group

cts_prel <- count.GSE109142 %>%
  column_to_rownames(var = "Gene.name.ID")
cts <- as.matrix(cts_prel)

cts <- cts[, rownames(coldata)]
all(rownames(coldata) %in% colnames(cts)) #compare row/column names
all(rownames(coldata) == colnames(cts)) #compare row/column name order

# Input data already normalized? Check using boxplot (looks at max 10 random samples)
Boxplot_prior(cts, "02_boxplot_raw_GSE109142")

# DESeq2 procedure
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)

print(nrow(dds))
dds <- dds[rowSums(counts(dds)) >= 1,] 
print(nrow(dds))

## PCA 
scores <- PCA_norm(dds, "Condition", "02_PCA_GSE109142")

DEseq.output <- DESeq(dds)

# Dispersion estimate QC
plotDispEsts(DEseq.output)

# Counts normalized for library size (NOT transformed [log2, vst, etc])
cts.norm = as.data.frame(counts(DEseq.output, normalized = TRUE))

## Boxplot of DESeq2-normalized data (looks at max 10 random samples)
Boxplot_post(cts.norm, "02_boxplot_norm_GSE109142")

#Pairs to analyze: 
#CTL_UC

DEG_CTL_UC <- as.data.frame(results(DEseq.output, contrast = c("Condition", "UC", "Control"))) %>% 
  rownames_to_column(var = "Gene.name.ID") 

DEG_CTL_UC_Proc <- merge(DEG_CTL_UC %>% dplyr::rename("ENSG.ID" = "Gene.name.ID"), Final_Annotation_List)
write_tsv(DEG_CTL_UC_Proc, "data/02_GSE109142_UC_CTL.tsv")

benchmark_plot_RNA(DEG_CTL_UC_Proc, "GSE109142: UC vs Healthy Controls")


# Data set GSE198520 (transcriptomics), RA, synovial tissue (paired & score) ------------


#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE198520

#Study compares RA patients before & after anti-inflammatory treatment (anti-TNF)

samples.GSE198520 <- read_tsv(paste("data/raw/", "GSE198520_sample_info.txt",sep=""), col_names = c("Sample Name", "Sample_ID"))

meta_raw.GSE198520 <- read_csv(paste("data/raw/", "SraRunTable.GSE198520.txt",sep=""))

meta.GSE198520 <- meta_raw.GSE198520 %>% 
  merge(samples.GSE198520, by = "Sample Name") %>% 
  mutate(Condition = Timepoint) %>% 
  mutate(Patient_ID = str_extract(Sample_ID, ".+(?=_p)"))

count.GSE198520 <- read_tsv(paste("data/raw/", "GSE198520_Raw_gene_count_matrix.txt", sep="")) %>% 
  dplyr::rename("Gene.name.ID" = "GeneSymbol")
colnames(count.GSE198520) <- c("Gene.name.ID", str_extract(colnames(count.GSE198520[-1]), ".+(?=_bx)"))

#As case-control study (paired, pre vs post immune-modulatory treatment)

coldata_prel <- meta.GSE198520 %>% 
  column_to_rownames(var = "Sample_ID") %>% 
  filter(status != "Non-Responder") #removed from paired design: non-responders showed low inflammation grade (which is why the immune-modulatory therapy did not work for them)

coldata_prel$Patient_ID <- factor(coldata_prel$Patient_ID)
coldata_prel$Condition <- factor(coldata_prel$Condition)
coldata <- coldata_prel 

table(coldata$Condition) #Get overview of sample numbers per group

cts_prel <- count.GSE198520 %>% 
  column_to_rownames(var = "Gene.name.ID")
cts <- as.matrix(cts_prel)

cts <- cts[, rownames(coldata)]
all(rownames(coldata) %in% colnames(cts)) #compare row/column names
all(rownames(coldata) == colnames(cts)) #compare row/column name order

# Input data already normalized? Check using boxplot (looks at max 10 random samples)
Boxplot_prior(cts, "02_boxplot_raw_GSE198520")

# DESeq2 procedure
#paired samples
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Patient_ID + Condition)

print(nrow(dds))
dds <- dds[rowSums(counts(dds)) >= 1,] 
print(nrow(dds))

## PCA 
scores <- PCA_norm(dds, "Condition", "02_PCA_GSE198520")

DEseq.output <- DESeq(dds)

# Dispersion estimate QC
plotDispEsts(DEseq.output)

# Counts normalized for library size
cts.norm = as.data.frame(counts(DEseq.output, normalized = TRUE))

## Boxplot of DESeq2-normalized data (looks at max 10 random samples)
Boxplot_post(cts.norm, "02_boxplot_norm_GSE198520")

# Differential expression analysis

DEG_RA_pre_post <- as.data.frame(results(DEseq.output, contrast = c("Condition", "Pre", "Post"))) %>% 
  rownames_to_column(var = "Gene.name.ID") 

DEG_RA_pre_post_Proc <- merge(DEG_RA_pre_post %>% dplyr::rename("Gene.name" = "Gene.name.ID"), Final_Annotation_List)
write_tsv(DEG_RA_pre_post_Proc, "data/02_GSE198520_RA_post_treatment.tsv")

benchmark_plot_RNA(DEG_RA_pre_post_Proc, "GSE198520: RA pre vs post treatment (paired)")


# Data set GSE137344 (transcriptomics), IBD, ileum tissue -----------------


#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137344

samples.GSE137344 <- read_tsv(paste("data/raw/", "GSE137344_sample_info.txt", sep = ""), col_names = c("GEO_Accession (exp)", "Sample"))

meta_raw.GSE137344 <- read_csv(paste("data/raw/", "SraRunTable.GSE137344.txt" , sep = "")) 

meta.GSE137344 <- meta_raw.GSE137344 %>% 
  left_join(samples.GSE137344, by = "GEO_Accession (exp)") %>% 
  dplyr::rename(GEO_Accession = `GEO_Accession (exp)`) %>% 
  mutate(Condition = case_when(
    Diagnosis == "Crohn's Disease" ~ "Crohns_Disease",
    Diagnosis == "Ulcerative Colitis" ~ "Ulcerative_Colitis",
    TRUE ~ as.character(Diagnosis)
  )) %>% 
  filter(Condition %in% c("Crohns_Disease", "Ulcerative_Colitis", "Control")) #removed unrelated diseases/samples

rawcount.GSE137344 <- read_tsv(paste("data/raw/", "GSE137344_Raw_Counts.txt", sep = "")) %>% 
  dplyr::rename(Gene.name.ID = `...1`)

#Analyze (diff expression: logFold changes)
# Raw counts fed into DESeq2

# cts = downloaded raw data matrix
# coldata = data.frame with rownames same as colnames in cts; and in the same order

coldata_prel <- meta.GSE137344 %>%
  column_to_rownames(var = "Sample")
coldata_prel$Condition <- factor(coldata_prel$Condition)

coldata <- coldata_prel 

table(coldata$Condition) #Get overview of sample numbers per group

cts_prel <- rawcount.GSE137344 %>% 
  filter(!str_detect(Gene.name.ID, "^\\d+(?=-)")) %>% #removes dates & NA values in gene names
  column_to_rownames(var = "Gene.name.ID")
cts <- as.matrix(cts_prel)

cts <- cts[, rownames(coldata)]
all(rownames(coldata) %in% colnames(cts)) #compare row/column names
all(rownames(coldata) == colnames(cts)) #compare row/column name order

# Input data already normalized? Check using boxplot (looks at max 10 random samples)
Boxplot_prior(cts, "02_boxplot_raw_GSE137344")

# DESeq2 procedure
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)

print(nrow(dds))
dds <- dds[rowSums(counts(dds)) >= 1,] 
print(nrow(dds))

## PCA 
scores <- PCA_norm(dds, "Condition", "02_PCA_GSE137344")

DEseq.output <- DESeq(dds)

# Dispersion estimate QC
plotDispEsts(DEseq.output)

# Counts normalized for library size
cts.norm = as.data.frame(counts(DEseq.output, normalized = TRUE))

## Boxplot of DESeq2-normalized data (looks at max 10 random samples)
Boxplot_post(cts.norm, "02_boxplot_norm_GSE137344")

# log-fold changes (disease-specific)

DEG_CD <- as.data.frame(results(DEseq.output, contrast = c("Condition", "Crohns_Disease", "Control"))) %>% 
  rownames_to_column(var = "Gene.name.ID")

DEG_CD_Proc <- merge(DEG_CD %>% dplyr::rename("Gene.name" = "Gene.name.ID"), Final_Annotation_List)
write_tsv(DEG_CD_Proc, "data/02_GSE137344_CD_CTL.tsv")

benchmark_plot_RNA(DEG_CD_Proc, "GSE137344: CD vs Healthy Controls")

hist(DEG_CD_Proc$log2FoldChange)
hist(DEG_CD_Proc$stat)
plot(1:17512, sort(DEG_CD_Proc$stat))

DEG_UC <- as.data.frame(results(DEseq.output, contrast = c("Condition", "Ulcerative_Colitis", "Control"))) %>% 
  rownames_to_column(var = "Gene.name.ID")

DEG_UC_Proc <- merge(DEG_UC %>% dplyr::rename("Gene.name" = "Gene.name.ID"), Final_Annotation_List)
write_tsv(DEG_UC_Proc, "data/02_GSE137344_UC_CTL.tsv")

benchmark_plot_RNA(DEG_UC_Proc, "GSE137344: UC vs Healthy Controls")

hist(DEG_UC_Proc$log2FoldChange)
hist(DEG_UC_Proc$stat)
plot(1:17512, sort(DEG_UC_Proc$stat))
plot(sort(DEG_UC_Proc$stat, decreasing = T), 1:17512)

# Data set GSE166925 (transcriptomics), IBD, intestinal tissue -----------------
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166925

meta_raw.GSE166925 <- read_csv(paste("data/raw/", "SraRunTable.GSE166925.txt",sep=""))

meta_raw.GSE166925 %>% dplyr::count(PATIENT_ID, Inflammation_status, site_taken)

meta.GSE166925 <- meta_raw.GSE166925 %>% 
  filter(disease %in% c("CD", "UC", "Control") & qc == "include") %>% #removed unrelated diseases/samples & samples that did NOT pass QC
  group_by(Sample_ID) %>% 
  filter(Bases == max(Bases)) %>% #Two entries per patient -> selected the one with higher base number
  ungroup() %>% 
  unite("Condition", c(disease, Inflammation_status, site_taken), sep = "_", remove = FALSE)

#Previous tests: paired samples (CD_infl_uninfl_SI & CD_infl_uninfl_SI) provided weaker inflammatory signature for CD than comparing to healthy controls 
#INSTEAD: unpaired samples analyzed for CD & UC: CD_infl_CTL_LI, UC_infl_CTL_LI
#Removed other contrast because DESeq2 refits >15000 genes if their samples are included

# Load RNAseq counts
count.GSE166925 <- read_tsv(paste("data/raw/", "GSE166925_featurecounts_counts.txt", sep="")) %>% 
  dplyr::rename("Gene.name.ID" = "gene_id") 

#Analyze (diff expression: logFold changes)
# Raw counts fed into DESeq2

coldata_prel <- meta.GSE166925 %>% 
  column_to_rownames(var = "Sample_ID") #%>% 
#filter(Condition %in% c("CD_inflamed_large_intestine", "UC_inflamed_large_intestine", "Control_uninflamed_large_intestine")) #include only contrasts of interest (explanation above)

coldata_prel$Condition <- factor(coldata_prel$Condition)
coldata_prel$PATIENT_ID <- factor(coldata_prel$PATIENT_ID)
coldata <- coldata_prel 

table(coldata$Condition) #Get overview of sample numbers per group

cts_prel <- count.GSE166925 %>% 
  column_to_rownames(var = "Gene.name.ID")
colnames(cts_prel) <- str_extract(colnames(cts_prel), "\\d+")
cts <- as.matrix(cts_prel)

cts <- cts[, rownames(coldata)]
all(rownames(coldata) %in% colnames(cts)) #compare row/column names
all(rownames(coldata) == colnames(cts)) #compare row/column name order

# Input data already normalized? Check using boxplot (looks at max 10 random samples)
Boxplot_prior(cts, "02_boxplot_raw_GSE166925")

# DESeq2 procedure
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)

coldata <- coldata %>%
  group_by(disease) %>%
  mutate(patient = factor(as.integer(factor(PATIENT_ID, levels = unique(PATIENT_ID)))))

data.frame(coldata$Condition, coldata$disease, coldata$PATIENT_ID, coldata$patient)

dds_paired <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ patient + Condition)

print(nrow(dds))
dds <- dds[rowSums(counts(dds)) >= 1,] 
print(nrow(dds))
dds_paired <- dds_paired[rowSums(counts(dds_paired)) >= 1,] 

## PCA
scores <- PCA_norm(dds, "Condition", "02_PCA_GSE166925")

DEseq.output <- DESeq(dds) #replacing outliers and refitting for 15850 genes
# there are no outliers in PCA
DEseq.output <- DESeq(dds, minReplicatesForReplace = Inf)
#DEseq.paired.output <- DESeq(dds_paired) # 1067 rows did not converge in beta
#dds_paired <- estimateSizeFactors(dds_paired)
#dds_paired <- estimateDispersions(dds_paired)
#DEseq.paired.output <- nbinomWaldTest(dds_paired, maxit=1000) #1067 rows did not converge in beta
# leave paired out

# Dispersion estimate QC
plotDispEsts(DEseq.output)
#plotDispEsts(DEseq.paired.output)

# Counts normalized for library size -> split by disease
cts.norm = as.data.frame(counts(DEseq.output, normalized = TRUE))

## Boxplot of DESeq2-normalized data (looks at max 10 random samples)
Boxplot_post(cts.norm, "02_boxplot_norm_GSE166925")

# log-fold changes (condition-specific)
#Paired samples for CD perform worse than comparing to healthy controls (CD_infl_uninfl_SI & CD_infl_uninfl_SI)
#INSTEAD: unpaired samples analyzed for CD & UC: CD_infl_CTL_LI, UC_infl_CTL_LI

CD_CTLsmall <- as.data.frame(results(DEseq.output, contrast = c("Condition", "CD_inflamed_small_intestine", "Control_uninflamed_large_intestine"), cooksCutoff = F)) %>% 
  rownames_to_column(var = "Gene.name.ID") 
CD_CTLlarge <- as.data.frame(results(DEseq.output, contrast = c("Condition", "CD_inflamed_large_intestine", "Control_uninflamed_large_intestine"), cooksCutoff = F)) %>% 
  rownames_to_column(var = "Gene.name.ID")
UC_CTL <- as.data.frame(results(DEseq.output, contrast = c("Condition", "UC_inflamed_large_intestine", "Control_uninflamed_large_intestine"), cooksCutoff = F)) %>% 
  rownames_to_column(var = "Gene.name.ID")

GSE166925_DEG_sets <- c("CD_CTLsmall", "CD_CTLlarge", "UC_CTL")

#Adding identifiers

DEG_Proc_list <- list()

for(dataset in GSE166925_DEG_sets){
  DEG_Proc <- merge(get(dataset) %>% dplyr::rename("ENSG.ID" = "Gene.name.ID"), Final_Annotation_List)
  write_tsv(DEG_Proc, paste0("data/02_GSE166925_", dataset, ".tsv"))
  
  DEG_Proc_list[[dataset]] <- DEG_Proc
}

benchmark_plot_RNA(DEG_Proc_list[["CD_CTLsmall"]], "GSE166925: CD inflamed vs non-tumor tissue [small intestine]")
benchmark_plot_RNA(DEG_Proc_list[["CD_CTLlarge"]], "GSE166925: CD inflamed vs non-tumor tissue [large intestine]")
benchmark_plot_RNA(DEG_Proc_list[["UC_CTL"]], "GSE166925: UC inflamed vs non-tumor tissue [large intestine]")


# Data set GSE121212 (transcriptomics), AD & PSO, skin punch biopsies (paired) --------

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE121212

samples.GSE121212 <- read_tsv(paste("data/raw/", "GSE121212_sample_info.txt",sep=""), col_names = c("GEO_Accession (exp)", "Sample"))

meta_raw.GSE121212 <- read_csv(paste("data/raw/", "SraRunTable.GSE121212.txt", sep = ""))

meta_GSE121212 <- meta_raw.GSE121212 %>% 
  left_join(samples.GSE121212, by = "GEO_Accession (exp)") %>% 
  dplyr::rename("GEO_Accession" = "GEO_Accession (exp)") %>% 
  separate_wider_delim(cols = "Sample", names = c("Disease", "Rest"), cols_remove = FALSE, delim = "_", too_many = "merge") %>% 
  separate_wider_delim(cols = "Rest", names = c("Patient", "lesion"), cols_remove = FALSE, delim = "_", too_many = "merge") %>% 
  unite("Patient_ID", c(Disease, Patient), remove = FALSE, sep = "_") %>% 
  unite("Condition_detail", c(Disease, Skin_Type), remove = FALSE, sep = "_") %>% 
  mutate(Condition = case_when(
    Condition_detail == "AD_chronic_lesion" ~ "AD_lesional",
    TRUE ~ as.character(Condition_detail)
  )) %>% 
  dplyr::select(-c(Patient, lesion, Rest)) %>% 
  mutate(Condition = str_replace_all(Condition, "-", "_"))

#ONLY paired samples analyzed: individually for AD_L_NL & PSO_L_NL
table(meta_GSE121212$Condition)

count.GSE121212 <- read_tsv(paste("data/raw/", "GSE121212_readcount.txt", sep = "")) %>% 
  dplyr::rename("Gene.name.ID" = "...1") 

cts_prel <- count.GSE121212 %>% 
  filter(!str_detect(Gene.name.ID, "(?<!.)\\d+(?=-)")) %>% #removes dates & NA values in gene names
  column_to_rownames(var = "Gene.name.ID")
cts <- as.matrix(cts_prel)

cts <- cts[, meta_GSE121212$Sample]
rownames(meta_GSE121212) <- meta_GSE121212$Sample
all(rownames(meta_GSE121212) %in% colnames(cts)) #compare row/column names
all(rownames(meta_GSE121212) == colnames(cts)) #compare row/column name order

## create deseq2 object + QC ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = meta_GSE121212,
                              design = ~ Condition)
dds


meta_GSE121212$Patient_ID <- factor(as.numeric(gsub("(PSO_|CTRL_|AD_)", "", meta_GSE121212$Patient_ID)))
dds_paired <- DESeqDataSetFromMatrix(countData = cts,
                                     colData = meta_GSE121212,
                                     design = ~ Patient_ID + Condition)
dds_paired

# Boxplot of Raw counts: Checking if counts are already normalized 
Boxplot_prior(cts, "02_boxplot_raw_GSE121212")

# filter low counts 
nrow(dds)
dds <- dds[rowSums(counts(dds)) >= 1,] 
nrow(dds)

# PCA 
scores <- PCA_norm(dds, "Condition", "02_PCA_GSE121212")

# Boxplot of size factor normalized counts
dds <- estimateSizeFactors(dds)
count.sf_norm <- counts(dds, normalized = T)
Boxplot_post(count.sf_norm, "02_boxplot_norm_GSE121212")

# Plot dispersion estimates (WIP)
dds <- estimateDispersions(dds)
plotDispEsts(dds) # cant figure out how to save this plot

## dge with deseq --------------------------------------------------------------
dds <- DESeq(dds, minReplicatesForReplace=Inf)
res1 <- as.data.frame(results(dds, contrast = c("Condition", "AD_lesional", "CTRL_healthy"),
                              cooksCutoff=FALSE)) %>%
  rownames_to_column(var = "Gene.name.ID")

res2 <- as.data.frame(results(dds, contrast = c("Condition", "PSO_lesional", "CTRL_healthy"),
                              cooksCutoff=FALSE)) %>%
  rownames_to_column(var = "Gene.name.ID")

#dds_paired <- DESeq(dds_paired) # need to rerun: 35 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest
dds_paired <- estimateSizeFactors(dds_paired)
dds_paired <- estimateDispersions(dds_paired)
dds_paired <- nbinomWaldTest(dds_paired, maxit=700)

res3 <- as.data.frame(results(dds_paired, contrast = c("Condition", "AD_lesional", "AD_non_lesional")))%>% 
  rownames_to_column(var = "Gene.name.ID")

res4 <- as.data.frame(results(dds_paired, contrast = c("Condition", "PSO_lesional", "PSO_non_lesional"))) %>%
  rownames_to_column(var = "Gene.name.ID")

## benchmark -------------------------------------------------------------------
res1 <- res1 %>% dplyr::rename("Gene.name" = "Gene.name.ID") %>% left_join(Final_Annotation_List)
res2 <- res2 %>% dplyr::rename("Gene.name" = "Gene.name.ID") %>% left_join(Final_Annotation_List)
res3 <- res3 %>% dplyr::rename("Gene.name" = "Gene.name.ID") %>% left_join(Final_Annotation_List)
res4 <- res4 %>% dplyr::rename("Gene.name" = "Gene.name.ID") %>% left_join(Final_Annotation_List)

benchmark_plot_RNA(filter(res1, Gene.type == "protein_coding"), "GSE121212: AD vs ctl")
benchmark_plot_RNA(filter(res2, Gene.type == "protein_coding"), "GSE121212: Pso vs ctl")
benchmark_plot_RNA(filter(res3, Gene.type == "protein_coding"), "GSE121212: AD paired")
benchmark_plot_RNA(filter(res4, Gene.type == "protein_coding"), "GSE121212: Pso paired")

# don't write data because I am using results of salmon in the end
#write_tsv(res1, "data/02_GSE121212_AD_Ctl.tsv")
#write_tsv(res2, "data/02_GSE121212_PSO_Ctl.tsv")
#write_tsv(res3, "data/02_GSE121212_AD_paired.tsv")
#write_tsv(res4, "data/02_GSE121212_PSO_paired.tsv")

# Data set GSE193309 (transcriptomics), AD, skin punch biopsies (paired) --------

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193309

samples.GSE193309 <- read_tsv(paste("data/raw/", "GSE193309_sample_info.txt",sep=""), col_names = c("GEO_Accession (exp)", "Sample"))

meta_raw.GSE193309 <- read_csv(paste("data/raw/", "SraRunTable.GSE193309.txt", sep = ""))

meta_raw.GSE193309 %>% dplyr::count(Subject_ID, Skin_Type, Visit_Date)

meta.GSE193309 <- meta_raw.GSE193309 %>% 
  filter(visit_id == "01") %>% 
  group_by(Subject_ID, Skin_Type, Visit_Date) %>% 
  filter(Bases == max(Bases)) %>% #Some have multiple entries -> selected the one with highest base number
  ungroup() %>% 
  merge(samples.GSE193309, by = "GEO_Accession (exp)") %>% 
  dplyr::rename("GEO_Accession" = "GEO_Accession (exp)") %>% 
  separate_wider_delim(cols = "Subject_ID", names = c("Disease", "Number"), cols_remove = FALSE, delim = "_") %>% 
  unite("Condition", c("Disease", "Skin_Type"), remove = FALSE, sep = "_") %>% 
  dplyr::select(-Number) %>% 
  dplyr::rename("Patient_ID" = "Subject_ID")

meta.GSE193309 %>% dplyr::count(Patient_ID, Skin_Type) %>% filter(n>1)

#ONLY paired samples analyzed: individually for AD_L_NL

count.GSE193309 <- read_csv(paste("data/raw/", "GSE193309_count_matrix.csv", sep = "")) %>% 
  dplyr::rename("Gene.name.ID" = "gene_name")   

#Analyze (diff expression: logFold changes)
# Raw counts fed into DESeq2

coldata_prel <- meta.GSE193309 %>% 
  #filter(Skin_Type != "HC") %>% #keep paired samples only
  column_to_rownames(var = "Sample")
coldata_prel$Condition <- factor(coldata_prel$Condition)
coldata_prel$Patient_ID <- factor(coldata_prel$Patient_ID)
coldata <- coldata_prel 

table(coldata$Condition) #Get overview of sample numbers per group

cts_prel <- count.GSE193309 %>%
  column_to_rownames(var = "Gene.name.ID")
cts <- as.matrix(cts_prel)

cts <- cts[, rownames(coldata)]
all(rownames(coldata) %in% colnames(cts)) #compare row/column names
all(rownames(coldata) == colnames(cts)) #compare row/column name order

# Input data already normalized? Check using boxplot (looks at max 10 random samples)
Boxplot_prior(cts, "02_boxplot_raw_GSE193309")

# DESeq2 procedure
# paired samples
# re factoring patient to avoid matrix full rank
data.frame(coldata$Patient_ID , coldata$Condition, coldata$Disease)

coldata <- coldata %>%
  group_by(Disease) %>%
  mutate(Patient_Num = factor(as.integer(factor(Patient_ID, levels = unique(Patient_ID)))))

print(dplyr::count(coldata, Patient_Num, Disease), n =56)

data.frame(coldata$Patient_ID , coldata$Condition, coldata$Disease, coldata$Patient_Num)
dds_paired <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Patient_Num + Condition)

# unpaired design
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)

print(nrow(dds))
dds <- dds[rowSums(counts(dds)) >= 1,] 
print(nrow(dds))

# PCA 
scores <- PCA_norm(dds, "Condition", "02_PCA_GSE193309")

DEseq.output <- DESeq(dds, minReplicatesForReplace=Inf) 
DEseq.paired.output <- DESeq(dds_paired) # 29 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest

# Remove gene row(s) without beta convergence (increasing maxit argument did not help)
DEseq.paired.output <- DEseq.paired.output[which(mcols(DEseq.paired.output)$betaConv),]

# Dispersion estimate QC
plotDispEsts(DEseq.output)
plotDispEsts(DEseq.paired.output)
# disp plots look weird

# Counts normalized for library size
cts.norm = as.data.frame(counts(DEseq.output, normalized = TRUE))

## Boxplot of DESeq2-normalized data (looks at max 10 random samples)
Boxplot_post(cts.norm, "02_boxplot_norm_GSE193309")

# log-fold changes (condition-specific)
#ONLY paired samples analyzed: individually for AD_L_NL

DEG_AD_L_NL <- as.data.frame(results(DEseq.paired.output, contrast = c("Condition", "AD_LS", "AD_NL"))) %>% 
  rownames_to_column(var = "Gene.name.ID")

DEG_AD_L_NL_Proc <- merge(DEG_AD_L_NL %>% dplyr::rename("Gene.name" = "Gene.name.ID"), Final_Annotation_List)
#write_tsv(DEG_AD_L_NL_Proc, "data/02_GSE193309_AD_paired.tsv")

benchmark_plot_RNA(DEG_AD_L_NL_Proc, "GSE193309: AD lesional vs non-lesional (paired)")

AD_CTL <- as.data.frame(results(DEseq.output, contrast = c("Condition", "AD_LS", "CO_HC"),
                                cooksCutoff = F)) %>% 
  rownames_to_column(var = "Gene.name.ID")
AD_CTL <- merge(AD_CTL %>% dplyr::rename("Gene.name" = "Gene.name.ID"), Final_Annotation_List)
#write_tsv(AD_CTL, "data/02_GSE193309_AD_Ctl.tsv")

# don't write data because I am using results of salmon in the end

benchmark_plot_RNA(AD_CTL, "GSE193309: AD lesional vs Ctl")

# Data set GSE138614 (transcriptomics), MS, white brain matter (autopsies) --------

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138614

samples.GSE138614 <- read_tsv(paste("data/raw/", "GSE138614_sample_info.txt",sep=""), col_names = c("GEO_Accession (exp)", "Sample"))

meta_raw.GSE138614 <- read_csv(paste("data/raw/", "SraRunTable.GSE138614.txt", sep = ""))

meta_raw.GSE138614 %>% dplyr::count(Individual, Lesion_type)

meta.GSE138614 <- meta_raw.GSE138614 %>% 
  #filter(Lesion_type %in% c("White matter (WM)", "Chronic active (CA)", "Active (AL)")) %>% #Select lesions of interest (= controls & inflamed)
  group_by(Individual) %>% 
  filter(Bases == max(Bases)) %>% #multiple lesions per MS patient -> pick the one with most sequenced bases
  ungroup() %>% 
  merge(samples.GSE138614, by = "GEO_Accession (exp)") %>% 
  dplyr::rename("GEO_Accession" = "GEO_Accession (exp)") %>% 
  mutate(Condition = case_when(
    Diagnosis == "Multiple sclerosis" & !(Lesion_type %in% c("Inactive (IL)", "Normal appearing white matter (NAWM)")) ~ "Active_ML",
    TRUE ~ as.character(Diagnosis)
  )) %>% 
  separate_wider_delim(cols = "Sample", names = c("Disease", "Sample_ID"), cols_remove = FALSE, delim = "_", too_many = "merge")

meta.GSE138614 %>% dplyr::count(Individual, Lesion_type)
table(meta.GSE138614$Condition)

count.GSE138614 <- read_tsv(paste("data/raw/", "GSE138614_countMatrix.txt", sep = "")) %>% 
  dplyr::rename("Gene.name.ID" = "...1")

colnames(count.GSE138614) <- c("Gene.name.ID", paste0("Sample_", str_extract(colnames(count.GSE138614)[-1], "(?<=G58-)\\d+")))

#Analyze (diff expression: logFold changes)
# Raw counts fed into DESeq2

coldata_prel <- meta.GSE138614 %>% 
  column_to_rownames(var = "Sample_ID")
coldata_prel$Condition <- factor(coldata_prel$Condition)
coldata <- coldata_prel 

table(coldata$Condition) #Get overview of sample numbers per group

cts_prel <- count.GSE138614 %>%
  column_to_rownames(var = "Gene.name.ID")
cts <- as.matrix(cts_prel)

cts <- cts[, rownames(coldata)]
all(rownames(coldata) %in% colnames(cts)) #compare row/column names
all(rownames(coldata) == colnames(cts)) #compare row/column name order

# Input data already normalized? Check using boxplot (looks at max 10 random samples)
Boxplot_prior(cts, "02_boxplot_raw_GSE138614")

# DESeq2 procedure
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)

print(nrow(dds))
dds <- dds[rowSums(counts(dds)) >= 1,] 
print(nrow(dds))

# PCA 
scores <- PCA_norm(dds, "Condition", "02_PCA_GSE138614")

DEseq.output <- DESeq(dds)

# Dispersion estimate QC
plotDispEsts(DEseq.output)

# Counts normalized for library size
cts.norm = as.data.frame(counts(DEseq.output, normalized = TRUE))
rownames(cts.norm) <- str_extract(rownames(cts.norm), ".+(?=\\.\\d+)")

## Boxplot of DESeq2-normalized data (looks at max 10 random samples)
Boxplot_post(cts.norm, "02_boxplot_norm_GSE138614")

#Pairs to analyze: 
#CTL_MS

DEG_CTL_MS_raw <- as.data.frame(results(DEseq.output, contrast = c("Condition", "Active_ML", "Control"))) %>% 
  rownames_to_column(var = "Gene.name.ID_raw") 

DEG_CTL_MS <- DEG_CTL_MS_raw %>% 
  mutate(Gene.name.ID = str_extract(DEG_CTL_MS_raw$Gene.name.ID_raw, ".+(?=\\.\\d+)")) %>% 
  dplyr::select(-Gene.name.ID_raw)

DEG_CTL_MS_Proc <- merge(DEG_CTL_MS %>% dplyr::rename("ENSG.ID" = "Gene.name.ID"), Final_Annotation_List)
write_tsv(DEG_CTL_MS_Proc, "data/02_GSE138614_MS_CTL.tsv")

benchmark_plot_RNA(DEG_CTL_MS_Proc, "GSE138614: MS autopsies vs Healthy Controls")


# MS vs inactive
DEG_CTL_MS_raw <- as.data.frame(results(DEseq.output, contrast = c("Condition", "Active_ML", "Multiple sclerosis"))) %>% 
  rownames_to_column(var = "Gene.name.ID_raw") 

DEG_CTL_MS <- DEG_CTL_MS_raw %>% 
  mutate(Gene.name.ID = str_extract(DEG_CTL_MS_raw$Gene.name.ID_raw, ".+(?=\\.\\d+)")) %>% 
  dplyr::select(-Gene.name.ID_raw)

DEG_CTL_MS_Proc <- merge(DEG_CTL_MS %>% dplyr::rename("ENSG.ID" = "Gene.name.ID"), Final_Annotation_List)
write_tsv(DEG_CTL_MS_Proc, "data/02_GSE138614_MS_MSinactive.tsv")

benchmark_plot_RNA(DEG_CTL_MS_Proc, "GSE138614: MS autopsies vs Healthy Controls")
# actually shows better enrichment

# Data set GSE123496 (transcriptomics), MS, different CNS regions (autopsies) --------


#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123496

samples.GSE123496 <- read_tsv(paste("data/raw/", "GSE123496_sample_info.txt",sep=""), col_names = c("GEO_Accession (exp)", "Sample"))

meta_raw.GSE123496 <- read_csv(paste("data/raw/", "SraRunTable.GSE123496.txt", sep = ""))

meta.GSE123496 <- meta_raw.GSE123496 %>% 
  merge(samples.GSE123496, by = "GEO_Accession (exp)") %>% 
  dplyr::rename("GEO_Accession" = "GEO_Accession (exp)") %>% 
  mutate(Condition = case_when(
    disease_state == "MS" ~ "Multiple_Sclerosis",
    disease_state == "healthy control" ~ "Control"
  )) 

count.GSE123496 <- read_csv(paste("data/raw/", "GSE123496_Human_MSNL_counts.m.csv", sep = "")) %>% 
  dplyr::rename("Gene.name.ID" = "...1")

#Analyze (diff expression: logFold changes)
# Raw counts fed into DESeq2

coldata_prel <- meta.GSE123496 %>% 
  column_to_rownames(var = "Sample")
coldata_prel$Condition <- factor(coldata_prel$Condition)
coldata <- coldata_prel 

table(coldata$Condition) #Get overview of sample numbers per group

cts_prel <- count.GSE123496 %>%
  column_to_rownames(var = "Gene.name.ID")
cts <- as.matrix(cts_prel)

cts <- cts[, rownames(coldata)]
all(rownames(coldata) %in% colnames(cts)) #compare row/column names
all(rownames(coldata) == colnames(cts)) #compare row/column name order

# Input data already normalized? Check using boxplot (looks at max 10 random samples)
Boxplot_prior(cts, "02_boxplot_raw_GSE123496")

# DESeq2 procedure
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)

print(nrow(dds))
dds <- dds[rowSums(counts(dds)) >= 1,] 
print(nrow(dds))

# PCA 
scores <- PCA_norm(dds, "Condition", "02_PCA_GSE123496")

DEseq.output <- DESeq(dds, minReplicatesForReplace=Inf)

# Dispersion estimate QC
plotDispEsts(DEseq.output)

# Counts normalized for library size
cts.norm = as.data.frame(counts(DEseq.output, normalized = TRUE))

## Boxplot of DESeq2-normalized data (looks at max 10 random samples)
Boxplot_post(cts.norm, "02_boxplot_norm_GSE123496")

#Pairs to analyze: 
#CTL_MS

DEG_CTL_MS <- as.data.frame(results(DEseq.output, contrast = c("Condition", "Multiple_Sclerosis", "Control"),
                                    cooksCutoff = F)) %>% 
  rownames_to_column(var = "Gene.name.ID") 

DEG_CTL_MS_Proc <- merge(DEG_CTL_MS %>% dplyr::rename("ENSG.ID" = "Gene.name.ID"), Final_Annotation_List)
write_tsv(DEG_CTL_MS_Proc, "data/02_GSE123496_MS_CTL.tsv")

benchmark_plot_RNA(DEG_CTL_MS_Proc, "GSE123496: MS autopsies vs Healthy Controls")

# Data set GSE126848 (transcriptomics), NASH, liver tissue --------


#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126848

samples.GSE126848 <- read_tsv(paste("data/raw/", "GSE126848_sample_info.txt",sep=""), col_names = c("GEO_Accession (exp)", "Sample"))

meta_raw.GSE126848 <- read_csv(paste("data/raw/", "SraRunTable.GSE126848.txt", sep = ""))

meta.GSE126848 <- meta_raw.GSE126848 %>% 
  merge(samples.GSE126848, by = "GEO_Accession (exp)") %>% 
  dplyr::rename("GEO_Accession" = "GEO_Accession (exp)") %>% 
  mutate(Condition = case_when(
    disease == "NASH" ~ "NASH",
    disease == "healthy" ~ "Control",
    disease == "NAFLD" ~ "NAFLD",
    disease == "obese" ~ "Obese"
  ))

count.GSE126848 <- read_tsv(paste("data/raw/", "GSE126848_Gene_counts_raw.txt", sep = "")) %>% 
  dplyr::rename("Gene.name.ID" = "key")

#Analyze (diff expression: logFold changes)
# Raw counts fed into DESeq2

coldata_prel <- meta.GSE126848 %>% 
  column_to_rownames(var = "Sample")
coldata_prel$Condition <- factor(coldata_prel$Condition)
coldata <- coldata_prel 

table(coldata$Condition) #Get overview of sample numbers per group

cts_prel <- count.GSE126848 %>%
  column_to_rownames(var = "Gene.name.ID")
cts <- as.matrix(cts_prel)

cts <- cts[, rownames(coldata)]
all(rownames(coldata) %in% colnames(cts)) #compare row/column names
all(rownames(coldata) == colnames(cts)) #compare row/column name order

# Input data already normalized? Check using boxplot (looks at max 10 random samples)
Boxplot_prior(cts, "02_boxplot_raw_GSE126848")

# DESeq2 procedure
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)

print(nrow(dds))
dds <- dds[rowSums(counts(dds)) >= 1,] 
print(nrow(dds))

# PCA 
scores <- PCA_norm(dds, "Condition", "02_PCA_GSE126848")
# very obvious batch effects
scores <- PCA_norm(dds, "gender", "02_PCA_gender_GSE126848")
# gender doesn't correspond with any clusters

DEseq.output <- DESeq(dds, minReplicatesForReplace = Inf)

# Dispersion estimate QC
plotDispEsts(DEseq.output)

# Counts normalized for library size
cts.norm = as.data.frame(counts(DEseq.output, normalized = TRUE))

## Boxplot of DESeq2-normalized data (looks at max 10 random samples)
Boxplot_post(cts.norm, "02_boxplot_norm_GSE126848")

#Pairs to analyze: 
#CTL_NASH

DEG_CTL_NASH <- as.data.frame(results(DEseq.output, contrast = c("Condition", "NASH", "Control"),
                                      cooksCutoff = F)) %>% 
  rownames_to_column(var = "Gene.name.ID") 

DEG_CTL_NASH_Proc <- merge(DEG_CTL_NASH %>% dplyr::rename("ENSG.ID" = "Gene.name.ID"), Final_Annotation_List)
write_tsv(DEG_CTL_NASH_Proc, "data/02_GSE126848_NASH_CTL.tsv")

benchmark_plot_RNA(DEG_CTL_NASH_Proc, "GSE126848: NASH vs Healthy Controls")

res <- as.data.frame(results(DEseq.output, contrast = c("Condition", "NAFLD", "Control"),
                             cooksCutoff = F)) %>% 
  rownames_to_column(var = "ENSG.ID") %>%
  left_join(Final_Annotation_List)
benchmark_plot_RNA(res, "NAFLD vs CTL")

res1 <- as.data.frame(results(DEseq.output, contrast = c("Condition", "Obese", "Control"),
                              cooksCutoff = F)) %>% 
  rownames_to_column(var = "ENSG.ID") %>%
  left_join(Final_Annotation_List)
benchmark_plot_RNA(res1, "Obese vs CTL")

write_tsv(res, "data/02_GSE126848_NAFLD_CTL.tsv")
write_tsv(res1, "data/02_GSE126848_OB_CTL.tsv")

# Data set GSE142530 (transcriptomics), AH, liver tissue -----------------------


#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142530

samples.GSE142530 <- read_tsv(paste("data/raw/", "GSE142530_sample_info.txt",sep=""), col_names = c("GEO_Accession (exp)", "Sample"))
samples.GSE142530$Sample <- str_replace_all(samples.GSE142530$Sample, " ", "_")

meta_raw.GSE142530 <- read_csv(paste("data/raw/", "SraRunTable.GSE142530.txt", sep = ""))

meta.GSE142530 <- meta_raw.GSE142530 %>% 
  merge(samples.GSE142530, by = "GEO_Accession (exp)") %>% 
  dplyr::rename("GEO_Accession" = "GEO_Accession (exp)") %>% 
  mutate(Condition = case_when(
    disease_state == "Alcoholic Hepatitis" ~ "Alcoholic_Hepatitis",
    disease_state == "Normal" ~ "Control",
    disease_state == "Alcohol-related Cirrhosis" ~ "Alcohol_related_Cirrhosis"
  ))

count.GSE142530_raw <- read_csv(paste("data/raw/", "GSE142530_Annoted-RNAseq-with-SampleIDs.csv", sep = ""))

count.GSE142530 <- count.GSE142530_raw %>% 
  dplyr::select(-gene_name) %>% 
  row_to_names(row_number = 1)
names(count.GSE142530)[1] <- "Gene.name.ID"
colnames(count.GSE142530) <- str_replace(colnames(count.GSE142530), "Lille\\s+", "")
colnames(count.GSE142530) <- str_replace(colnames(count.GSE142530), "JH-", "")
colnames(count.GSE142530) <- str_replace(colnames(count.GSE142530), "TPF\\s+", "")
colnames(count.GSE142530) <- str_replace(colnames(count.GSE142530), "not.used", "Control_120396")
colnames(count.GSE142530) <- str_replace_all(colnames(count.GSE142530), " ", "_")

all(samples.GSE142530$Sample %in% colnames(count.GSE142530))

#Analyze (diff expression: logFold changes)
# Raw counts fed into DESeq2

coldata_prel <- meta.GSE142530 %>% 
  column_to_rownames(var = "Sample")
coldata_prel$Condition <- factor(coldata_prel$Condition)
coldata <- coldata_prel 

table(coldata$Condition) #Get overview of sample numbers per group

cts_raw <- count.GSE142530 %>%
  filter(str_starts(Gene.name.ID, "ENSG")) %>%
  column_to_rownames(var = "Gene.name.ID") 

cts_prel <- mutate_all(cts_raw, as.numeric)

cts <- as.matrix(cts_prel)

cts <- cts[, rownames(coldata)]
all(rownames(coldata) %in% colnames(cts)) #compare row/column names
all(rownames(coldata) == colnames(cts)) #compare row/column name order

# Input data already normalized? Check using boxplot (looks at max 10 random samples)
Boxplot_prior(cts, "02_boxplot_raw_GSE142530")

# DESeq2 procedure
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)

print(nrow(dds))
dds <- dds[rowSums(counts(dds)) >= 1,] 
print(nrow(dds))

#smallestGroupSize <- min(count(coldata, Condition)$n)
#keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
#dds <- dds[keep,]
#nrow(dds)

counts(dds)

# PCA 
scores <- PCA_norm(dds, "Condition", "02_PCA_GSE142530")

DEseq.output <- DESeq(dds, minReplicatesForReplace = Inf)

# Dispersion estimate QC
plotDispEsts(DEseq.output)

# Counts normalized for library size
cts.norm = as.data.frame(counts(DEseq.output, normalized = TRUE))

## Boxplot of DESeq2-normalized data (looks at max 10 random samples)
Boxplot_post(cts.norm, "02_boxplot_norm_GSE142530")

#Pairs to analyze: 
#CTL_AH

DEG_CTL_AH <- as.data.frame(results(DEseq.output, contrast = c("Condition", "Alcoholic_Hepatitis", "Control"),
                                    cooksCutoff = F)) %>% 
  rownames_to_column(var = "Gene.name.ID") 

DEG_CTL_AH_Proc <- merge(DEG_CTL_AH %>% dplyr::rename("ENSG.ID" = "Gene.name.ID"), Final_Annotation_List)
write_tsv(DEG_CTL_AH_Proc, "data/02_GSE142530_AH_CTL.tsv")

benchmark_plot_RNA(DEG_CTL_AH_Proc, "GSE142530: AH vs Healthy Controls")

res <- as.data.frame(results(DEseq.output, contrast = c("Condition", "Alcohol_related_Cirrhosis", "Control"),
                             cooksCutoff = F)) %>% 
  rownames_to_column(var = "ENSG.ID") %>%
  left_join(Final_Annotation_List)

benchmark_plot_RNA(filter(res, Gene.type=="protein_coding"), "alcohol related cirrhosis")
write_tsv(res, "data/02_GSE142530_Cirr_CTL.tsv")

# Data set GSE124180 (transcriptomics), COPD, airway epithelium ----------------


#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124180

meta_raw.GSE124180 <- read_tsv(paste("data/raw/", "GSE124180_metadata.txt", sep = ""))

meta.GSE124180.proc <- as.data.frame(t(meta_raw.GSE124180)) %>% 
  rownames_to_column() %>% 
  row_to_names(row_number = 1)

for(n in c(colnames(meta.GSE124180.proc)[-c(1:6)])){
  meta.GSE124180.proc[, n] <- gsub(paste0(n, ": "), "", meta.GSE124180.proc[, n])
}

rownames(meta.GSE124180.proc)

meta.GSE124180 <- meta.GSE124180.proc %>% 
  mutate(Condition = paste(copd, `cell type`, sep=".")) %>%
  rownames_to_column()
#filter(`cell type` == "bronchial epithelium") %>% #removed unrelated cell types (= different batches)
#dplyr::rename("Condition" = "copd")

rownames(meta.GSE124180)

count.GSE124180_raw <- read_tsv(paste("data/raw/", "GSE124180_gene_count_table.tsv", sep = "")) 

count.GSE124180 <- count.GSE124180_raw %>% 
  filter(!row_number() %in% 1) %>% 
  dplyr::rename("Gene.name.ID" = "ENSEMBL_GENEID")

#Analyze (diff expression: logFold changes)
# Raw counts fed into DESeq2

coldata_prel <- meta.GSE124180 %>% 
  column_to_rownames(var = "Sample_ID")
coldata_prel$Condition <- factor(coldata_prel$Condition)
coldata <- coldata_prel 

table(coldata$Condition) #Get overview of sample numbers per group

cts_raw <- count.GSE124180 %>%
  column_to_rownames(var = "Gene.name.ID")

cts_prel <- mutate_all(cts_raw, as.numeric)

cts <- as.matrix(cts_prel)

cts <- cts[, rownames(coldata)]
all(rownames(coldata) %in% colnames(cts)) #compare row/column names
all(rownames(coldata) == colnames(cts)) #compare row/column name order

# Input data already normalized? Check using boxplot (looks at max 10 random samples)
Boxplot_prior(cts, "02_boxplot_raw_GSE124180")

# DESeq2 procedure
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)

print(nrow(dds))
dds <- dds[rowSums(counts(dds)) >= 1,] 
print(nrow(dds))

# PCA 
scores <- PCA_norm(dds, "Condition", "02_PCA_GSE124180")

DEseq.output <- DESeq(dds, minReplicatesForReplace = Inf)

# Dispersion estimate QC
plotDispEsts(DEseq.output)

# Counts normalized for library size
cts.norm = as.data.frame(counts(DEseq.output, normalized = TRUE))

## Boxplot of DESeq2-normalized data (looks at max 10 random samples)
Boxplot_post(cts.norm, "02_boxplot_norm_GSE124180")

#Pairs to analyze: 
#CTL_COPD

DEG_CTL_COPD <- as.data.frame(results(DEseq.output, contrast = c("Condition", "case.bronchial epithelium", "cont.bronchial epithelium"),
                                      cooksCutoff = Inf)) %>% 
  rownames_to_column(var = "Gene.name.ID") 

DEG_CTL_COPD_Proc <- merge(DEG_CTL_COPD %>% dplyr::rename("ENSG.ID" = "Gene.name.ID"), Final_Annotation_List)
write_tsv(DEG_CTL_COPD_Proc, "data/02_GSE124180_COPD_CTL.tsv")

benchmark_plot_RNA(DEG_CTL_COPD_Proc, "GSE124180: COPD vs Healthy Controls")


# Data set GSE100297 (transcriptomics), MS, optic chiasm ------------


#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100297

samples.GSE100297 <- read_tsv(paste("data/raw/", "GSE100297_sample_info.txt",sep=""), col_names = c("GEO_Accession (exp)", "Sample_ID"))

meta_raw.GSE100297 <- read_csv(paste("data/raw/", "SraRunTable.GSE100297.txt", sep = ""))

meta.GSE100297 <- meta_raw.GSE100297 %>% 
  merge(samples.GSE100297, by = "GEO_Accession (exp)") %>% 
  dplyr::rename("GEO_Accession" = "GEO_Accession (exp)") %>% 
  mutate(Condition = case_when(
    disease_state == "healthy control" ~ "Control",
    TRUE ~ as.character(disease_state)
  ))

count.GSE100297_CTL_raw <- read_tsv(paste("data/raw/", "GSE100297_hg19.gene_NL_Opt.txt", sep = ""))

count.GSE100297_MS_raw <- read_tsv(paste("data/raw/", "GSE100297_hg19.gene_MS_Opt.txt", sep = ""))

count.GSE100297 <- merge(count.GSE100297_CTL_raw, count.GSE100297_MS_raw, by = c("Gene.name.ID", "width")) %>% 
  dplyr::select(-width)

#Analyze (diff expression: logFold changes)
# Raw counts fed into DESeq2

coldata_prel <- meta.GSE100297 %>% 
  column_to_rownames(var = "Sample_ID")
coldata_prel$Condition <- factor(coldata_prel$Condition)
coldata <- coldata_prel 

table(coldata$Condition) #Get overview of sample numbers per group

cts_prel <- count.GSE100297 %>%
  column_to_rownames(var = "Gene.name.ID") 

cts <- as.matrix(cts_prel)

cts <- cts[, rownames(coldata)]
all(rownames(coldata) %in% colnames(cts)) #compare row/column names
all(rownames(coldata) == colnames(cts)) #compare row/column name order

# Input data already normalized? Check using boxplot (looks at max 10 random samples)
Boxplot_prior(cts, "02_boxplot_raw_GSE100297")

# DESeq2 procedure
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)

print(nrow(dds))
dds <- dds[rowSums(counts(dds)) >= 1,] 
print(nrow(dds))

# PCA 
scores <- PCA_norm(dds, "Condition", "02_PCA_GSE100297")

DEseq.output <- DESeq(dds)

# Dispersion estimate QC
plotDispEsts(DEseq.output)

# Counts normalized for library size
cts.norm = as.data.frame(counts(DEseq.output, normalized = TRUE))

## Boxplot of DESeq2-normalized data (looks at max 10 random samples)
Boxplot_post(cts.norm, "02_boxplot_norm_GSE100297")

#Pairs to analyze: 
#CTL_MS

DEG_CTL_MS <- as.data.frame(results(DEseq.output, contrast = c("Condition", "MS", "Control"))) %>% 
  rownames_to_column(var = "Gene.name.ID") 

DEG_CTL_MS_Proc <- merge(DEG_CTL_MS %>% dplyr::rename("ENSG.ID" = "Gene.name.ID"), Final_Annotation_List)
write_tsv(DEG_CTL_MS_Proc, "data/02_GSE100297_MS_CTL.tsv")

benchmark_plot_RNA(DEG_CTL_MS_Proc, "GSE100297: MS vs Healthy Controls")


# Data set GSE66511 (transcriptomics), PSO, skin punch biopsies (paired) --------


#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66511

samples.GSE66511 <- read_tsv(paste("data/raw/", "GSE66511_sample_info.txt",sep=""), col_names = c("GEO_Accession (exp)", "Sample_ID"))

meta_raw.GSE66511 <- read_csv(paste("data/raw/", "SraRunTable.GSE66511.txt", sep = ""))

meta.GSE66511 <- meta_raw.GSE66511 %>% 
  merge(samples.GSE66511, by = "GEO_Accession (exp)") %>% 
  dplyr::rename("GEO_Accession" = "GEO_Accession (exp)") %>% 
  dplyr::rename("Condition" = "disease") %>% 
  dplyr::rename("Patient_ID" = "Subject") %>%
  mutate(disease = if_else(Condition == "C", "Ctl", "Pso"))
#filter(Condition != "C") #remove healthy controls form paired design

count.GSE66511_raw <- read_tsv(paste("data/raw/", "GSE66511_Psoriasis_counts.txt", sep = ""))

count.GSE66511_temp <- count.GSE66511_raw %>% 
  dplyr::rename("Gene.name.ID" = "Symbol")
count.GSE66511_temp <- count.GSE66511_temp[rowSums(count.GSE66511_temp[2:37])>0,] #remove empty rows (responsible for many duplicate ENSG.IDs)

## Filter duplicates
count.GSE66511_dup <- count.GSE66511_temp %>% 
  group_by(Gene.name.ID) %>% 
  filter(n() > 1) %>% 
  ungroup() %>% 
  mutate(Count_sum = rowSums(.[2:37])) %>% 
  group_by(Gene.name.ID) %>% 
  filter(Count_sum == max(Count_sum)) %>% #select duplicate ENSG.ID with higher count sum across samples
  ungroup() %>% 
  dplyr::select(-Count_sum)

## Create final count table
count.GSE66511 <- count.GSE66511_temp %>% 
  group_by(Gene.name.ID) %>% 
  filter(n() == 1) %>% 
  ungroup() %>% 
  rbind(count.GSE66511_dup)

#Analyze (diff expression: logFold changes)
# Raw counts fed into DESeq2

coldata_prel <- meta.GSE66511 %>% 
  column_to_rownames(var = "Sample_ID")
coldata_prel$Condition <- factor(coldata_prel$Condition)
coldata_prel$Patient_ID <- factor(coldata_prel$Patient_ID)
coldata <- coldata_prel 

table(coldata$Condition) #Get overview of sample numbers per group

cts_prel <- count.GSE66511 %>%
  column_to_rownames(var = "Gene.name.ID") 

cts <- as.matrix(cts_prel)

cts <- cts[, rownames(coldata)]
all(rownames(coldata) %in% colnames(cts)) #compare row/column names
all(rownames(coldata) == colnames(cts)) #compare row/column name order

# Input data already normalized? Check using boxplot (looks at max 10 random samples)
Boxplot_prior(cts, "02_boxplot_raw_GSE66511")

# DESeq2 procedure
coldata <- coldata %>%
  group_by(disease) %>%
  mutate(Patient_Num = factor(as.integer(factor(Patient_ID, levels = unique(Patient_ID)))))

data.frame(coldata$Condition, coldata$Patient_ID, coldata$Patient_Num)

dds.paired <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Patient_Num + Condition)
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)

print(nrow(dds))
dds <- dds[rowSums(counts(dds)) >= 1,] 
print(nrow(dds))

# PCA 
scores <- PCA_norm(dds, "Condition", "02_PCA_GSE66511")

DEseq.output <- DESeq(dds, minReplicatesForReplace = Inf)
DEseq.paired.output <- DESeq(dds.paired)

# Dispersion estimate QC
plotDispEsts(DEseq.output)

# Counts normalized for library size
cts.norm = as.data.frame(counts(DEseq.output, normalized = TRUE))

## Boxplot of DESeq2-normalized data (looks at max 10 random samples)
Boxplot_post(cts.norm, "02_boxplot_norm_GSE66511")

#Pairs to analyze: 
#PSO_L_NL

DEG_PSO_L_NL <- as.data.frame(results(DEseq.paired.output, contrast = c("Condition", "LP", "NLP"))) %>% 
  rownames_to_column(var = "Gene.name.ID") 

DEG_PSO_L_NL_Proc <- merge(DEG_PSO_L_NL %>% dplyr::rename("Gene.name" = "Gene.name.ID"), Final_Annotation_List)
write_tsv(DEG_PSO_L_NL_Proc, "data/02_GSE66511_PSO_paired.tsv")

benchmark_plot_RNA(DEG_PSO_L_NL_Proc, "GSE66511: PSO lesional vs non-lesional (paired)")

# Unpaired
res <- as.data.frame(results(DEseq.output, contrast = c("Condition", "LP", "C"),
                             cooksCutoff = F)) %>% 
  rownames_to_column(var = "Gene.name.ID") 

res <- merge(res %>% dplyr::rename("Gene.name" = "Gene.name.ID"), Final_Annotation_List)
write_tsv(res, "data/02_GSE66511_PSO_Ctl.tsv")

benchmark_plot_RNA(res, "GSE66511: PSO vs control")

# Data set GSE197307 (transcriptomics), FSG/MN/MCD, glomeruli ---------


#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197307

meta_raw.GSE197307 <- read_tsv(paste("data/raw/", "GSE197307_series_matrix.txt", sep = ""), skip = 31)

meta_temp.GSE197307 <- meta_raw.GSE197307 %>% 
  mutate(`!Sample_title` = paste0(`!Sample_title`, "_", seq(1, nrow(.), 1))) %>% #assign unique names to columns
  rotate_df(cn = TRUE) %>% 
  rownames_to_column(var = "Sample_ID") %>% 
  dplyr::rename("GEO_Accession" = "!Sample_geo_accession_1") %>% 
  dplyr::rename("Disease" = "!Sample_characteristics_ch1_10")

# Remove "!" in column names
colnames(meta_temp.GSE197307) <- ifelse(str_detect(colnames(meta_temp.GSE197307), "!"), str_extract(colnames(meta_temp.GSE197307), "(?<=!).+"), colnames(meta_temp.GSE197307))

meta_temp2.GSE197307 <- meta_temp.GSE197307 %>% 
  mutate(Placeholder = str_replace(meta_temp.GSE197307$Disease, "disease: ", "")) %>% 
  mutate(Condition = case_when(
    Placeholder == "Healthy living transplant donor" ~ "Control",
    TRUE ~ as.character(Placeholder)
  ))

meta.GSE197307 <- meta_temp2.GSE197307 %>%  
  mutate(Condition = str_replace_all(meta_temp2.GSE197307$Condition, " ", "_")) #%>% 
#filter(!Condition %in% c("Nephrotic_Syndrome,_unspecified", "Other_Nephrotic_Syndrome")) #removed unrelated diseases/samples


count.GSE197307_raw <- read_tsv(paste("data/raw/", "GSE197307_NEPTUNE_Glom_RNA_counts_2022.txt", sep = ""))

count.GSE197307 <- count.GSE197307_raw %>% 
  mutate(Gene.name.ID = str_replace(count.GSE197307_raw$ENSEMBL_ene_ID, "ENS", "ENSG")) %>% 
  relocate(Gene.name.ID, .before = ENSEMBL_ene_ID) %>% 
  dplyr::select(-ENSEMBL_ene_ID)

#Analyze (diff expression: logFold changes)
# Raw counts fed into DESeq2

coldata_prel <- meta.GSE197307 %>% 
  column_to_rownames(var = "Sample_ID")
coldata_prel$Condition <- factor(coldata_prel$Condition)
coldata <- coldata_prel 

table(coldata$Condition) #Get overview of sample numbers per group

cts_prel <- count.GSE197307 %>%
  column_to_rownames(var = "Gene.name.ID") 

cts <- as.matrix(cts_prel)

cts <- cts[, rownames(coldata)]
all(rownames(coldata) %in% colnames(cts)) #compare row/column names
all(rownames(coldata) == colnames(cts)) #compare row/column name order

# Input data already normalized? Check using boxplot (looks at max 10 random samples)
Boxplot_prior(cts, "02_boxplot_raw_GSE197307")

# DESeq2 procedure
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)

print(nrow(dds))
dds <- dds[rowSums(counts(dds)) >= 1,] 
print(nrow(dds))

# PCA 
scores <- PCA_norm(dds, "Condition", "02_PCA_GSE197307")

DEseq.output <- DESeq(dds, minReplicatesForReplace = Inf)

# Dispersion estimate QC
plotDispEsts(DEseq.output)

# Counts normalized for library size
cts.norm = as.data.frame(counts(DEseq.output, normalized = TRUE))

## Boxplot of DESeq2-normalized data (looks at max 10 random samples)
Boxplot_post(cts.norm, "02_boxplot_norm_GSE197307")

# many 0 counts in normalized data
sum(rowSums(cts.norm) >= 282) # about half the genes have 0 counts in at least one sample so makes sense

# according to bioconductor post, might make sense to keep only normalized count of at least 10 in two or more rows
filter <- rowSums(cts.norm >= 10) >= 2
dds <- dds[filter,]

DEseq.output <- DESeq(dds, minReplicatesForReplace = Inf)

#Pairs to analyze: 
#CTL_Focal_Segmental_Glomerulosclerosis, CTL_Membranous_Nephropathy, CTL_Minimal_Change_Disease

FSGS_CTL <- as.data.frame(results(DEseq.output, 
                                  contrast = c("Condition", "Focal_Segmental_Glomerulosclerosis", "Control"),
                                  cooksCutoff = F)) %>% 
  rownames_to_column(var = "Gene.name.ID") 

MN_CTL <- as.data.frame(results(DEseq.output, contrast = c("Condition", "Membranous_Nephropathy", "Control"),
                                cooksCutoff = F)) %>% 
  rownames_to_column(var = "Gene.name.ID") 

MCD_CTL <- as.data.frame(results(DEseq.output, contrast = c("Condition", "Minimal_Change_Disease", "Control"),
                                 cooksCutoff = F)) %>% 
  rownames_to_column(var = "Gene.name.ID") 

GSE197307_DEG_sets <- c("FSGS_CTL", "MN_CTL", "MCD_CTL")

#Adding identifiers

DEG_Proc_list <- list()

for(dataset in GSE197307_DEG_sets){
  DEG_Proc <- merge(get(dataset) %>% dplyr::rename("ENSG.ID" = "Gene.name.ID"), Final_Annotation_List)
  write_tsv(DEG_Proc, paste0("data/02_GSE197307_", dataset, ".tsv"))
  
  DEG_Proc_list[[dataset]] <- DEG_Proc
}

benchmark_plot_RNA(DEG_Proc_list[["FSGS_CTL"]], "GSE197307: FSG vs Healthy Controls")
benchmark_plot_RNA(DEG_Proc_list[["MN_CTL"]], "GSE197307: MN vs Healthy Controls")
benchmark_plot_RNA(DEG_Proc_list[["MCD_CTL"]], "GSE197307: MCD vs Healthy Controls")

# Data set GSE84346 (transcriptomics), HCV, liver biopsies --------


#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84346

samples.GSE84346 <- read_tsv(paste("data/raw/", "GSE84346_sample_info.txt",sep=""), col_names = c("GEO_Accession (exp)", "Sample_ID"))

meta_raw.GSE84346 <- read_csv(paste("data/raw/", "SraRunTable.GSE84346.txt", sep = ""))

meta_raw.GSE84346 %>% dplyr::count(patient_identifier) %>% filter(n>1)
meta_raw.GSE84346 %>% dplyr::count(infection_status)

meta.GSE84346 <- meta_raw.GSE84346 %>% 
  merge(samples.GSE84346, by = "GEO_Accession (exp)") %>% 
  dplyr::rename("GEO_Accession" = "GEO_Accession (exp)") %>% 
  mutate(infection_status_mod = case_when(
    infection_status == "non-HCV" ~ "Control",
    TRUE ~ as.character(infection_status)
  )) 

meta.GSE84346$treatment_mod <- str_replace_all(meta.GSE84346$treatment, "/", "_")
meta.GSE84346$treatment_mod <- str_replace_all(meta.GSE84346$treatment_mod, " ", "_")
meta.GSE84346$treatment_bin <- if_else(meta.GSE84346$treatment_mod == "None", "non_treat", "treat")

meta.GSE84346 <- meta.GSE84346 %>% 
  mutate(Condition = paste0(infection_status_mod, "_", treatment_bin)) 

meta.GSE84346 %>% dplyr::count(Condition, patient_identifier)
# there are paired treated - untreated samples
meta.GSE84346 %>% dplyr::count(treatment_bin, patient_identifier)

count.GSE84346 <- read_tsv(paste("data/raw/", "GSE84346_UniqueReadCounts.txt", sep = "")) %>% 
  dplyr::rename("Gene.name.ID" = "GeneID")

#Analyze (diff expression: logFold changes)
# Raw counts fed into DESeq2

coldata_prel <- meta.GSE84346
coldata_prel$Condition <- factor(coldata_prel$Condition)
# refactor patient info
coldata_prel <- coldata_prel %>%
  group_by(infection_status) %>%
  mutate(Patient_Num = factor(as.integer(factor(patient_identifier, levels = unique(patient_identifier)))))
coldata <- coldata_prel %>% 
  column_to_rownames(var = "Sample_ID")

table(coldata$Condition) #Get overview of sample numbers per group
data.frame(coldata$infection_status, coldata$Condition, coldata$patient_identifier, coldata$Patient_Num)

cts_prel <- count.GSE84346 %>%
  column_to_rownames(var = "Gene.name.ID")
cts <- as.matrix(cts_prel)

cts <- cts[, rownames(coldata)]
all(rownames(coldata) %in% colnames(cts)) #compare row/column names
all(rownames(coldata) == colnames(cts)) #compare row/column name order

# Input data already normalized? Check using boxplot (looks at max 10 random samples)
Boxplot_prior(cts, "02_boxplot_raw_GSE84346")

# DESeq2 procedure
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)
dds.paired <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Patient_Num + Condition)

print(nrow(dds))
dds <- dds[rowSums(counts(dds)) >= 1,] 
print(nrow(dds))

dds.paired <- dds.paired[rowSums(counts(dds.paired)) >= 1,] 

# PCA 
scores <- PCA_norm(dds, "Condition", "02_PCA_GSE84346")

DEseq.output <- DESeq(dds, minReplicatesForReplace = Inf)

# DEseq.paired.output <- DESeq(dds.paired) #7 rows did not converge in beta
dds.paired <- estimateSizeFactors(dds.paired)
dds.paired <- estimateDispersions(dds.paired)
dds.paired <- nbinomWaldTest(dds.paired, maxit=1000)

# Remove gene row(s) without beta convergence (increasing maxit argument did not help)
dds.paired <- dds.paired[which(mcols(dds.paired)$betaConv),]

# Dispersion estimate QC
plotDispEsts(DEseq.output)
plotDispEsts(dds.paired)

# Counts normalized for library size
cts.norm = as.data.frame(counts(DEseq.output, normalized = TRUE))

## Boxplot of DESeq2-normalized data (looks at max 10 random samples)
Boxplot_post(cts.norm, "02_boxplot_norm_GSE84346")

#Pairs to analyze: 
#CTL_HCV

DEG_CTL_HCV <- as.data.frame(results(DEseq.output,
                                     contrast = c("Condition", "HCV_non_treat", "Control_non_treat"),
                                     cooksCutoff = F)) %>% 
  rownames_to_column(var = "Gene.name.ID") 

DEG_CTL_HCV_Proc <- merge(DEG_CTL_HCV %>% dplyr::rename("ENSG.ID" = "Gene.name.ID"), Final_Annotation_List)
write_tsv(DEG_CTL_HCV_Proc, "data/02_GSE84346_HCV_CTL.tsv")

benchmark_plot_RNA(DEG_CTL_HCV_Proc, "GSE84346: HCV (no treatment) vs Healthy Controls")

# paired
res <- as.data.frame(results(dds.paired,
                             contrast = c("Condition", "HCV_treat", "HCV_non_treat"),
                             cooksCutoff = F)) %>% 
  rownames_to_column(var = "Gene.name.ID") 

res <- merge(res %>% dplyr::rename("ENSG.ID" = "Gene.name.ID"), Final_Annotation_List)

benchmark_plot_RNA(res, "GSE84346: HCV non treated vs treated")
write_tsv(res, "data/02_GSE84346_HCV_paired.tsv")

# Data set GSE148036 (transcriptomics), TB/AC/SD, lung -------------------


#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148036

meta_raw.GSE148036 <- read_tsv(paste("data/raw/", "GSE148036_series_matrix.txt", sep = ""), skip = 38)

meta_temp.GSE148036 <- meta_raw.GSE148036 %>% 
  mutate(`!Sample_title` = paste0(`!Sample_title`, "_", seq(1, nrow(.), 1))) %>% #assign unique names to columns
  rotate_df(cn = TRUE) %>% 
  rownames_to_column(var = "Description") %>% 
  dplyr::rename("GEO_Accession" = "!Sample_geo_accession_1") %>% 
  dplyr::rename("Sample_ID" = "!Sample_source_name_ch1_7") %>% 
  dplyr::rename("Disease" = "!Sample_characteristics_ch1_12") %>% 
  relocate(Sample_ID, .after = Description) %>% 
  relocate(Disease, .after = Sample_ID)

# Remove "!" in column names
colnames(meta_temp.GSE148036) <- ifelse(str_detect(colnames(meta_temp.GSE148036), "!"), str_extract(colnames(meta_temp.GSE148036), "(?<=!).+"), colnames(meta_temp.GSE148036))

meta.GSE148036 <- meta_temp.GSE148036 %>% 
  mutate(Condition = str_replace(meta_temp.GSE148036$Disease, "disease: ", "")) %>% 
  relocate(Condition, .after = Sample_ID)

count.GSE148036_raw <- read_tsv(paste("data/raw/", "GSE148036_PRJNA609278count_matrix.txt", sep = ""))

count.GSE148036 <- count.GSE148036_raw %>% 
  dplyr::rename("Gene.name.ID" = "Geneid") %>% 
  filter(!is.na(Gene.name.ID))

#Analyze (diff expression: logFold changes)
# Raw counts fed into DESeq2

coldata_prel <- meta.GSE148036 %>% 
  column_to_rownames(var = "Sample_ID")
coldata_prel$Condition <- factor(coldata_prel$Condition)
coldata <- coldata_prel 

table(coldata$Condition) #Get overview of sample numbers per group

cts_prel <- count.GSE148036 %>%
  column_to_rownames(var = "Gene.name.ID") 
cts <- as.matrix(cts_prel)

cts <- cts[, rownames(coldata)]
all(rownames(coldata) %in% colnames(cts)) #compare row/column names
all(rownames(coldata) == colnames(cts)) #compare row/column name order

# Input data already normalized? Check using boxplot (looks at max 10 random samples)
Boxplot_prior(cts, "02_boxplot_raw_GSE148036")

# DESeq2 procedure
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)

print(nrow(dds))
dds <- dds[rowSums(counts(dds)) >= 1,] 
print(nrow(dds))

# PCA 
scores <- PCA_norm(dds, "Condition", "02_PCA_GSE148036")

DEseq.output <- DESeq(dds)

# Dispersion estimate QC
plotDispEsts(DEseq.output)

# Counts normalized for library size
cts.norm = as.data.frame(counts(DEseq.output, normalized = TRUE))

## Boxplot of DESeq2-normalized data (looks at max 10 random samples)
Boxplot_post(cts.norm, "02_boxplot_norm_GSE148036")

# filter lowly expressed genes
#filter <- rowSums(cts.norm >= 10) >= 2
#print(nrow(dds))
#dds <- dds[filter,]
#print(nrow(dds))

DEseq.output <- DESeq(dds)

#Pairs to analyze: 
#CTL_Tuberculosis, CTL_Sacrodosis

TB_CTL <- as.data.frame(results(DEseq.output, contrast = c("Condition", "Tuberculosis", "Normal"))) %>% 
  rownames_to_column(var = "Gene.name.ID") 

SD_CTL <- as.data.frame(results(DEseq.output, contrast = c("Condition", "Sacrodosis", "Normal"))) %>% 
  rownames_to_column(var = "Gene.name.ID") 

AC_CTL <- as.data.frame(results(DEseq.output, contrast = c("Condition", "Adenocarcinoma", "Normal"))) %>% 
  rownames_to_column(var = "Gene.name.ID") 

GSE148036_DEG_sets <- c("TB_CTL", "SD_CTL", "AC_CTL")

#Adding identifiers

DEG_Proc_list <- list()

for(dataset in GSE148036_DEG_sets){
  DEG_Proc <- merge(get(dataset) %>% dplyr::rename("ENSG.ID" = "Gene.name.ID"), Final_Annotation_List)
  write_tsv(DEG_Proc, paste0("data/02_GSE148036_", dataset, ".tsv"))
  
  DEG_Proc_list[[dataset]] <- DEG_Proc
}

benchmark_plot_RNA(DEG_Proc_list[["TB_CTL"]], "GSE148036: TB vs Healthy Controls")
benchmark_plot_RNA(DEG_Proc_list[["SD_CTL"]], "GSE148036: SD vs Healthy Controls")
benchmark_plot_RNA(DEG_Proc_list[["AC_CTL"]], "GSE148036: Adenocarcinoma vs Healthy Controls")


# Data set GSE150910 (transcriptomics), CHP/IPF, lung ---------------


#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150910

meta_raw.GSE150910 <- read_tsv(paste("data/raw/", "GSE150910_series_matrix.txt", sep = ""), skip = 30)

meta_temp.GSE150910 <- meta_raw.GSE150910 %>% 
  mutate(`!Sample_title` = paste0(`!Sample_title`, "_", seq(1, nrow(.), 1))) %>% #assign unique names to columns
  rotate_df(cn = TRUE) %>% 
  rownames_to_column(var = "Sample_ID") %>% 
  dplyr::rename("GEO_Accession" = "!Sample_geo_accession_1", "Disease" = "!Sample_characteristics_ch1_12") %>% 
  relocate(Disease, .before = GEO_Accession)

# Remove "!" in column names
colnames(meta_temp.GSE150910) <- ifelse(str_detect(colnames(meta_temp.GSE150910), "!"), str_extract(colnames(meta_temp.GSE150910), "(?<=!).+"), colnames(meta_temp.GSE150910))

meta.GSE150910 <- meta_temp.GSE150910 %>% 
  mutate(Condition = str_replace(meta_temp.GSE150910$Disease, "diagnosis: ", "")) %>% 
  relocate(Condition, .before = Disease) %>% 
  filter(Sample_ID != "ipf_H03") #different batch

count.GSE150910_raw <- read_csv(paste("data/raw/", "GSE150910_gene_level_count_file.csv", sep = ""))

count.GSE150910 <- count.GSE150910_raw %>% 
  dplyr::rename("Gene.name.ID" = "symbol") %>% 
  filter(!is.na(Gene.name.ID))

#Analyze (diff expression: logFold changes)
# Raw counts fed into DESeq2

coldata_prel <- meta.GSE150910 %>% 
  column_to_rownames(var = "Sample_ID")
coldata_prel$Condition <- factor(coldata_prel$Condition)
coldata <- coldata_prel 

table(coldata$Condition) #Get overview of sample numbers per group

cts_prel <- count.GSE150910 %>%
  column_to_rownames(var = "Gene.name.ID") 

cts <- as.matrix(cts_prel)

cts <- cts[, rownames(coldata)]
all(rownames(coldata) %in% colnames(cts)) #compare row/column names
all(rownames(coldata) == colnames(cts)) #compare row/column name order

# Input data already normalized? Check using boxplot (looks at max 10 random samples)
Boxplot_prior(cts, "02_boxplot_raw_GSE150910")

# DESeq2 procedure
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)

print(nrow(dds))
dds <- dds[rowSums(counts(dds)) >= 1,] 
print(nrow(dds))

# PCA 
scores <- PCA_norm(dds, "Condition", "02_PCA_GSE150910")

DEseq.output <- DESeq(dds, minReplicatesForReplace = Inf)

# Dispersion estimate QC
plotDispEsts(DEseq.output)

# Counts normalized for library size
cts.norm = as.data.frame(counts(DEseq.output, normalized = TRUE))

## Boxplot of DESeq2-normalized data (looks at max 10 random samples)
Boxplot_post(cts.norm, "02_boxplot_norm_GSE150910")

#Pairs to analyze: 
#CTL_Chronic_hypersensitivity_pneumonitis, CTL_Idiopathic_pulmonary_fibrosis

CHP_CTL <- as.data.frame(results(DEseq.output, contrast = c("Condition", "chp", "control"),
                                 cooksCutoff = F)) %>% 
  rownames_to_column(var = "Gene.name.ID") 

IPF_CTL <- as.data.frame(results(DEseq.output, contrast = c("Condition", "ipf", "control"),
                                 cooksCutoff = F)) %>% 
  rownames_to_column(var = "Gene.name.ID") 

GSE150910_DEG_sets <- c("CHP_CTL", "IPF_CTL")

#Adding identifiers

DEG_Proc_list <- list()

for(dataset in GSE150910_DEG_sets){
  DEG_Proc <- merge(get(dataset) %>% dplyr::rename("Gene.name" = "Gene.name.ID"), Final_Annotation_List)
  write_tsv(DEG_Proc, paste0("data/02_GSE150910_", dataset, ".tsv"))
  
  DEG_Proc_list[[dataset]] <- DEG_Proc
}

benchmark_plot_RNA(DEG_Proc_list[["CHP_CTL"]], "GSE150910: CHP vs Healthy Controls")
benchmark_plot_RNA(DEG_Proc_list[["IPF_CTL"]], "GSE150910: IPF vs Healthy Controls")


# Data set GSE231693 (transcriptomics), IPF/SSc-ILD, lung -----------

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE231693

gset <- getGEO("GSE231693")
meta <- pData(gset[[1]])

meta <- meta %>%
  mutate(Condition = gsub(" ", "_", `disease state:ch1`))

count.GSE231693_raw <- read_tsv(paste("data/raw/GSE231693_count.tsv", sep = ""))

count.GSE231693 <- count.GSE231693_raw %>% 
  dplyr::rename("Gene.name.ID" = "...1") %>% 
  filter(!is.na(Gene.name.ID))

#Analyze (diff expression: logFold changes)
# Raw counts fed into DESeq2

coldata_prel <- meta %>% 
  rownames_to_column() %>%
  column_to_rownames(var = "title")
coldata_prel$Condition <- factor(coldata_prel$Condition)
coldata <- coldata_prel 

table(coldata$Condition) #Get overview of sample numbers per group

cts_prel <- count.GSE231693 %>%
  column_to_rownames(var = "Gene.name.ID") 

cts <- as.matrix(cts_prel)

cts <- cts[, rownames(coldata)]
all(rownames(coldata) %in% colnames(cts)) #compare row/column names
all(rownames(coldata) == colnames(cts)) #compare row/column name order

# Input data already normalized? Check using boxplot (looks at max 10 random samples)
Boxplot_prior(cts, "02_boxplot_raw_GSE231693")

# DESeq2 procedure
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)

print(nrow(dds))
dds <- dds[rowSums(counts(dds)) >= 1,] 
print(nrow(dds))

# PCA 
scores <- PCA_norm(dds, "Condition", "02_PCA_GSE231693")

DEseq.output <- DESeq(dds)

# Dispersion estimate QC
plotDispEsts(DEseq.output)

# Counts normalized for library size
cts.norm = as.data.frame(counts(DEseq.output, normalized = TRUE))

## Boxplot of DESeq2-normalized data (looks at max 10 random samples)
Boxplot_post(cts.norm, "02_boxplot_norm_GSE231693")

#Pairs to analyze: 
#CTL_Idiopathic_pulmonary_fibrosis, CTL_Scleroderma_associated_interstitial_lung_disease

IPF_CTL <- as.data.frame(results(DEseq.output, contrast = c("Condition", "idiopathic_pulmonary_fibrosis", "normal"))) %>% 
  rownames_to_column(var = "Gene.name.ID") 

SSc_ILD_CTL <- as.data.frame(results(DEseq.output, contrast = c("Condition", "systemic_sclerosis", "normal"))) %>% 
  rownames_to_column(var = "Gene.name.ID") 

GSE231693_DEG_sets <- c("IPF_CTL", "SSc_ILD_CTL")

#Adding identifiers

DEG_Proc_list <- list()

for(dataset in GSE231693_DEG_sets){
  DEG_Proc <- merge(get(dataset) %>% dplyr::rename("ENSG.ID" = "Gene.name.ID"), Final_Annotation_List)
  write_tsv(DEG_Proc, paste0("data/02_GSE231693_", dataset, ".tsv"))
  
  DEG_Proc_list[[dataset]] <- DEG_Proc
}


benchmark_plot_RNA(DEG_Proc_list[["IPF_CTL"]], "GSE231693: IPF vs Healthy Controls")
benchmark_plot_RNA(DEG_Proc_list[["SSc_ILD_CTL"]], "GSE231693: SSc-ILD vs Healthy Controls")

# NEW --------------------------------------------------------------------------
# GSE154926 (Sjögren’s syndrome, salivary glands) ------------------------------
## get count data --------------------------------------------------------------
count <- read.csv(gzfile("data/raw/GSE154926_Raw_gene_counts_43pSS+7HVs.csv.gz")) %>%
  column_to_rownames(var = "X")
## get series matrix -----------------------------------------------------------
gset <- getGEO("GSE154926")
meta <- pData(gset[[1]])

## deal with semitechnical reps ------------------------------------------------
# there is no patient ID 

## make groups  ----------------------------------------------------------------
meta <- meta %>% 
  mutate(group = if_else(`diagnosis:ch1` == "Healthy volunteer", "Ctl", "SS"))

# How many samples per group
table(meta$group) 

## arrange samples -------------------------------------------------------------
count <- count[, meta$title]
all(meta$title == colnames(count))
# rename samples
colnames(count) <- meta$geo_accession

## create deseq2 object + QC ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = meta,
                              design = ~ group)
dds

# Boxplot of Raw counts: Checking if counts are already normalized 
Boxplot_prior(count, "02_boxplot_raw_GSE154926")

# filter low counts 
nrow(dds)
dds <- dds[rowSums(counts(dds)) >= 1,] 
nrow(dds)

# PCA 
scores <- PCA_norm(dds, "group", "02_PCA_GSE154926")

# Boxplot of size factor normalized counts
dds <- estimateSizeFactors(dds)
count.sf_norm <- counts(dds, normalized = T)
Boxplot_post(count.sf_norm, "02_boxplot_norm_GSE154926")

# Plot dispersion estimates (WIP)
dds <- estimateDispersions(dds)
plotDispEsts(dds) # cant figure out how to save this plot

## dge with deseq --------------------------------------------------------------
dds <- DESeq(dds)
res <- as.data.frame(results(dds, contrast = c("group", "SS", "Ctl"))) %>%
  rownames_to_column(var = "ENSG.ID")

## benchmark -------------------------------------------------------------------
benchmark_plot_RNA(res, "GSE154926: SS vs Healthy Ctls")
# not enriched

res <- res %>% left_join(Final_Annotation_List)
write_tsv(res, "data/02_GSE154926_SS_Ctl.tsv")

# GSE244120, obesity (thigh muscle) --------------------------------------------
## get count data --------------------------------------------------------------
count <- read_tsv(gzfile("data/raw/GSE244120_thigh.muscle_all.gene_counts.txt.gz")) %>%
  dplyr::select(-c(entrezgene, external_gene_name, gene_biotype, external_gene_source, transcript_count, description)) %>%
  distinct() %>%
  column_to_rownames(var = "ensembl_gene_id") 

## get series matrix -----------------------------------------------------------
gset <- getGEO("GSE244120")
meta <- pData(gset[[1]])

## make groups  ----------------------------------------------------------------
meta <- meta %>% 
  dplyr::rename(group = `metabolic group:ch1`)

# How many samples per group
table(meta$group) 

## arrange samples -------------------------------------------------------------
meta <- meta %>%
  mutate(title = gsub("tissue biopsy of thigh muscle \\[(.+)\\]", "\\1", title))

count <- count[, meta$title]
all(meta$title == colnames(count))
# rename samples
colnames(count) <- meta$geo_accession

## create deseq2 object + QC ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = meta,
                              design = ~ group)
dds

# Boxplot of Raw counts: Checking if counts are already normalized 
Boxplot_prior(count, "02_boxplot_raw_GSE244120")

# filter low counts 
nrow(dds)
dds <- dds[rowSums(counts(dds)) >= 1,] 
nrow(dds)

# PCA 
scores <- PCA_norm(dds, "group", "02_PCA_GSE244120")

# Boxplot of size factor normalized counts
dds <- estimateSizeFactors(dds)
count.sf_norm <- counts(dds, normalized = T)
Boxplot_post(count.sf_norm, "02_boxplot_norm_GSE244120")

# Plot dispersion estimates (WIP)
dds <- estimateDispersions(dds)
plotDispEsts(dds) # cant figure out how to save this plot

## dge with deseq --------------------------------------------------------------
dds <- DESeq(dds)
res1 <- as.data.frame(results(dds, contrast = c("group", "MUO", "MHO"))) %>%
  rownames_to_column(var = "ENSG.ID")

res2 <- as.data.frame(results(dds, contrast = c("group", "MUO", "MHL"))) %>%
  rownames_to_column(var = "ENSG.ID")
## benchmark -------------------------------------------------------------------
benchmark_plot_RNA(res1, "GSE244120: MUO vs MHO")
benchmark_plot_RNA(res2, "GSE244120: MUO vs MHL")

res1 <- res1 %>% left_join(Final_Annotation_List)
res2 <- res2 %>% left_join(Final_Annotation_List)
write_tsv(res1, "data/02_GSE244120_OB_MHO.tsv")
write_tsv(res2, "data/02_GSE244120_OB_MHL.tsv")

# GSE244118 (obesity, abdominal fat) -------------------------------------------
## get count data --------------------------------------------------------------
count <- read_tsv(gzfile("data/raw/GSE244118_abdominal.fat_all.gene_counts.txt.gz")) %>%
  dplyr::select(-c(entrezgene, external_gene_name, gene_biotype, external_gene_source, transcript_count, description)) %>%
  distinct() %>%
  column_to_rownames(var = "ensembl_gene_id") 

## get series matrix -----------------------------------------------------------
gset <- getGEO("GSE244118")
meta <- pData(gset[[1]])

## make groups  ----------------------------------------------------------------
meta <- meta %>% 
  dplyr::rename(group = `metabolic group:ch1`)

# How many samples per group
table(meta$group) 

## arrange samples -------------------------------------------------------------
meta <- meta %>%
  mutate(title = gsub("tissue biopsy of abdominal fat \\[(.+)\\]", "\\1", title))

count <- count[, meta$title]
all(meta$title == colnames(count))
# rename samples
colnames(count) <- meta$geo_accession

## create deseq2 object + QC ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = meta,
                              design = ~ group)
dds

# Boxplot of Raw counts: Checking if counts are already normalized 
Boxplot_prior(count, "02_boxplot_raw_GSE244118")

# filter low counts 
nrow(dds)
dds <- dds[rowSums(counts(dds)) >= 1,] 
nrow(dds)

# PCA 
scores <- PCA_norm(dds, "group", "02_PCA_GSE244118")

# Boxplot of size factor normalized counts
dds <- estimateSizeFactors(dds)
count.sf_norm <- counts(dds, normalized = T)
Boxplot_post(count.sf_norm, "02_boxplot_norm_GSE244118")

# Plot dispersion estimates (WIP)
dds <- estimateDispersions(dds)
plotDispEsts(dds) # cant figure out how to save this plot

## dge with deseq --------------------------------------------------------------
dds <- DESeq(dds)
res1 <- as.data.frame(results(dds, contrast = c("group", "MUO", "MHO"))) %>%
  rownames_to_column(var = "ENSG.ID")

res2 <- as.data.frame(results(dds, contrast = c("group", "MUO", "MHL"))) %>%
  rownames_to_column(var = "ENSG.ID")
## benchmark -------------------------------------------------------------------
benchmark_plot_RNA(res1, "GSE244118: MUO vs MHO")
benchmark_plot_RNA(res2, "GSE244118: MUO vs MHL")

res1 <- res1 %>% left_join(Final_Annotation_List)
res2 <- res2 %>% left_join(Final_Annotation_List)
write_tsv(res1, "data/02_GSE244118_OB_MHO.tsv")
write_tsv(res2, "data/02_GSE244118_OB_MHL.tsv")

# GSE141295 (kidney, IgA nephropathy) ------------------------------------------
## get counts ------------------------------------------------------------------
count <- read_tsv("data/raw/PRJNA593015_counts.tsv") 
count <- count %>%
  mutate(`...1` = gsub("(.+)\\.\\d+", "\\1", `...1`)) %>%
  column_to_rownames(var = "...1")

## get series matrix -----------------------------------------------------------
meta <- read_csv("data/raw/SraRunTable .GSE141295.txt")

## make groups  ----------------------------------------------------------------
meta <- meta %>% 
  mutate(group = if_else(tissue1 == "Kidney diagnosed for immunoglobulin A nephropathy", "IgA.neph", "Ctl"))

# How many samples per group
table(meta$group) 

## arrange samples -------------------------------------------------------------
count <- count[, meta$Run]
all(meta$Run == colnames(count))
meta <- meta %>% 
  column_to_rownames(var="Run")

## create deseq2 object + QC ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = meta,
                              design = ~ group)
dds

# Boxplot of Raw counts: Checking if counts are already normalized 
Boxplot_prior(count, "02_boxplot_raw_GSE141295")

# filter low counts 
nrow(dds)
dds <- dds[rowSums(counts(dds)) >= 1,] 
nrow(dds)

# PCA 
scores <- PCA_norm(dds, "group", "02_PCA_GSE141295")

# Boxplot of size factor normalized counts
dds <- estimateSizeFactors(dds)
count.sf_norm <- counts(dds, normalized = T)
Boxplot_post(count.sf_norm, "02_boxplot_norm_GSE141295")

# Plot dispersion estimates (WIP)
dds <- estimateDispersions(dds)
plotDispEsts(dds) # cant figure out how to save this plot

## dge with deseq --------------------------------------------------------------
dds <- DESeq(dds)
res <- as.data.frame(results(dds, contrast = c("group", "IgA.neph", "Ctl"))) %>%
  rownames_to_column(var = "ENSG.ID")

## benchmark -------------------------------------------------------------------
benchmark_plot_RNA(res, "GSE141295: IgA nephropathy vs Ctl")

res <- res %>% left_join(Final_Annotation_List)
write_tsv(res, "data/02_GSE141295_IgAN_Ctl.tsv")

# GSE125583 (Alzheimer) --------------------------------------------------------
## get counts ------------------------------------------------------------------
count <- read_tsv("data/raw/PRJNA516886_counts.tsv") 
count <- count %>%
  mutate(`...1` = gsub("(.+)\\.\\d+", "\\1", `...1`)) %>%
  column_to_rownames(var = "...1")

## get series matrix -----------------------------------------------------------
meta <- read_csv("data/raw/SraRunTable.PRJNA516886.txt")

## make groups  ----------------------------------------------------------------
meta <- meta %>% 
  mutate(group = if_else(diagnosis == "Alzheimer's disease", "Alz", "Ctl"))

# How many samples per group
table(meta$group) 

## arrange samples -------------------------------------------------------------
meta <- meta %>% filter(Run %in% colnames(count))
count <- count[, meta$Run]
all(meta$Run == colnames(count))
meta <- meta %>% 
  column_to_rownames(var="Run")

## create deseq2 object + QC ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = meta,
                              design = ~ group)
dds

# Boxplot of Raw counts: Checking if counts are already normalized 
Boxplot_prior(count, "02_boxplot_raw_GSE125583")

# filter low counts 
nrow(dds)
dds <- dds[rowSums(counts(dds)) >= 1,] 
nrow(dds)

# PCA 
scores <- PCA_norm(dds, "group", "02_PCA_GSE125583")

# Boxplot of size factor normalized counts
dds <- estimateSizeFactors(dds)
count.sf_norm <- counts(dds, normalized = T)
Boxplot_post(count.sf_norm, "02_boxplot_norm_GSE125583")

# Plot dispersion estimates 
dds <- estimateDispersions(dds)
plotDispEsts(dds) # cant figure out how to save this plot

## dge with deseq --------------------------------------------------------------
dds <- DESeq(dds)
res <- as.data.frame(results(dds, contrast = c("group", "Alz", "Ctl"))) %>%
  rownames_to_column(var = "ENSG.ID")

## benchmark -------------------------------------------------------------------
benchmark_plot_RNA(res, "GSE125583: Alzheimer's vs Ctl")

res <- res %>% left_join(Final_Annotation_List)
write_tsv(res, "data/02_GSE125583_Alz_Ctl.tsv")

# GSE236758 (Allergic reaction) ------------------------------------------------
## get counts ------------------------------------------------------------------
count <- read_tsv("data/raw/PRJNA992371_counts.tsv") 
count <- count %>%
  mutate(`...1` = gsub("(.+)\\.\\d+", "\\1", `...1`)) %>%
  column_to_rownames(var = "...1")

## get series matrix -----------------------------------------------------------
meta <- read_csv("data/raw/SraRunTable.PRJNA992371.txt")

## make groups  ----------------------------------------------------------------
meta <- meta %>% 
  mutate(group = if_else(condition %in% c("Allergic contact dermatitis\\, mild", "Allergic contact dermatitis\\, severe") , "ACD", "Ctl"),
         severity = case_when(condition == "Allergic contact dermatitis\\, mild" ~ "mild",
                              condition == "Allergic contact dermatitis\\, severe" ~ "severe"))

# How many samples per group
table(meta$group) 

## arrange samples -------------------------------------------------------------
meta <- meta %>% filter(Run %in% colnames(count))
count <- count[, meta$Run]
all(meta$Run == colnames(count))
meta <- meta %>% 
  column_to_rownames(var="Run")

## create deseq2 object + QC ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = meta,
                              design = ~ subject_id + treatment)
dds

# Boxplot of Raw counts: Checking if counts are already normalized 
Boxplot_prior(count, "02_boxplot_raw_GSE236758")

# filter low counts 
nrow(dds)
dds <- dds[rowSums(counts(dds)) >= 1,] 
nrow(dds)

# PCA 
scores <- PCA_norm(dds, "treatment", "02_PCA_GSE236758")

# Boxplot of size factor normalized counts
dds <- estimateSizeFactors(dds)
count.sf_norm <- counts(dds, normalized = T)
Boxplot_post(count.sf_norm, "02_boxplot_norm_GSE236758")

# Plot dispersion estimates 
dds <- estimateDispersions(dds)
plotDispEsts(dds) 

## dge with deseq --------------------------------------------------------------
dds <- DESeq(dds)
res <- as.data.frame(results(dds, contrast = c("treatment", "PPD", "Vaseline/Petrolatum"))) %>%
  rownames_to_column(var = "ENSG.ID")

## benchmark -------------------------------------------------------------------
benchmark_plot_RNA(res, "GSE236758: PPD vs Ctl")

res <- res %>% left_join(Final_Annotation_List)
write_tsv(res, "data/02_GSE236758_Acontact_Ctl.tsv")

# GSE57148 (COPD) --------------------------------------------------------------
## get counts ------------------------------------------------------------------
count <- read_tsv("data/raw/GSE57148_counts.tsv") 
count <- count %>%
  mutate(`...1` = gsub("(.+)\\.\\d+", "\\1", `...1`)) %>%
  column_to_rownames(var = "...1")

## get series matrix -----------------------------------------------------------
meta <- read_csv("data/raw/SraRunTable.GSE57148.txt")

## make groups  ----------------------------------------------------------------
meta <- meta %>%
  mutate(group = if_else(source_name == "COPD lung tissue", "COPD", "Ctl"))

## arrange samples -------------------------------------------------------------
meta <- meta %>% filter(Run %in% colnames(count))
count <- count[, meta$Run]
all(meta$Run == colnames(count))
meta <- meta %>% 
  column_to_rownames(var="Run")

## create deseq2 object + QC ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = meta,
                              design = ~ group)
dds

dds

# Boxplot of Raw counts
Boxplot_prior(count, "02_boxplot_raw_GSE57148")

# filter low counts 
nrow(dds)
dds <- dds[rowSums(counts(dds)) >= 1,] 
nrow(dds)

# PCA 
scores <- PCA_norm(dds, "group", "02_PCA_GSE57148")

# Boxplot of size factor normalized counts
dds <- estimateSizeFactors(dds)
count.sf_norm <- counts(dds, normalized = T)
Boxplot_post(count.sf_norm, "02_boxplot_norm_GSE57148")

# Plot dispersion estimates 
dds <- estimateDispersions(dds)
plotDispEsts(dds) 

## dge with deseq --------------------------------------------------------------
dds <- DESeq(dds, minReplicatesForReplace=Inf)
res <- as.data.frame(results(dds, contrast = c("group", "COPD", "Ctl"), cooksCutoff = F)) %>%
  rownames_to_column(var = "ENSG.ID")

## benchmark -------------------------------------------------------------------
res <- res %>% left_join(Final_Annotation_List)

benchmark_plot_RNA(filter(res, Gene.type == "protein_coding"), "GSE57148: COPD vs Ctl")

write_tsv(res, "data/02_GSE57148_COPD_Ctl.tsv")

# GSE121212 reanalyzed with salmon ---------------------------------------------
## get counts ------------------------------------------------------------------
count <- read_tsv("data/raw/GSE121212_counts.tsv") 
count <- count %>%
  mutate(`...1` = gsub("(.+)\\.\\d+", "\\1", `...1`)) %>%
  column_to_rownames(var = "...1")

## get series matrix -----------------------------------------------------------
gset <- getGEO("GSE121212")
series <- pData(gset[[1]])
meta <- read_csv("data/raw/SraRunTable.GSE121212.txt")

meta <- meta %>%
  dplyr::rename(geo_accession = "GEO_Accession (exp)") %>%
  inner_join(series)

## make groups  ----------------------------------------------------------------
meta <- meta %>% 
  mutate(Skin_Type = if_else(Skin_Type == "chronic_lesion", "lesional", Skin_Type))

meta <- meta %>% 
  group_by(`patient's condition:ch1`) %>%
  mutate(patient = factor(as.integer(factor(gsub("_lesional|_non-lesional|_chronic_lesion", "", title), levels = unique(gsub("_lesional|_non-lesional", "", title))))),
         group = paste(`patient's condition:ch1`, Skin_Type, sep="_")) %>%
  ungroup()

# How many samples per group
table(meta$group) 

data.frame(meta$title, meta$patient, meta$group, meta$Skin_Type, meta$`patient's condition:ch1`)

## arrange samples -------------------------------------------------------------
meta <- meta %>% filter(Run %in% colnames(count))
count <- count[, meta$Run]
all(meta$Run == colnames(count))
meta <- meta %>% 
  column_to_rownames(var="Run")
## create deseq2 object + QC ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = meta,
                              design = ~ group)
dds

dds_paired <- DESeqDataSetFromMatrix(countData = count,
                                     colData = meta,
                                     design = ~ patient + group)
dds_paired

# Boxplot of Raw counts: Checking if counts are already normalized 
Boxplot_prior(count, "02_boxplot_raw_GSE121212_salmon")

# filter low counts 
nrow(dds)
dds <- dds[rowSums(counts(dds)) >= 1,] 
nrow(dds)

dds_paired <- dds_paired[rowSums(counts(dds_paired)) >= 1,] 

# PCA 
scores <- PCA_norm(dds, "group", "02_PCA_GSE121212_salmon")

# Boxplot of size factor normalized counts
dds <- estimateSizeFactors(dds)
count.sf_norm <- counts(dds, normalized = T)
Boxplot_post(count.sf_norm, "02_boxplot_norm_GSE121212_salmon")

# Plot dispersion estimates 
dds <- estimateDispersions(dds)
plotDispEsts(dds) # cant figure out how to save this plot

dds_paired <- estimateSizeFactors(dds_paired)
dds_paired <- estimateDispersions(dds_paired)
plotDispEsts(dds_paired)

## dge with deseq --------------------------------------------------------------
# dds <- DESeq(dds) #replacing outliers and refitting for 122 genes
dds <- DESeq(dds, minReplicatesForReplace=Inf)
res1 <- as.data.frame(results(dds, contrast = c("group", "AD_lesional", "CTRL_healthy"),
                              cooksCutoff=FALSE)) %>%
  rownames_to_column(var = "ENSG.ID")

res2 <- as.data.frame(results(dds, contrast = c("group", "PSO_lesional", "CTRL_healthy"),
                              cooksCutoff=FALSE)) %>%
  rownames_to_column(var = "ENSG.ID")

dds_paired <- nbinomWaldTest(dds_paired, maxit=1000) # 71 rows did not converge in beta
dds_paired <- dds_paired[which(mcols(dds_paired)$betaConv),]

res3 <- as.data.frame(results(dds_paired, contrast = c("group", "AD_lesional", "AD_non-lesional")))%>% 
  rownames_to_column(var = "ENSG.ID")

res4 <- as.data.frame(results(dds_paired, contrast = c("group", "PSO_lesional", "PSO_non-lesional"))) %>%
  rownames_to_column(var = "ENSG.ID")

## benchmark -------------------------------------------------------------------
res1 <- res1  %>% left_join(Final_Annotation_List)
res2 <- res2  %>% left_join(Final_Annotation_List)
res3 <- res3  %>% left_join(Final_Annotation_List)
res4 <- res4  %>% left_join(Final_Annotation_List)

benchmark_plot_RNA(filter(res1, Gene.type == "protein_coding"), "GSE121212: AD vs ctl")
benchmark_plot_RNA(filter(res2, Gene.type == "protein_coding"), "GSE121212: Pso vs ctl")
benchmark_plot_RNA(filter(res3, Gene.type == "protein_coding"), "GSE121212: AD paired")
benchmark_plot_RNA(filter(res4, Gene.type == "protein_coding"), "GSE121212: Pso paired")

write_tsv(res1, "data/02_GSE121212_AD_Ctlsalmon.tsv")
write_tsv(res2, "data/02_GSE121212_PSO_Ctlsalmon.tsv")
write_tsv(res3, "data/02_GSE121212_AD_pairedsalmon.tsv")
write_tsv(res4, "data/02_GSE121212_PSO_pairedsalmon.tsv")

# GSE104704 (Alzheimer) --------------------------------------------------------
## get counts ------------------------------------------------------------------
count <- read_tsv("data/raw/GSE104704_counts.tsv") 
count <- count %>%
  mutate(`...1` = gsub("(.+)\\.\\d+", "\\1", `...1`)) %>%
  column_to_rownames(var = "...1")

## get series matrix -----------------------------------------------------------
meta <- read_csv("data/raw/SraRunTable.GSE104704.txt")

## make groups  ----------------------------------------------------------------
meta <- meta %>% 
  mutate(group = if_else(study_group == "Aged, diseased", "Alz", study_group))

# How many samples per group
table(meta$group) 

## arrange samples -------------------------------------------------------------
meta <- meta %>% filter(Run %in% colnames(count))
count <- count[, meta$Run]
all(meta$Run == colnames(count))
meta <- meta %>% 
  column_to_rownames(var="Run")

## create deseq2 object + QC ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = meta,
                              design = ~ group)
dds

# Boxplot of Raw counts: Checking if counts are already normalized 
Boxplot_prior(count, "02_boxplot_raw_GSE104704_salmon")

# filter low counts 
nrow(dds)
dds <- dds[rowSums(counts(dds)) >= 1,] 
nrow(dds)

# PCA 
scores <- PCA_norm(dds, "group", "02_PCA_GSE104704_salmon")

# Boxplot of size factor normalized counts
dds <- estimateSizeFactors(dds)
count.sf_norm <- counts(dds, normalized = T)
Boxplot_post(count.sf_norm, "02_boxplot_norm_GSE104704_salmon")

# Plot dispersion estimates 
dds <- estimateDispersions(dds)
plotDispEsts(dds) # cant figure out how to save this plot

## dge with deseq --------------------------------------------------------------
#dds <- DESeq(dds) # replacing outliers and refitting for 355 genes
dds <- DESeq(dds, minReplicatesForReplace=Inf)
res1 <- as.data.frame(results(dds, contrast = c("group", "Alz", "Old"),
                              cooksCutoff=FALSE)) %>%
  rownames_to_column(var = "ENSG.ID")

res2 <- as.data.frame(results(dds, contrast = c("group", "Alz", "Young"),
                              cooksCutoff=FALSE)) %>%
  rownames_to_column(var = "ENSG.ID")

## benchmark -------------------------------------------------------------------
res1 <- res1  %>% left_join(Final_Annotation_List)
res2 <- res2  %>% left_join(Final_Annotation_List)

benchmark_plot_RNA(filter(res1, Gene.type == "protein_coding"), "GSE104704: Alz vs Old ctl")
benchmark_plot_RNA(filter(res2, Gene.type == "protein_coding"), "GSE104704: Alz vs Young ctl")

write_tsv(res1, "data/02_GSE104704_Alz_Old.tsv")
write_tsv(res2, "data/02_GSE104704_Alz_Young.tsv")

# GSE203206 (Alzheimer) --------------------------------------------------------
## get counts ------------------------------------------------------------------
count <- read_tsv("data/raw/GSE203206_counts.tsv") 
count <- count %>%
  mutate(`...1` = gsub("(.+)\\.\\d+", "\\1", `...1`)) %>%
  column_to_rownames(var = "...1")

## get series matrix -----------------------------------------------------------
meta <- read_csv("data/raw/SraRunTable.GSE203206.txt")

## make groups  ----------------------------------------------------------------
meta <- meta %>% 
  mutate(group = case_when(Group == "Early-onset sporadic AD" ~ "early_AD",
                           Group == "Late-onset sporadic AD" ~ "late_AD",
                           Group == "control" ~"ctl"))

# How many samples per group
table(meta$group) 

# check in PCA that early and late onset look similar

## arrange samples -------------------------------------------------------------
meta <- meta %>% filter(Run %in% colnames(count))
count <- count[, meta$Run]
all(meta$Run == colnames(count))
meta <- meta %>% 
  column_to_rownames(var="Run")

## create deseq2 object + QC ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = meta,
                              design = ~ group)
dds

# Boxplot of Raw counts: Checking if counts are already normalized 
Boxplot_prior(count, "02_boxplot_raw_GSE203206_salmon")

# filter low counts 
nrow(dds)
dds <- dds[rowSums(counts(dds)) >= 1,] 
nrow(dds)

# PCA 
scores <- PCA_norm(dds, "group", "02_PCA_GSE203206_salmon")

# Boxplot of size factor normalized counts
dds <- estimateSizeFactors(dds)
count.sf_norm <- counts(dds, normalized = T)
Boxplot_post(count.sf_norm, "02_boxplot_norm_GSE203206_salmon")

# Plot dispersion estimates 
dds <- estimateDispersions(dds)
plotDispEsts(dds) # cant figure out how to save this plot

## dge with deseq --------------------------------------------------------------
dds <- DESeq(dds) # replacing outliers and refitting for 153 genes
dds <- DESeq(dds, minReplicatesForReplace=Inf)
res1 <- as.data.frame(results(dds, contrast = c("group", "early_AD", "ctl"),
                              cooksCutoff=FALSE)) %>%
  rownames_to_column(var = "ENSG.ID")

res2 <- as.data.frame(results(dds, contrast = c("group", "late_AD", "ctl"),
                              cooksCutoff=FALSE)) %>%
  rownames_to_column(var = "ENSG.ID")

## benchmark -------------------------------------------------------------------
res1 <- res1  %>% left_join(Final_Annotation_List)
res2 <- res2  %>% left_join(Final_Annotation_List)

benchmark_plot_RNA(filter(res1, Gene.type == "protein_coding"), "GSE203206: Early Alz vs  ctl")
benchmark_plot_RNA(filter(res2, Gene.type == "protein_coding"), "GSE203206: Late Alz vs ctl")

# trying with all alz as one group
meta$group2 <- if_else(meta$group == "ctl", "ctl", "AD")
dds2 <- DESeqDataSetFromMatrix(countData = count,
                               colData = meta,
                               design = ~ group2)
dds2 <- DESeq(dds2, minReplicatesForReplace=Inf)
res <- as.data.frame(results(dds2, contrast = c("group2", "AD", "ctl"),
                             cooksCutoff=FALSE)) %>%
  rownames_to_column(var = "ENSG.ID")

res <- res  %>% left_join(Final_Annotation_List)

benchmark_plot_RNA(filter(res, Gene.type == "protein_coding"), "GSE203206: Alz vs ctl")

# best enrichment with late alz vs ctl, so keep separate
write_tsv(res1, "data/02_GSE203206_Alz_ctlearly.tsv")
write_tsv(res2, "data/02_GSE203206_Alz_ctllate.tsv")

# GSE193309 reanalyzed with salmon ---------------------------------------------
## get counts ------------------------------------------------------------------
count <- read_tsv("data/raw/GSE193309_counts.tsv") 
count <- count %>%
  mutate(`...1` = gsub("(.+)\\.\\d+", "\\1", `...1`)) %>%
  column_to_rownames(var = "...1")

## get series matrix -----------------------------------------------------------
meta <- read_csv("data/raw/SraRunTable.GSE193309.txt")

## make groups  ----------------------------------------------------------------
meta %>% dplyr::count(Subject_ID, Skin_Type, Visit_Date)
# keep only baseline visit
meta <- meta %>% 
  filter(visit_id == "01") 

print(meta %>% dplyr::count(Subject_ID, Skin_Type, visit_id, anatomic_region), n=83)
# semitechnical reps
# sum reps with collapseReplicates, need to paste LS/NL to Subject_ID

meta <- meta %>% 
  mutate(group = Skin_Type,
         disease = gsub("(CO|AD)_\\d+", "\\1", Subject_ID),
         rep_group = paste(Subject_ID, group, sep = "_")) %>%
  group_by(disease) %>%
  mutate(patient = factor(as.integer(factor(Subject_ID), 
                                     levels = unique(Subject_ID)))) %>%
  ungroup()

data.frame(meta$patient, meta$Subject_ID, meta$rep_group, meta$group, meta$disease)

print(dplyr::count(meta, patient, disease), n=56)
table(meta$group)

## arrange samples -------------------------------------------------------------
meta <- meta %>% filter(Run %in% colnames(count))
count <- count[, meta$Run]
all(meta$Run == colnames(count))

# I want to check if dispersion plot was weird bc choosing highest base number
# as opposed to summin tech reps
meta$depth <- colSums(count)
meta_depth <- meta %>% 
  group_by(rep_group) %>%
  filter(depth == max(depth)) %>%
  ungroup()

count_depth <- count[, meta_depth$Run]
all(meta_depth$Run == colnames(count_depth))

meta <- meta %>% 
  column_to_rownames(var="Run")
meta_depth <- meta_depth %>% 
  column_to_rownames(var="Run")

## create deseq2 object + QC ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = meta,
                              design = ~ group)
length(meta$rep_group)
ncol(dds)
dds <- collapseReplicates(dds, groupby = dds$rep_group)
ncol(dds)

dds_depth <- DESeqDataSetFromMatrix(countData = count_depth,
                                    colData = meta_depth,
                                    design = ~ group)

# disp plots
# filter low counts 
nrow(dds)
dds <- dds[rowSums(counts(dds)) >= 1,] 
nrow(dds)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
plotDispEsts(dds)

dds_depth <- dds_depth[rowSums(counts(dds_depth)) >= 1,] 
dds_depth <- estimateSizeFactors(dds_depth)
dds_depth <- estimateDispersions(dds_depth)
plotDispEsts(dds_depth)
# not so different, but better summing reps

# paired design
dds_paired <- DESeqDataSetFromMatrix(countData = count,
                                     colData = meta,
                                     design = ~ patient + group)

dds_paired <- collapseReplicates(dds_paired, groupby = dds_paired$rep_group)
dds_paired <- dds_paired[rowSums(counts(dds_paired)) >= 1,] 

# Boxplot of Raw counts: Checking if counts are already normalized 
Boxplot_prior(count, "02_boxplot_raw_GSE193309_salmon")

# PCA 
scores <- PCA_norm(dds, "group", "02_PCA_GSE193309_salmon")

# Boxplot of size factor normalized counts
count.sf_norm <- counts(dds, normalized = T)
Boxplot_post(count.sf_norm, "02_boxplot_norm_GSE193309_salmon")

# Plot dispersion estimates 
dds_paired <- estimateSizeFactors(dds_paired)
dds_paired <- estimateDispersions(dds_paired)
plotDispEsts(dds_paired) 

## dge with deseq --------------------------------------------------------------
# dds <- DESeq(dds) #replacing outliers and refitting for 2771 genes
dds <- DESeq(dds, minReplicatesForReplace = Inf)
res1 <- as.data.frame(results(dds, contrast = c("group", "LS", "HC"),
                              cooksCutoff=FALSE)) %>%
  rownames_to_column(var = "ENSG.ID")

# dds_paired <- DESeq(dds_paired) # 32 rows did not converge in beta
dds_paired <- nbinomWaldTest(dds_paired, maxit=1500) # 31 rows did not converge in beta
dds_paired <- dds_paired[which(mcols(dds_paired)$betaConv),]

res2 <- as.data.frame(results(dds_paired, contrast = c("group", "LS", "NL")))%>% 
  rownames_to_column(var = "ENSG.ID")

## benchmark -------------------------------------------------------------------
res1 <- res1  %>% left_join(Final_Annotation_List)
res2 <- res2  %>% left_join(Final_Annotation_List)

benchmark_plot_RNA(filter(res1, Gene.type == "protein_coding"), "GSE193309: AD vs ctl")
benchmark_plot_RNA(filter(res2, Gene.type == "protein_coding"), "GSE193309: AD paired")

write_tsv(res1, "data/02_GSE193309_AD_Ctlsalmon.tsv")
write_tsv(res2, "data/02_GSE193309_AD_pairedsalmon.tsv")


# GSE175759 (kidney, various diseases) -----------------------------------------
## get counts ------------------------------------------------------------------
count <- read_tsv("data/raw/GSE175759_counts.tsv") 
count <- count %>%
  mutate(`...1` = gsub("(.+)\\.\\d+", "\\1", `...1`)) %>%
  column_to_rownames(var = "...1")

## get series matrix -----------------------------------------------------------
meta <- read_csv("data/raw/SraRunTable.GSE175759.txt")

## make groups  ----------------------------------------------------------------
meta <- meta %>% 
  mutate(group = gsub(" ", "_", diagnosis))

# How many samples per group
table(meta$group) 

# check in PCA that early and late onset look similar

## arrange samples -------------------------------------------------------------
meta <- meta %>% filter(Run %in% colnames(count))
count <- count[, meta$Run]
all(meta$Run == colnames(count))
meta <- meta %>% 
  column_to_rownames(var="Run")

#test <- count 
#colnames(test) <- paste(meta$group, seq(1,90))
#test <- test %>% rownames_to_column("gene")


## create deseq2 object + QC ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = meta,
                              design = ~ group)
dds

# Boxplot of Raw counts: Checking if counts are already normalized 
Boxplot_prior(count, "02_boxplot_raw_GSE175759")

# filter low counts 
nrow(dds)
dds <- dds[rowSums(counts(dds)) >= 1,] 
nrow(dds)

#smallestGroupSize <- min(count(meta, group)$n)
#keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
#dds <- dds[keep,]
#nrow(dds)

# PCA 
scores <- PCA_norm(dds, "group", "02_PCA_GSE175759")

# Boxplot of size factor normalized counts
dds <- estimateSizeFactors(dds)
count.sf_norm <- counts(dds, normalized = T)

#test.norm <- count.sf_norm 
#colnames(test.norm) <- paste(meta$group, seq(1,90))
#test.norm <- data.frame(test.norm) %>% rownames_to_column("gene")

Boxplot_post(count.sf_norm, "02_boxplot_norm_GSE175759")

# Plot dispersion estimates 
dds <- estimateDispersions(dds)
plotDispEsts(dds) # cant figure out how to save this plot

## dge with deseq --------------------------------------------------------------
dds <- DESeq(dds) # replacing outliers and refitting for 426 genes --> fine bc outliers in pca

table(meta$group)

res1 <- as.data.frame(results(dds, contrast = c("group", "IgAN", "Control"))) %>%
  rownames_to_column(var = "ENSG.ID")

res2 <- as.data.frame(results(dds, contrast = c("group", "minimal_change_disease", "Control"))) %>%
  rownames_to_column(var = "ENSG.ID")

res3 <- as.data.frame(results(dds, contrast = c("group", "Diabetic_nephropathy", "Control"))) %>%
  rownames_to_column(var = "ENSG.ID")

res4 <- as.data.frame(results(dds, contrast = c("group", "FSGS", "Control"))) %>%
  rownames_to_column(var = "ENSG.ID")

res5 <- as.data.frame(results(dds, contrast = c("group", "Lupus_nephritis", "Control"))) %>%
  rownames_to_column(var = "ENSG.ID")

hist(res5$log2FoldChange)
hist(res5$stat)

res6 <- as.data.frame(results(dds, contrast = c("group", "Membranous_nephropathy", "Control"))) %>%
  rownames_to_column(var = "ENSG.ID")
## benchmark -------------------------------------------------------------------
res1 <- res1  %>% left_join(Final_Annotation_List)
res2 <- res2  %>% left_join(Final_Annotation_List)
res3 <- res3  %>% left_join(Final_Annotation_List)
res4 <- res4  %>% left_join(Final_Annotation_List)
res5 <- res5  %>% left_join(Final_Annotation_List)
res6 <- res6  %>% left_join(Final_Annotation_List)

benchmark_plot_RNA(filter(res1, Gene.type == "protein_coding"), "GSE175759: IgAN vs  ctl")
benchmark_plot_RNA(filter(res2, Gene.type == "protein_coding"), "GSE175759: MCD vs ctl")
benchmark_plot_RNA(filter(res3, Gene.type == "protein_coding"), "GSE175759: Diabetic neph vs ctl")
benchmark_plot_RNA(filter(res4, Gene.type == "protein_coding"), "GSE175759: FSGS vs ctl")
benchmark_plot_RNA(filter(res5, Gene.type == "protein_coding"), "GSE175759: Lupus vs ctl")
benchmark_plot_RNA(filter(res6, Gene.type == "protein_coding"), "GSE175759: MN vs ctl")

write_tsv(res1, "data/02_GSE175759_IgAN_ctl.tsv")
write_tsv(res2, "data/02_GSE175759_MCD_ctl.tsv")
write_tsv(res3, "data/02_GSE175759_DiabNeph_ctl.tsv")
write_tsv(res4, "data/02_GSE175759_FSGS_ctl.tsv")
write_tsv(res5, "data/02_GSE175759_LupusNeph_ctl.tsv")
write_tsv(res6, "data/02_GSE175759_MN_ctl.tsv")
