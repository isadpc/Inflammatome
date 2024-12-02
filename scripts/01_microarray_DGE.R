# Clear workspace --------------------------------------------------------------
rm(list = ls())

# Load packages ----------------------------------------------------------------
library(limma)
library(GEOquery)
library(tidyverse)
library(umap)

# Define functions -------------------------------------------------------------
source(file = "scripts/99_project_functions.R")

# GSE66407 (CD inflammed/non, UC inflammed/non, Ctls) --------------------------
## Load series and platform data from GEO --------------------------------------
gset <- getGEO("GSE66407", GSEMatrix =TRUE, AnnotGPL=FALSE) 
if (length(gset) > 1) idx <- grep("GPL19833", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# preprocess data 
ex <- geo2r_preprocess(gset)

## Make groups -----------------------------------------------------------------
# get phenodata
pd <- pData(gset)

pd <- dplyr::rename(pd, 
                    patient="patient:ch1",
                    diagnosis="diagnosis:ch1", 
                    inflammation="inflammation:ch1", 
                    tissue="tissue:ch1")

# Drop NAs in tissue, diagnosis, inflammation and corresponding expression values 
pd <- pd %>%
  drop_na(tissue, diagnosis, inflammation)

ex <- ex[,pd$geo_accession]
colnames(ex) == pd$geo_accession # right order

pd$tissue <- as.factor(pd$tissue)
pd$patient <- as.factor(pd$patient)

# new variable containing groups
pd$group <- as.factor(paste(pd$diagnosis, pd$inflammation))
levels(pd$group) <- gsub("\\s", "_", levels(pd$group))

## Check expression distribution -------------------------------------------------
as.data.frame(ex) %>%
  pivot_longer(everything(), names_to = "geo_accession", values_to = "expr") %>%
  left_join(dplyr::select(pd, geo_accession, group, diagnosis)) %>%
  ggplot(mapping=aes(diagnosis, expr, color=diagnosis)) +
  geom_violin()

plotDensities(ex, group=pd$group, legend ="topright")
plotDensities(normalizeBetweenArrays(ex, method="quantile"), group=pd$group, legend ="topright")
# needs normalization
ex <- normalizeBetweenArrays(ex, method="quantile")

##  Average semi-tech replicates -----------------------------------------------
pd <- pd %>%
  mutate(rep = paste(tissue, patient, diagnosis, inflammation, sep = "_"))

colnames(ex) == pd$geo_accession # right order
colnames(ex) <- pd$rep
colnames(ex)

ex <- as.data.frame(avearrays(ex))

pd <- pd[!duplicated(pd$rep), ]
pd$rep == colnames(ex) # right order

## UMAP plot (dimensionality reduction) ----------------------------------------
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 7, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)

cc <- palette()
palette(c(cc,"purple","brown"))

plot(ump$layout, main="UMAP plot, nbrs=7", xlab="", ylab="", col=pd$group, pch=20, cex=1.5)
legend("topright", legend=levels(pd$group), pch=20,
       col=1:nlevels(pd$group), title="Group", pt.cex=1.5)

## Create design matrix --------------------------------------------------------
options(na.action="na.fail")
# single design matrix containing patient bc there are some paired samples, and
# same patient has more than one sample in each group 
paired_design <- model.matrix(~ 0 + pd$group +  pd$tissue + pd$patient)
colnames(paired_design) <- gsub("pd\\$group|pd\\$tissue|pd\\$|\\s",
                                "",
                                colnames(paired_design))

paired_design

unpaired_design <- model.matrix(~ 0 + pd$group +  pd$tissue)
colnames(unpaired_design) <- gsub("pd\\$group|pd\\$tissue|pd\\$|\\s",
                                  "",
                                  colnames(unpaired_design))

colSums(unpaired_design)

## Define contrasts ------------------------------------------------------------
paired_contrasts_mat <- makeContrasts(CD_CD_non = CD_yes - CD_non,
                                      UC_UC_non = UC_yes - UC_non,
                                      levels=paired_design)

unpaired_contrasts_mat <- makeContrasts(CD_Ctl = CD_yes - Control_non,
                                        UC_Ctl = UC_yes - Control_non,
                                        levels = unpaired_design)

## Fit model -------------------------------------------------------------------
paired_fit <- lmFit(ex,paired_design)
paired_fit2 <- contrasts.fit(paired_fit, contrasts = paired_contrasts_mat)
paired_fit2 <- eBayes(paired_fit2)
summary(decideTests(paired_fit2))

unpaired_fit <- lmFit(ex,unpaired_design)
unpaired_fit2 <- contrasts.fit(unpaired_fit, contrasts = unpaired_contrasts_mat)
unpaired_fit2 <- eBayes(unpaired_fit2)
summary(decideTests(unpaired_fit2))

## Extract DEG and merge with annotation ---------------------------------------
# CD inflamed vs healthy
CD_ctl <- topTable(unpaired_fit2, 
                   number=Inf, 
                   coef = "CD_Ctl")

row.names(CD_ctl) <- gsub("(\\w+)_at",
                          "\\1",
                          rownames(CD_ctl))
CD_ctl_coding <- merge_genes_metadata(CD_ctl, Final_Annotation_List, T)
CD_ctl_all <- merge_genes_metadata(CD_ctl, Final_Annotation_List, F)

# Paired CD inflamed non inflamed
CD_CD_non <- topTable(paired_fit2, 
                      number=Inf, 
                      coef = "CD_CD_non")
row.names(CD_CD_non) <- gsub("(\\w+)_at",
                             "\\1",
                             rownames(CD_CD_non))
CD_CD_non_coding <- merge_genes_metadata(CD_CD_non, Final_Annotation_List, T)
CD_CD_non_all <- merge_genes_metadata(CD_CD_non, Final_Annotation_List, F)

# UC inflamed vs healthy
UC_ctl <- topTable(unpaired_fit2, 
                   number=Inf, 
                   coef = "UC_Ctl")
row.names(UC_ctl) <- gsub("(\\w+)_at",
                          "\\1",
                          rownames(UC_ctl))
UC_ctl_coding <- merge_genes_metadata(UC_ctl, Final_Annotation_List, T)
UC_ctl_all <- merge_genes_metadata(UC_ctl, Final_Annotation_List, F)

# Paired UC inflamed vs non inflamed
UC_UC_non <- topTable(paired_fit2, 
                      number=Inf, 
                      coef = "UC_UC_non")
row.names(UC_UC_non) <- gsub("(\\w+)_at",
                             "\\1",
                             rownames(UC_UC_non))
UC_UC_non_coding <- merge_genes_metadata(UC_UC_non, Final_Annotation_List, T)
UC_UC_non_all <- merge_genes_metadata(UC_UC_non, Final_Annotation_List, F)

benchmark_plot(CD_ctl_coding, "GSE66407: CD vs Ctl")
benchmark_plot(CD_CD_non_coding, "GSE66407: CD vs CD non inflammed")
benchmark_plot(UC_ctl_coding, "GSE66407: UC vs Ctl")
benchmark_plot(UC_UC_non_coding, "GSE66407: UC vs UC non inflammed")

#write_tsv(CD_ctl_coding, "DGE_microarray_res/GSE66407_CD_ctl_coding.tsv")
#write_tsv(CD_CD_non_coding, "DGE_microarray_res/GSE66407_CD_CD_non_coding.tsv")
#write_tsv(UC_ctl_coding, "DGE_microarray_res/GSE66407_UC_ctl_coding.tsv")
#write_tsv(UC_UC_non_coding, "DGE_microarray_res/GSE66407_UC_UC_non_coding.tsv")

write_tsv(CD_ctl_all, "data/01_GSE66407_CD_ctl.tsv")
write_tsv(CD_CD_non_all, "data/01_GSE66407_CD_paired.tsv")
write_tsv(UC_ctl_all, "data/01_GSE66407_UC_ctl.tsv")
write_tsv(UC_UC_non_all, "data/01_GSE66407_UC_paired.tsv")

# GSE28619 Alcoholic hepatitis (AH, Ctl) ---------------------------------------
## Load series and platform data from GEO --------------------------------------
gset <- getGEO("GSE28619", GSEMatrix =TRUE, AnnotGPL=FALSE) # annotGPL ?
# grep string == platform name 
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# preprocess data 
ex <- geo2r_preprocess(gset)

## Make groups -----------------------------------------------------------------
# get phenodata
pd <- pData(gset)
pd <- pd %>%
  mutate(group = as.factor(if_else(characteristics_ch1=="liver sample group: : Control",
                                   "Ctl", "AH")))

## UMAP plot (dimensionality reduction) ----------------------------------------
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 9, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=9", xlab="", ylab="", col=pd$group, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(pd$group), pch=20,
       col=1:nlevels(pd$group), title="Group", pt.cex=1.5)
# looks different from geo2r

## Merge with annotation here -------------------------------------------------
ex <- ps_to_ensg(ex, pd, hgu133plus2.db)

## Check expression distribution -------------------------------------------------
as.data.frame(ex) %>%
  pivot_longer(everything(), names_to = "geo_accession", values_to = "expr") %>%
  left_join(dplyr::select(pd, geo_accession, group)) %>%
  ggplot(mapping=aes(geo_accession, expr, color=group)) +
  geom_boxplot()

plotDensities(ex, group=pd$group, legend ="topright")
plotDensities(normalizeBetweenArrays(ex, method="quantile"), group=pd$group, legend ="topright")
# no change really

## Create design matrix --------------------------------------------------------
options(na.action="na.fail")
design <- model.matrix(~ 0 + pd$group)
colnames(design)<- levels(pd$group)
colSums(design)

## Define contrasts ------------------------------------------------------------
contrasts_mat <- makeContrasts(AH_Ctl = AH - Ctl,
                               levels=design)

## Fit model -------------------------------------------------------------------
fit <- lmFit(ex, design)
fit2 <- contrasts.fit(fit, contrasts = contrasts_mat)
fit2 <- eBayes(fit2)
summary(decideTests(fit2))

## Extract DEG and merge with annotation ---------------------------------------
AH_Ctl <- topTable(fit2, 
                   number=Inf, 
                   coef = "AH_Ctl")

# Function to annotate GPL570 platform (HG-U133_Plus_2 annotation)
AH_Ctl <- merge_genes_metadata(AH_Ctl, Final_Annotation_List, F)
# all genes are coding 

count(AH_Ctl, ENSG.ID) %>% filter(n>1)

benchmark_plot(AH_Ctl, "GSE28619")

write_tsv(AH_Ctl, "data/01_GSE28619_AH_ctl.tsv")

# GSE186582 (CD inflammed/ileal margin/ post-op, Ctls) -------------------------
## Load series and platform data from GEO --------------------------------------
gset <- getGEO("GSE186582", GSEMatrix =TRUE, AnnotGPL=FALSE) 
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# preprocess data 
ex <- geo2r_preprocess(gset)

## Make groups ------------------------------------------------------------------
pheno <- pData(gset)
pheno <- dplyr::rename(pheno,
                       group = "characteristics_ch1") %>%
  mutate(patient = gsub("_M6|_M0|_MI|_CTRL|_IC",
                        "",
                        title))

pheno$group <- as.factor(pheno$group)
levels(pheno$group) <- c("Ctl", "inflammed", "ileal_margin", "postop")

pheno$patient <- as.factor(pheno$patient)

## Check expression values -----------------------------------------------------
plotDensities(ex, group=pheno$group, legend ="topright")

as.data.frame(ex) %>%
  pivot_longer(everything(), names_to = "geo_accession", values_to = "expr") %>%
  left_join(dplyr::select(pheno, geo_accession, group)) %>%
  ggplot(mapping=aes(group, expr)) +
  geom_violin()

as.data.frame(ex) %>%
  pivot_longer(everything(), names_to = "geo_accession", values_to = "expr") %>%
  left_join(dplyr::select(pheno, geo_accession, group)) %>%
  ggplot(mapping=aes(geo_accession, expr, color=group)) +
  geom_boxplot()

plotDensities(normalizeBetweenArrays(ex, method="quantile"), group=pheno$group, legend ="topright")

ex <- normalizeBetweenArrays(ex, method="quantile")

## Merge with annotation here -------------------------------------------------
ex <- ps_to_ensg(ex, pheno, hgu133plus2.db)

##  Average semi-tech replicates -----------------------------------------------
test <- pheno %>% 
  group_by(patient, group) %>%
  count()

# there are no semi technical replicates, some paired samples are missing but 
# paired design takes this into account

## UMAP plot (dimensionality reduction) ----------------------------------------
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 7, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)

cc <- palette()
palette(c(cc,"purple","brown"))

plot(ump$layout, main="UMAP plot, nbrs=7", xlab="", ylab="", col=pheno$group, pch=20, cex=1.5)
legend("topright", legend=levels(pheno$group), pch=20,
       col=1:nlevels(pheno$group), title="Group", pt.cex=1.5)


## Create design matrix --------------------------------------------------------
design_paired <- model.matrix(~ 0 + pheno$group + pheno$patient)
design_unpaired <- model.matrix(~ 0 + pheno$group)
colnames(design_paired) <- make.names(colnames(design_paired))
colSums(design_paired)
colnames(design_unpaired) <- make.names(colnames(design_unpaired))
colSums(design_unpaired)

## Define contrasts ------------------------------------------------------------
inflamed_ctl_contrast <- makeContrasts(inflamed_ctl = pheno.groupinflammed - pheno.groupCtl,
                                       levels=design_unpaired)

paired_contrasts <- makeContrasts(inflamed_margin = pheno.groupinflammed - pheno.groupileal_margin,
                                  inflamed_postop = pheno.groupinflammed - pheno.grouppostop,
                                  levels = design_paired)

## Fit model -------------------------------------------------------------------
inf_ctl_fit <- lmFit(ex, design_unpaired)
inf_ctl_fit_2 <- contrasts.fit(inf_ctl_fit, contrasts = inflamed_ctl_contrast)
inf_ctl_fit_2 <- eBayes(inf_ctl_fit_2)
summary(decideTests(inf_ctl_fit_2))

paired_fit <- lmFit(ex, design_paired)
paired_fit_2 <- contrasts.fit(paired_fit, contrasts = paired_contrasts)
paired_fit_2 <- eBayes(paired_fit_2)
summary(decideTests(paired_fit_2))

## Extract  DEG and merge with annotation --------------------------------------
# Inflammed vs control
inf_ctl <- topTable(inf_ctl_fit_2, 
                    number=Inf, 
                    coef = "inflamed_ctl")

inf_ctl_all <- merge_genes_metadata(inf_ctl, Final_Annotation_List, F)
inf_ctl_coding <- merge_genes_metadata(inf_ctl, Final_Annotation_List, T)
# all genes are coding
benchmark_plot(inf_ctl_all, "GSE186582: inflammed CD vs healthy controls")

# inflamed vs ileal margin
inf_margin <- topTable(paired_fit_2, 
                       number=Inf, 
                       coef = "inflamed_margin")

inf_margin <- merge_genes_metadata(inf_margin, Final_Annotation_List, F)
benchmark_plot(inf_margin, "GSE186582: paired inflammed CD vs ileal margin")

# inflamed vs post operation
inf_post <- topTable(paired_fit_2, 
                     number=Inf, 
                     coef = "inflamed_postop")

inf_post <- merge_genes_metadata(inf_post, Final_Annotation_List, F)
benchmark_plot(inf_post, "GSE186582: paired inflammed CD vs 6 months post-op")

write_tsv(inf_ctl_all, "data/01_GSE186582_CD_ctl.tsv")
write_tsv(inf_margin, "data/01_GSE186582_CD_marginpaired.tsv")
write_tsv(inf_post, "data/01_GSE186582_CD_postpaired.tsv")


# GSE126802 (MS) ---------------------------------------------------------------
## Load series and platform data from GEO --------------------------------------
gset <- getGEO("GSE126802", GSEMatrix =TRUE, AnnotGPL=FALSE) 
if (length(gset) > 1) idx <- grep("GPL13497", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# preprocess data 
ex <- geo2r_preprocess(gset)

## Make groups -----------------------------------------------------------------
pd <- pData(gset) 

pd <- pd %>%
  mutate(group = as.factor(`ms/control:ch1`)) 

## Check expression distribution -------------------------------------------------
as.data.frame(ex) %>%
  pivot_longer(everything(), names_to = "geo_accession", values_to = "expr") %>%
  left_join(dplyr::select(pd, geo_accession, group)) %>%
  ggplot(mapping=aes(geo_accession, expr, color=group)) +
  geom_boxplot()

plotDensities(ex, group=pd$group, legend ="topright")
plotDensities(normalizeBetweenArrays(ex, method="quantile"), group=pd$group, legend ="topright")
ex <- normalizeBetweenArrays(ex, method="quantile")

## Merge with annotation here -------------------------------------------------
ex <- ps_to_ensg(ex, pd, HsAgilentDesign026652.db)

## UMAP plot (dimensionality reduction) ----------------------------------------
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 9, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)

cc <- palette()
palette(c(cc,"purple","brown"))

plot(ump$layout, main="UMAP plot, nbrs=7", xlab="", ylab="", col=pd$group, pch=20, cex=1.5)
legend("topright", legend=levels(pd$group), pch=20,
       col=1:nlevels(pd$group), title="Group", pt.cex=1.5)

## Drop NAs in group and corresponding expression values ----------------------
pd <- pd %>%
  drop_na(group)

ex <- ex[,row.names(pd)]

## Create design matrix --------------------------------------------------------
options(na.action="na.fail")
design <- model.matrix(~ 0 + pd$group)
colnames(design) <- gsub("pd\\$group",
                         "",
                         colnames(design))
colSums(design)

## Define contrasts ------------------------------------------------------------
MS_ctl_contrast <- makeContrasts(MS_ctl = MS - CTR,
                                 levels=design)


## Fit model -------------------------------------------------------------------
MS_ctl_fit <- lmFit(ex, design)
MS_ctl_fit_2 <- contrasts.fit(MS_ctl_fit, contrasts = MS_ctl_contrast)
MS_ctl_fit_2 <- eBayes(MS_ctl_fit_2)
summary(decideTests(MS_ctl_fit_2))

# no differentially expressed genes, in the publication also they subset MS according
# to cortisol levels to find DEG 

# GSE206848 (OA) ---------------------------------------------------------------
## Load series and platform data from GEO --------------------------------------
gset <- getGEO("GSE206848", GSEMatrix =TRUE, AnnotGPL=FALSE) 
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# preprocess data 
ex <- geo2r_preprocess(gset)

## Make groups -----------------------------------------------------------------
pd <- pData(gset) 

pd <- pd %>%
  mutate(group = as.factor(gsub("\\d+", "", title))) 

## Check expression distribution -------------------------------------------------
as.data.frame(ex) %>%
  pivot_longer(everything(), names_to = "geo_accession", values_to = "expr") %>%
  left_join(dplyr::select(pd, geo_accession, group)) %>%
  ggplot(mapping=aes(geo_accession, expr, color=group)) +
  geom_boxplot()

plotDensities(ex, group=pd$group, legend ="topright")
plotDensities(normalizeBetweenArrays(ex, method="quantile"), group=pd$group, legend ="topright")
#ex <- normalizeBetweenArrays(ex, method="quantile")

## Merge with annotation here -------------------------------------------------
ex <- ps_to_ensg(ex, pd, hgu133plus2.db)

## UMAP plot (dimensionality reduction) ----------------------------------------
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 9, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)

cc <- palette()
palette(c(cc,"purple","brown"))

plot(ump$layout, main="UMAP plot, nbrs=7", xlab="", ylab="", col=pd$group, pch=20, cex=1.5)
legend("topright", legend=levels(pd$group), pch=20,
       col=1:nlevels(pd$group), title="Group", pt.cex=1.5)


## Create design matrix --------------------------------------------------------
options(na.action="na.fail")
design <- model.matrix(~ 0 + pd$group)
colnames(design) <- gsub("pd\\$group",
                         "",
                         colnames(design))
colSums(design)

## Define contrasts ------------------------------------------------------------
OA_ctl_contrast <- makeContrasts(OA_ctl = OAS - NS,
                                 levels=design)


## Fit model -------------------------------------------------------------------
OA_ctl_fit <- lmFit(ex, design)
OA_ctl_fit_2 <- contrasts.fit(OA_ctl_fit, contrasts = OA_ctl_contrast)
OA_ctl_fit_2 <- eBayes(OA_ctl_fit_2)
summary(decideTests(OA_ctl_fit_2))

## Benchmarking ----------------------------------------------------------------
# Osteoarthritis vs control
OA_ctl <- topTable(OA_ctl_fit_2, 
                   number=Inf, 
                   coef = "OA_ctl")

OA_ctl_coding <- merge_genes_metadata(OA_ctl, Final_Annotation_List, T)
# all coding 

benchmark_plot(OA_ctl_coding, "GSE206848: Osteoarthritis vs healthy controls") # mislabeling?

write_tsv(OA_ctl_coding, "data/01_GSE206848_OA_Ctl.tsv")


# GSE153007 (AD, contact AD, PSORIASIS) ----------------------------------------
## Load data -------------------------------------------------------------------
gset <- getGEO("GSE153007", GSEMatrix =TRUE, AnnotGPL=FALSE) 
if (length(gset) > 1) idx <- grep("GPL6480", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# preprocess data 
ex <- geo2r_preprocess(gset)

## Make groups -----------------------------------------------------------------
pd <- pData(gset) 

pd <- pd %>%
  mutate(lesion = if_else(`lesion:ch1` == "TRUE",
                          "LS",
                          "NL")) %>%
  mutate(group = as.factor(if_else(`diagnosis:ch1` != "Control",
                                   paste(`diagnosis:ch1`, lesion, sep = "_"),
                                   `diagnosis:ch1`))) %>%
  dplyr::rename(patient = `patient_id:ch1`)

levels(pd$group) <- c("AD_LS", "AD_NL", "CD_LS", "CD_NL","Ctl",            
                      "Ps_LS", "Ps_NL" )

## Check expression distribution -----------------------------------------------
as.data.frame(ex) %>%
  pivot_longer(everything(), names_to = "geo_accession", values_to = "expr") %>%
  left_join(dplyr::select(pd, geo_accession, group)) %>%
  ggplot(mapping=aes(geo_accession, expr, color=group)) +
  geom_boxplot()

as.data.frame(ex) %>%
  pivot_longer(everything(), names_to = "geo_accession", values_to = "expr") %>%
  left_join(dplyr::select(pd, geo_accession, group)) %>%
  ggplot(mapping=aes(geo_accession, expr, color=group)) +
  geom_violin()

plotDensities(ex, group=pd$group, legend ="topright")
plotDensities(normalizeBetweenArrays(ex, method="quantile"), group=pd$group, legend ="topright")
ex <- normalizeBetweenArrays(ex, method="quantile")

## Merge with annotation here -------------------------------------------------
ex <- ps_to_ensg(ex, pd, hgug4112a.db)

##  Average semi-tech replicates -----------------------------------------------
test <- pd %>% 
  group_by(patient, group) %>%
  count()
# no reps

## UMAP plot (dimensionality reduction) ----------------------------------------
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 9, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)

cc <- palette()
palette(c(cc,"purple","brown"))

plot(ump$layout, main="UMAP plot, nbrs=7", xlab="", ylab="", col=pd$group, pch=20, cex=1.5)
legend("topright", legend=levels(pd$grou), pch=20,
       col=1:nlevels(pd$group), title="Group", pt.cex=1.5)
# there are 5 AD_LS and 2 CD_LS that cluster with the NL and controls

## Create design matrix --------------------------------------------------------
design_paired <- model.matrix(~ 0 + pd$group + pd$patient)
design_unpaired <- model.matrix(~ 0 + pd$group)
colnames(design_paired) <- gsub("pd\\$group|pd\\$patient",
                                "",
                                colnames(design_paired))
colSums(design_paired)
colnames(design_unpaired) <- gsub("pd\\$group",
                                  "",
                                  colnames(design_unpaired))
colSums(design_unpaired)

## Define contrasts ------------------------------------------------------------
lesion_ctl_contrast <- makeContrasts(AD_ctl = AD_LS - Ctl,
                                     Ps_ctl = Ps_LS - Ctl,
                                     Contact_ctl = CD_LS - Ctl,
                                     levels=design_unpaired)

paired_contrasts <- makeContrasts(AD = AD_LS - AD_NL,
                                  Ps = Ps_LS - Ps_NL,
                                  Contact = CD_LS - CD_NL,
                                  levels = design_paired)

## Fit model -------------------------------------------------------------------
LS_ctl_fit <- lmFit(ex, design_unpaired)
LS_ctl_fit_2 <- contrasts.fit(LS_ctl_fit, contrasts = lesion_ctl_contrast)
LS_ctl_fit_2 <- eBayes(LS_ctl_fit_2)
summary(decideTests(LS_ctl_fit_2))

paired_fit <- lmFit(ex, design_paired)
paired_fit_2 <- contrasts.fit(paired_fit, contrasts = paired_contrasts)
paired_fit_2 <- eBayes(paired_fit_2)
summary(decideTests(paired_fit_2))

## Benchmarking  ---------------------------------------------------------------
# case controls
AD_ctl <- topTable(LS_ctl_fit_2, 
                   number=Inf, 
                   coef = "AD_ctl") 
Ps_ctl <- topTable(LS_ctl_fit_2, 
                   number=Inf, 
                   coef = "Ps_ctl")
Contact_ctl <- topTable(LS_ctl_fit_2, 
                        number=Inf, 
                        coef = "Contact_ctl")

AD_ctl_coding <- merge_genes_metadata(AD_ctl, Final_Annotation_List, TRUE)
# all genes are coding
Ps_ctl <- merge_genes_metadata(Ps_ctl, Final_Annotation_List, TRUE)
Contact_ctl <-merge_genes_metadata(Contact_ctl, Final_Annotation_List, TRUE)

benchmark_plot(AD_ctl_coding, "GSE153007: AD lesional vs Ctl")
benchmark_plot(Ps_ctl, "GSE153007: Ps lesional vs Ctl")
benchmark_plot(Contact_ctl, "GSE153007: Contact Dermatitis lesional vs Ctl")

write_tsv(AD_ctl_coding, "data/01_GSE153007_AD_ctl.tsv")
write_tsv(Ps_ctl, "data/01_GSE153007_PSO_ctl.tsv")
write_tsv(Contact_ctl, "data/01_GSE153007_Contact_ctl.tsv")

# Paired contrasts
AD <- topTable(paired_fit_2, 
               number=Inf, 
               coef = "AD") 
Ps <- topTable(paired_fit_2, 
               number=Inf, 
               coef = "Ps")
Contact <- topTable(paired_fit_2, 
                    number=Inf, 
                    coef = "Contact")

AD <- merge_genes_metadata(AD, Final_Annotation_List, TRUE)
Ps <- merge_genes_metadata(Ps, Final_Annotation_List, TRUE)
Contact <-merge_genes_metadata(Contact, Final_Annotation_List, TRUE)

benchmark_plot(AD, "GSE153007: AD lesional vs AD non-lesional")
benchmark_plot(Ps, "GSE153007: Ps lesional vs Ps non-lesional")
benchmark_plot(Contact, "GSE153007: Contact Dermatitis lesional vs non-lesional")

write_tsv(AD, "data/01_GSE153007_AD_paired.tsv")
write_tsv(Ps, "data/01_GSE153007_PSO_paired.tsv")
write_tsv(Contact, "data/01_GSE153007_Contact_paired.tsv")

# GSE185764 (AD, psoriasis) ----------------------------------------------------
## Load ------------------------------------------------------------------------
gset <- getGEO("GSE185764", GSEMatrix =TRUE, AnnotGPL=FALSE) 
if (length(gset) > 1) idx <- grep("GPL17586", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# preprocess data 
ex <- geo2r_preprocess(gset)

## Make groups -----------------------------------------------------------------
pd <- pData(gset) 

pd$group <- as.factor(gsub("disease: ", "", pd$characteristics_ch1.1))
levels(pd$group) <- make.names(levels(pd$group)) 
levels(pd$group)

## Check expression distribution -----------------------------------------------
as.data.frame(ex) %>%
  pivot_longer(everything(), names_to = "geo_accession", values_to = "expr") %>%
  left_join(dplyr::select(pd, geo_accession, group)) %>%
  ggplot(mapping=aes(geo_accession, expr, color=group)) +
  geom_boxplot()

as.data.frame(ex) %>%
  pivot_longer(everything(), names_to = "geo_accession", values_to = "expr") %>%
  left_join(dplyr::select(pd, geo_accession, group)) %>%
  ggplot(mapping=aes(geo_accession, expr, color=group)) +
  geom_violin()

plotDensities(ex, group=pd$group, legend ="topright")
plotDensities(normalizeBetweenArrays(ex, method="quantile"), group=pd$group, legend ="topright")
ex <- normalizeBetweenArrays(ex, method="quantile")

## UMAP plot (dimensionality reduction) ----------------------------------------
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 9, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)

cc <- palette()
palette(c(cc,"purple","brown"))

plot(ump$layout, main="UMAP plot, nbrs=7", xlab="", ylab="", col=pd$group, pch=20, cex=1.5)
legend("topright", legend=levels(pd$group), pch=20,
       col=1:nlevels(pd$group), title="Group", pt.cex=1.5)

## Merge with annotation here -------------------------------------------------
ex <- ps_to_ensg(ex, pd, hta20transcriptcluster.db)

## Create design matrix --------------------------------------------------------
options(na.action="na.fail")
design <- model.matrix(~ 0 + pd$group)
colnames(design) <- gsub("pd\\$group",
                         "",
                         colnames(design))
colSums(design)

## Define contrasts ------------------------------------------------------------
contrast <- makeContrasts(AD_ctl = Atopic.Dermatitis - Healthy.control,
                          Ps_ctl = Psoriasis - Healthy.control,
                          Induced_Ps_ctl = Dupilumab.induced.psoriasis - Healthy.control,
                          levels = design)


## Fit model -------------------------------------------------------------------
fit <- lmFit(ex, design)
fit_2 <- contrasts.fit(fit, contrasts = contrast)
fit_2 <- eBayes(fit_2)
summary(decideTests(fit_2))

# no differentially expressed genes


# GSE131282, (MS lesion, MS non-lesion, Ctls) ----------------------------------
## Load series and platform data from GEO --------------------------------------
gsm_ms <- getGEO("GSE131282", GSEMatrix =TRUE, AnnotGPL=FALSE) # annotGPL ?
# grep string == platform name 
if (length(gsm_ms) > 1) idx <- grep("GPL10558", attr(gsm_ms, "names")) else idx <- 1
gsm_ms <- gsm_ms[[idx]]

# preprocess data 
ex <- geo2r_preprocess(gsm_ms)

## Make groups -----------------------------------------------------------------
pheno <- pData(gsm_ms)
pheno <- dplyr::rename(pheno,
                       patient_id = "characteristics_ch1",
                       condition = "description",
                       source = "source_name_ch1")

pheno$condition <- as.factor(pheno$condition)
levels(pheno$condition) <- c("Ctl", "MS_nonlesion", "MS_lesion")

pheno$source <- as.factor(pheno$source)
levels(pheno$source) <- c("frontal", "parietal")

pheno <- pheno %>% 
  mutate(patient_id = as.factor(gsub("patient id: ",
                                     "",
                                     patient_id)))

## Check expression distribution -----------------------------------------------
plotDensities(ex, group=pheno$condition, legend ="topright")

as.data.frame(ex) %>%
  pivot_longer(everything(), names_to = "geo_accession", values_to = "expr") %>%
  left_join(dplyr::select(pheno, geo_accession, condition)) %>%
  ggplot(mapping=aes(condition, expr)) +
  geom_violin()

as.data.frame(ex) %>%
  pivot_longer(everything(), names_to = "geo_accession", values_to = "expr") %>%
  left_join(dplyr::select(pheno, geo_accession, condition)) %>%
  ggplot(mapping=aes(geo_accession, expr, color=condition)) +
  geom_boxplot()

plotDensities(normalizeBetweenArrays(ex, method="quantile"), group=pheno$condition, legend ="topright")
ex <- normalizeBetweenArrays(ex, method="quantile")

## UMAP plot (dimensionality reduction) ----------------------------------------
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 7, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)

cc <- palette()
palette(c(cc,"purple","brown"))

plot(ump$layout, main="UMAP plot, nbrs=7", xlab="", ylab="", col=interaction(pheno$source, pheno$condition), pch=20, cex=1.5)
legend("topright", legend=levels(interaction(pheno$source, pheno$condition)), pch=20,
       col=1:nlevels(interaction(pheno$source, pheno$condition)), title="Group", pt.cex=1.5)
# source doesnt seem to matter 

## Merge with annotation here -------------------------------------------------
# before averaging technical reps bc the I'm changing sample names
# ensg as rownames instead of ps, then after deg merge with gene metadata
ex <- ps_to_ensg(ex, pheno, illuminaHumanv4.db)

## Average semitech reps -------------------------------------------------------
pheno <- pheno %>%
  mutate(rep = paste(patient_id, condition, sep = "_"))

colnames(ex) == pheno$geo_accession # right order
colnames(ex) <- pheno$rep
colnames(ex)

ex <- as.data.frame(avearrays(ex))

pheno <- pheno[!duplicated(pheno$rep), ]
pheno$rep == colnames(ex) # right order

plotDensities(ex, group=pheno$condition, legend ="topright") # normalize again?

# try normalizing again 
# ex <- normalizeBetweenArrays(ex, method="quantile")
# benchmark looks pretty much the same

## repeat UMAP plot (dimensionality reduction) ---------------------------------
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 7, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)

cc <- palette()
palette(c(cc,"purple","brown"))

plot(ump$layout, main="UMAP plot, nbrs=7", xlab="", ylab="", col=pheno$condition, pch=20, cex=1.5)
legend("topright", legend=levels(pheno$condition), pch=20,
       col=1:nlevels(pheno$condition), title="Group", pt.cex=1.5)

## Create design matrix --------------------------------------------------------
options(na.action="na.fail")
# several patients with more than one sample in each group + some paired 
# lesion non-lesion samples
paired_design <- model.matrix(~ 0 + pheno$condition + pheno$patient_id)
colnames(paired_design) <- make.names(colnames(paired_design))
colSums(paired_design)

unpaired_design <- model.matrix(~ 0 + pheno$condition)
colnames(unpaired_design) <- make.names(colnames(unpaired_design))
colSums(unpaired_design)

## Define contrasts ------------------------------------------------------------
paired_contrasts_mat <- makeContrasts(lesion_nonlesion = pheno.conditionMS_lesion - pheno.conditionMS_nonlesion,
                                      levels=paired_design)

unpaired_contrasts_mat <- makeContrasts(ctl_lesion = pheno.conditionMS_lesion - pheno.conditionCtl,
                                        ctl_nonlesion = pheno.conditionMS_nonlesion - pheno.conditionCtl,
                                        levels = unpaired_design)

## Fit model -------------------------------------------------------------------
paired_fit <- lmFit(ex, paired_design)
paired_fit2 <- contrasts.fit(paired_fit, contrasts = paired_contrasts_mat)
paired_fit2 <- eBayes(paired_fit2)
summary(decideTests(paired_fit2))

unpaired_fit <- lmFit(ex, unpaired_design)
unpaired_fit2 <- contrasts.fit(unpaired_fit, contrasts = unpaired_contrasts_mat)
unpaired_fit2 <- eBayes(unpaired_fit2)
summary(decideTests(unpaired_fit2))

## Benchmarking ----------------------------------------------------------------

MS_paired <- topTable(paired_fit2, coef="lesion_nonlesion", Inf)
# MS_paired <- AnnotationDbi_mapping(MS_paired, illuminaHumanv4.db)

MS_paired_coding <- merge_genes_metadata(MS_paired,
                                         Final_Annotation_List,
                                         TRUE)
MS_paired_all <- merge_genes_metadata(MS_paired,
                                      Final_Annotation_List,
                                      FALSE)
# same number of genes i.e. all genes are coding in this case
benchmark_plot(MS_paired_coding, "GSE131282 MS lesion vs MS non-lesion")


MS_LS_ctl <- topTable(unpaired_fit2, coef="ctl_lesion", Inf)
MS_LS_ctl_all <- merge_genes_metadata(MS_LS_ctl,
                                      Final_Annotation_List,
                                      FALSE)
benchmark_plot(MS_LS_ctl_all, "GSE131282 MS lesion vs Ctl")

MS_NL_ctl <- topTable(unpaired_fit2, coef="ctl_nonlesion", Inf)
MS_NL_ctl_all <- merge_genes_metadata(MS_NL_ctl, 
                                      Final_Annotation_List,
                                      FALSE)
benchmark_plot(MS_NL_ctl_all, "GSE131282 MS non-lesion vs Ctl")

write_tsv(MS_paired_all, "data/01_GSE131282_MS_paired.tsv")
write_tsv(MS_LS_ctl_all, "data/01_GSE131282_MS_ctl.tsv")

