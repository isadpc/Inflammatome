# Clear workspace --------------------------------------------------------------
rm(list = ls())

# Load packages-----------------------------------------------------------------
#install.packages("ncdf4", configure.args = "--with-nc-config=/opt/homebrew/bin/nc-config")
#BiocManager::install("DEP")
list.of.packages <- c("ggplot2","dplyr","tidyr","limma","DEP","biomaRt", "ggrepel")

lapply(list.of.packages, library, character.only=TRUE)

# Read the ranked list obtained by aggregation with Birra ----------------------
ranked.list <- read.csv("data/04_rank_agg_list.tsv",sep="\t",header=TRUE)

# Preprocessing ----------------------------------------------------------------

# Read initial data
df.prot = read.table("data/raw/proteinGroups.txt", header=T, sep="\t",
                     stringsAsFactors = F, comment.char = "", quote ="")
df.prot<-df.prot[-1,]

# Remove contaminants 
df.prot = df.prot[!df.prot$Reverse=="+",]
df.prot = df.prot[!df.prot$Contaminant=="+",]

# Extract Uniprot ID and Gene name from Majority.protein.IDs column
df.prot$Majority.protein.IDs=as.character(df.prot$Majority.protein.IDs)

df.prot <- df.prot %>% 
  separate_longer_delim(cols=Majority.protein.IDs, delim = ";") %>%
  separate_wider_delim(cols=Majority.protein.IDs, names=c('sp', 'UniprotID','Gene.name.full'), 
                       delim="|", cols_remove=FALSE, too_many = "error")

#df.prot <- df.prot %>% 
#  separate_wider_delim(cols=Majority.protein.IDs, names=c('sp', 'UniprotID','Gene.name.full'), 
#                       delim="|", cols_remove=FALSE, too_many = "drop")

df.prot <- df.prot %>% 
  separate_wider_delim(Gene.name.full, names=c('Gene.name', 'Gene.description'), 
                       delim="_", cols_remove=FALSE)

#Extract Uniprot - Ensembl mapping
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
genes.ensembl <- getBM(attributes = c("uniprotswissprot", "ensembl_gene_id", "hgnc_symbol","gene_biotype"),
                       filters = "uniprotswissprot", values = df.prot$UniprotID, mart = ensembl)

#write.table(genes.ensembl, paste(datadir,"genes.ensembl.UC.Andersen.tsv",sep=""), sep = "\t", row.names = F, quote=F)
# check if needed, if not delete

genes.ensembl.in.list <- genes.ensembl[genes.ensembl$ensembl_gene_id %in% ranked.list$ENSG.ID, ]

# Remove duplicates; for each uniprot ID keep the first ensembl ID, alphabetically
genes.ensembl.in.list <- genes.ensembl.in.list[order(genes.ensembl.in.list$uniprotswissprot, genes.ensembl.in.list$ensembl_gene_id), ]
genes.ensembl.in.list.dedup <- genes.ensembl.in.list[!duplicated(genes.ensembl.in.list$uniprotswissprot), ]

df.prot.dedup=merge(df.prot, genes.ensembl.in.list.dedup[c(1,2)],by.x="UniprotID",by.y="uniprotswissprot")
length(unique(df.prot.dedup$UniprotID))



# Prepare data for DEP ---------------------------------------------------------

# Extract LFQ columns 
df.LFQ <- df.prot.dedup %>% dplyr::select(starts_with("UC_") | starts_with("Ctrl_")  )
df.LFQ[df.LFQ==0] <- NA
df.LFQ <- data.frame(lapply(df.LFQ, function(x) as.numeric(as.character(x))))
rownames(df.LFQ)=df.prot.dedup$UniprotID

# Extract sample info based on column names
coldata=as.data.frame(colnames(df.LFQ))
colnames(coldata)=c("sample")
coldata <- coldata %>%
  mutate(
    condition = ifelse(grepl("UC", sample), "UC", "Ctrl"),
    patient = sub("(_[0-9]+)$", "", sample))  # Removes the last part like "_1", "_2"
all(coldata$sample==colnames(df.LFQ))

# Check if duplicated names 
df.prot.dedup$ensembl_gene_id %>% duplicated() %>% any()
df.prot.dedup$UniprotID %>% duplicated() %>% any()
uniprot.ens = df.prot.dedup[c("UniprotID", "ensembl_gene_id")]

# Prepare input for make_se function (data_unique, LFQ columns and experimental_design)
data_unique <- make_unique(df.prot.dedup, "UniprotID", "ensembl_gene_id", delim = ";")
LFQ_columns <- grep("^(UC|Ctrl)", colnames(data_unique)) # get LFQ column numbers
data_unique[, LFQ_columns] <- lapply(data_unique[,LFQ_columns], as.numeric)
experimental_design <- coldata
colnames(experimental_design)=c("label","condition.2","condition")
experimental_design$replicate <- as.numeric(sub(".*_(\\d)$", "\\1", experimental_design$label))

# Make SummarizedExperiment object
data_se <- make_se(data_unique, LFQ_columns, experimental_design)
plot_frequency(data_se)

# Filter rows where more than 50% NA values  
data_filt <- filter_proteins(data_se, "fraction", min = 0.5)

# Normalize using vsn 
data_norm <- normalize_vsn(data_filt)

# Plots
plot_normalization(data_se, data_norm)
plot_normalization(data_norm)
plot_numbers(data_filt)
plot_missval(data_norm)
plot_detect(data_norm)

# Imputation of left-censored missing data
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
plot_imputation(data_norm, data_imp)
data.vsn.imp = as.data.frame(data_imp@assays@data)
data.vsn.imp <- data.vsn.imp[,3:62]
boxplot(data.vsn.imp,las=2,main="")

protein.matrix = data.vsn.imp

# PCA --------------------------------------------------------------------------
# Standardize the data: calculate PCA on the transposed matrix so proteins are rows and samples are columns
protein_scaled <- scale(t(protein.matrix))

# Perform PCA
pca_result <- prcomp(protein_scaled, center = TRUE, scale. = TRUE)

# Extract PCA scores for the first two components
pca_data <- as.data.frame(pca_result$x[, 1:2])
colnames(pca_data) <- c("PC1", "PC2")

# If you have metadata for each sample, merge it with the PCA data for coloring
# For example, if you want to color by 'condition' from the experimental_design data:
pca_data <- pca_data %>%
  mutate(label = colnames(protein.matrix)) %>%
  left_join(experimental_design, by = "label")

# Plot the PCA results
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 2) +
  geom_text_repel(aes(label = label), size = 1,max.overlaps =20) +
  labs(title = "PCA Plot of Protein Matrix",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_classic() +
  theme(legend.title = element_blank())

ggsave("figures/05_PCA_UC_andersen.png", device = "png", width = 6, height = 5)



# DE analysis limma ------------------------------------------------------------
condition <- as.factor(coldata$condition)
design <- model.matrix(~0 + condition)
colnames(design) <- levels(condition)
block <- as.factor(coldata$patient)  # This identifies the replicates by patient ID
corfit <- duplicateCorrelation(protein.matrix, design, block = block)
fit <- lmFit(protein.matrix,design,block=block,correlation=corfit$consensus)
head(coef(fit))
contrast.matrix <- makeContrasts(UC-Ctrl,  levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2, trend = T, robust = T)
results=decideTests(fit2)
summary(results)
res<-topTable(fit2, coef=1, n=Inf, adjust="BH")

res$UniprotID=rownames(res)
res=merge(res, uniprot.ens,by="UniprotID")
filtered.res <- res %>%
  group_by(ensembl_gene_id) %>%
  slice_min(P.Value, n = 1, with_ties = FALSE) %>%
  ungroup()

write.table(filtered.res, "data/05_DE_UC_andersen.tsv", sep = "\t", row.names = F, quote=F)

# Correlation with severity ----------------------------------------------------

metadata <- data.frame(
  Patient = c("UC_1", "UC_2", "UC_3", "UC_4", "UC_5", "UC_6", "UC_7", "UC_8", "UC_9", "UC_10"),
  Colon_Inflammation_Grade_Score = c(36, 1, 25, 1, 0, 4, 36, 4, 4, 20.25),
  F_Cal = c(1285, 889, 1343, 3600, 30, 705, 575, 163, 30, 2060)
)

protein.matrix.filtered = protein.matrix[rownames(protein.matrix) %in% filtered.res$UniprotID,]
protein.matrix.filtered = merge(res[c("UniprotID","ensembl_gene_id")],protein.matrix.filtered,by.x="UniprotID",by.y="row.names")
rownames(protein.matrix.filtered)=protein.matrix.filtered$ensembl_gene_id

protein.matrix.filtered=protein.matrix.filtered[-c(1,2)]

# Reshape protein.matrix to long format
protein_long <- protein.matrix.filtered %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "protein") %>%
  gather(key = "label", value = "value", -protein)
merged_df <- protein_long %>%
  left_join(experimental_design[, c("label", "condition")], by = "label")

# Average replicates
averaged_df <- merged_df %>%
  group_by(protein, condition) %>%
  summarise(mean_value = mean(value)) %>%
  ungroup()
averaged_wide_df <- averaged_df %>%
  spread(key = condition, value = mean_value)
protein.matrix.mean=as.data.frame(averaged_wide_df)
rownames(protein.matrix.mean)=protein.matrix.mean$protein
protein.matrix.mean=protein.matrix.mean[-1]

# Reorder such as to fit the metadata
expr.UC=protein.matrix.mean[colnames(protein.matrix.mean) %in% metadata$Patient]
expr.UC=expr.UC[,metadata$Patient]

# added by me: defining expr.UC.markers.100
ranked.list$ENSG.ID[1:100]
expr.UC.markers.100 <- expr.UC[ranked.list$ENSG.ID[1:100],] %>% drop_na()

expr.UC.100.mean <- apply(expr.UC.markers.100, 2, mean)
expr.UC.100.median <- apply(expr.UC.markers.100, 2, median)

expr.UC.100.m = expr.UC.100.mean; type="mean"
#expr.UC.100.m = expr.UC.100.median; type="median"

# also defined by me: score
score <- metadata$Colon_Inflammation_Grade_Score

data <- data.frame(expr.UC.100.m, score); correlation <- cor(expr.UC.100.m, score)

write.table(data, "data/05_score_UC_andersen.tsv", sep = "\t", row.names = T, quote=F)

ggplot(data, aes(y = expr.UC.100.m, x = score)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red", se = FALSE) + 
  #labs(x = paste(type,"inflammatome"), x = "Colon inflammation grade score", title = "") +
  labs(y = "Inflammation signature-based score", x = "Colon inflammation grade score", title = "") +
  theme_classic() +
  # Add correlation coefficient as text
  annotate("text", y = max(expr.UC.100.mean) + 0.3 , x = max(score), 
           label = paste("Pearson r =", round(correlation, 2)), 
           hjust = 1.5, vjust = 3, color = "red", size = 3.5)

ggsave("figures/05_scatter.cor.UC.Andersen.png", device = "png", width = 5, height = 3)

