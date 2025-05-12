# Clear workspace --------------------------------------------------------------
rm(list = ls())

# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(readr)
library(ggrepel) 
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(readxl)
library(msigdbr)
library(tidytext)

# Define functions -------------------------------------------------------------
source(file = "scripts/99_project_functions.R")

# Define plot theme ------------------------------------------------------------
plotTheme <- theme_minimal() + theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5), # Smaller title size, centered
    strip.text = element_text(size = 12), # Size for facet labels
    axis.text = element_text(size = 10), # Axis tick labels
    axis.title = element_text(size = 12), 
    legend.position = "right", 
    legend.text = element_text(size = 10), # Legend text size
    legend.title = element_text(size = 12), # Bold legend title
    #panel.grid.major = element_line(color = "grey80", size = 0.5), # Subtle grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), # No minor grid lines
    plot.margin = unit(c(1, 1, 1, 1), "cm"), # Adjust margins for a clean layout
    axis.ticks = element_line(linewidth = 0.5), # Thinner axis ticks
    axis.line = element_line(linewidth = 0.5) # Subtle axis lines
  )

# Read data --------------------------------------------------------------------
all.genes <- read_tsv("data/04_rank_agg_list.tsv") 


# Prepare ranked lists ---------------------------------------------------------
## contrasts that were excluded (should borderline be included?) ---------------
AH <- read_tsv("data/02_GSE142530_AH_CTL.tsv") %>%
  filter(ENSG.ID %in% Final_Annotation_List$ENSG.ID)

igan <- read_tsv("data/02_GSE175759_IgAN_ctl.tsv") %>% 
  filter(ENSG.ID %in% Final_Annotation_List$ENSG.ID)

MS <- read_tsv("data/02_GSE138614_MS_CTL.tsv") %>%
  filter(ENSG.ID %in% Final_Annotation_List$ENSG.ID)

## Proteomics ------------------------------------------------------------------
UC.hansen.pre <- read_excel("data/supp.table.processed.xlsx", sheet = "All proteins")
# UC.andersen <- read_tsv("data/raw/DE.res.UC.Andersen.tsv")  # DE by Oana, t stat
UC.andersen <- read_tsv("data/05_DE_UC_andersen.tsv")

UC.andersen <- UC.andersen %>%
  dplyr::rename(ENSG.ID = ensembl_gene_id,
         stat = t) %>%
  filter(ENSG.ID %in% all.genes$ENSG.ID)

### UC Hansen data processing --------------------------------------------------
# filter proteins in over 70% samples (they do this in the paper)
UC.hansen.pre <- filter(UC.hansen.pre, Number.of.samples.quantified.in >= 23)

# doing this to avoid weird duplicates
multiple.genes <- filter(UC.hansen.pre, str_detect(Gene.names, ";"))
single.gene <- filter(UC.hansen.pre, !(Gene.names %in% multiple.genes$Gene.names))

multiple.genes <- multiple.genes %>%  
  mutate(Gene.name = Gene.names) %>%
  separate_rows(Gene.name, sep = ";") %>%
  filter(!(Gene.name %in% single.gene$Gene.names))

single.gene <- single.gene %>%
  mutate(Gene.name = Gene.names)

UC.hansen <- rbind(multiple.genes, single.gene) %>%
  filter(Gene.name %in% all.genes$Gene.name)  %>% # after this there are still some duplicate Gene.name, then just keep the one with lowest pval like Oana did
  mutate(pvalue = 10 ^ (-minus.log.pvalue),
         adj_pvalue = p.adjust(pvalue, method = "BH"), 
         logFC = log2(Ratio.UC.H),
         sorting.value = logFC*-log(adj_pvalue)) %>%
  drop_na(Gene.name) %>% 
  group_by(Gene.name) %>% 
  filter(pvalue == min(pvalue)) %>%
  distinct(Gene.name, .keep_all = TRUE) %>% # Keep first occurrence bc there is one weird dup with same pval 
  arrange(desc(sorting.value)) %>%
  ungroup()
# create sorting value:  logFC*-log(pval) 

UC.hansen <- UC.hansen %>%
  rename(stat = sorting.value) %>%
  inner_join(Final_Annotation_List)

## test T2D obesity ------------------------------------------------------------
t2d <- read_excel("adi7548_Data_S1.xlsx", sheet = "B. Complete_data_matrix")

t2d <- dplyr::rename(t2d, Gene.name = Genename)

length(unique(t2d$Gene.name))
count(t2d, Gene.name) %>% filter(n>1)

dim(t2d %>% filter(Gene.name %in% Final_Annotation_List$Gene.name))
t2d <- t2d %>% filter(Gene.name %in% Final_Annotation_List$Gene.name)

length(unique(t2d$Gene.name))
count(t2d, Gene.name) %>% filter(n>1)
# just one gene that is duplicated, take the one with lowest p val across contrasts

t2d <- t2d %>%
  mutate(sort.t2d.lean = as.numeric(T2D_Lean_logFC)*-log(p.adjust(as.numeric(T2D_Lean_P.Value), method = "BH")),
         sort.t2d.ob = as.numeric(T2D_Obese_logFC)*-log(p.adjust(as.numeric(T2D_Obese_P.Value), method = "BH")),
         sort.ob.lean = as.numeric(Obese_Lean_logFC)*-log(p.adjust(as.numeric(Obese_Lean_P.Value), method = "BH"))) #%>%
   #group_by(Gene.name) %>% 
   #filter(sum(as.numeric(c(T2D_Lean_P.Value,T2D_Obese_P.Value, Obese_Lean_P.Value))) == min(sum(as.numeric(c(T2D_Lean_P.Value,T2D_Obese_P.Value, Obese_Lean_P.Value))))) %>%
   #ungroup()
  
t2d <- t2d %>% filter(Gene.name != "PALM2AKAP2")

# Prepare gene sets ------------------------------------------------------------
# top 100, msigdb hallmark inflammatory response, GO:BP inflammatory response

## GO sets ---------------------------------------------------------------------
## from msigdb
go <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_name, ensembl_gene, gene_symbol, gs_exact_source) %>%
  filter(gs_exact_source %in% c("GO:0002534", "GO:0006954", "GO:0002526", "GO:0002544")) %>%
  rename(gene = ensembl_gene,
         gs = gs_name)

go %>% 
  count(gs)

length(intersect(go$gene[go$gs == "GOBP_INFLAMMATORY_RESPONSE"], filter(all.genes, Position <= 2000)$ENSG.ID))

go %>% 
  filter(gene %in% Final_Annotation_List$ENSG.ID) %>%
  count(gs)

length(unique(filter(go, gs == "GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE")$gene_symbol))

setdiff(filter(go, gs == "GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE")$gene, Final_Annotation_List$ENSG.ID)
# some duplicate gene symbols with different ensembl id because they are not in chromosome
setdiff(filter(go, gs == "GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE")$gene_symbol, Final_Annotation_List$Gene.name)
# there are also some non-coding RNAs

go %>% 
  group_by(gs) %>%
  count(gene) %>%
  filter(n>1)

go %>% 
  filter(gene %in% Final_Annotation_List$ENSG.ID) %>%
  group_by(gs) %>%
  count(gene) %>%
  filter(n>1)
# no duplicates in ENSG.ID

go <- go %>% 
  filter(gene %in% Final_Annotation_List$ENSG.ID)

## KEGG ------------------------------------------------------------------------
# CP:KEGG_LEGACY: KEGG Legacy Pathways 
# these are not update since 2011, maybe exclude?
#hsa04062 KEGG_CHEMOKINE_SIGNALING_PATHWAY
#hsa04060 KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION

## Wikipathways ----------------------------------------------------------------
# WP4493 WP_CELLS_AND_MOLECULES_INVOLVED_IN_LOCAL_ACUTE_INFLAMMATORY_RESPONSE
# WP530 WP_CYTOKINES_AND_INFLAMMATORY_RESPONSE
# WP5198 WP_INFLAMMATORY_BOWEL_DISEASE_SIGNALING ?? not found
# WP453 WP_INFLAMMATORY_RESPONSE_PATHWAY

wp <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, ensembl_gene, gene_symbol, gs_exact_source) %>%
  filter(gs_exact_source %in% c("WP4493", "WP530", "WP5198", "WP453")) %>%
  rename(gene = ensembl_gene,
         gs = gs_name)

#wp <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
#  dplyr::select(gs_name, ensembl_gene, gene_symbol, gs_exact_source) %>%
#  filter(gs_name %in% c("WP_INFLAMMATORY_BOWEL_DISEASE_SIGNALING")) %>%
#  rename(gene = ensembl_gene,
#         gs = gs_name)

wp %>% 
  count(gs)

wp %>% 
  filter(gene %in% Final_Annotation_List$ENSG.ID) %>%
  count(gs)

setdiff(unique(wp$gene), Final_Annotation_List$ENSG.ID) # again located in scaffolds or non-coding

wp <- wp %>% 
  filter(gene %in% Final_Annotation_List$ENSG.ID)

## Reactome --------------------------------------------------------------------
# R-HSA-622312 REACTOME_INFLAMMASOMES
# R-HSA-913531 REACTOME_INTERFERON_SIGNALING
# R-HSA-9020702 REACTOME_INTERLEUKIN_1_SIGNALING
# R-HSA-75893 REACTOME_TNF_SIGNALING

reactome <- msigdbr(species = "Homo sapiens", category = "C2") %>% 
  dplyr::select(gs_name, ensembl_gene, gene_symbol, gs_exact_source) %>%
  filter(gs_exact_source %in% c("R-HSA-622312", "R-HSA-913531", "R-HSA-9020702", "R-HSA-75893")) %>%
  rename(gene = ensembl_gene,
         gs = gs_name)

reactome %>% 
  count(gs)

reactome %>% 
  filter(gene %in% Final_Annotation_List$ENSG.ID) %>%
  count(gs)

setdiff(unique(reactome$gene), Final_Annotation_List$ENSG.ID) 

reactome <- reactome %>% 
  filter(gene %in% Final_Annotation_List$ENSG.ID)

## Bind all sets ---------------------------------------------------------------
top.100 = filter(all.genes, Position <= 100)$ENSG.ID
go = dplyr::select(go, gs, gene)
wp = dplyr::select(wp, gs, gene)
reactome = dplyr::select(reactome, gs, gene)
msigdb = msigdb.markers$ENSG.ID

all.sets <- rbind(
  data.frame(gs = "top100", gene = top.100) ,
  data.frame(gs = "MSigDB hallmark inflammatory response", gene = msigdb),
  go,
  wp,
  reactome
)

all.sets %>% count(gs)
all.sets %>% count(gs) %>% 
  ggplot(aes(n)) + 
  geom_histogram() +
  theme_bw() +
  labs(title = "Distribution of gene set sizes")

# write_tsv(all.sets, "data/05_inflammationGeneSets.tsv") # (06/12/2024)

all.sets <- read_tsv("data/05_inflammationGeneSets.tsv")

## Check overlap ---------------------------------------------------------------

### upset plot -----------------------------------------------------------------
#install.packages("ComplexUpset")

# Convert TERM2GENE to a binary presence/absence matrix
binary_matrix <- all.sets %>%
  pivot_wider(names_from = gs, values_from = gs, values_fn = length, values_fill = 0) %>%
  mutate(across(where(is.numeric), ~ ifelse(. > 0, 1, 0))) %>%
  column_to_rownames("gene")

as.data.frame(binary_matrix) %>% dplyr::select(starts_with("GOBP"))

# Create UpSet plot
ComplexUpset::upset(
  as.data.frame(binary_matrix) %>% dplyr::select(starts_with("GOBP")),
  intersect = colnames(as.data.frame(binary_matrix) %>% dplyr::select(starts_with("GOBP"))), # Use gene sets as intersections
  width_ratio = 0.2
)

ComplexUpset::upset(
  as.data.frame(binary_matrix),
  intersect = colnames(binary_matrix), # Use gene sets as intersections
  width_ratio = 0.2,
  min_size = 1)

ggsave("figures/05_UpSet_geneSets.png", w = 25, h=12)

jaccard_similarity <- function(x, y) {
  intersection <- length(intersect(x, y))
  union <- length(union(x, y))
  return(intersection / union)
}

# Number of gene sets
n <- length(distinct(all.sets, gs)$gs)

# Initialize a matrix to store Jaccard similarities
jaccard_matrix <- matrix(0, nrow = n, ncol = n)

row = 0
# Calculate pairwise Jaccard similarities
for (set in distinct(all.sets, gs)$gs) {
  row <- row + 1
  col = 0 
  for (set2 in distinct(all.sets, gs)$gs) {
    col <- col + 1
    jaccard_similarity_value <- jaccard_similarity(all.sets$gene[all.sets$gs == set], all.sets$gene[all.sets$gs == set2])
    jaccard_matrix[row, col] <- jaccard_similarity_value
    jaccard_matrix[col, row] <- jaccard_similarity_value # because Jaccard is symmetric
  }
}

rownames(jaccard_matrix) <- distinct(all.sets, gs)$gs
colnames(jaccard_matrix) <- distinct(all.sets, gs)$gs
# not very high
# should I calculate this for the result of each dataset? 

# checking correct
length(intersect(all.sets$gene[all.sets$gs == "GOBP_INFLAMMATORY_RESPONSE"], all.sets$gene[all.sets$gs == "top100"]))
length(union(all.sets$gene[all.sets$gs == "GOBP_INFLAMMATORY_RESPONSE"], all.sets$gene[all.sets$gs == "top100"]))

# I think overlap coefficient makes more sense in this context
overlap_coef <- function(x, y) {
  intersection <- length(intersect(x, y))
  min <- min(c(length(x), length(y)))
  return(intersection / min)
}

# Number of gene sets
n <- length(distinct(all.sets, gs)$gs)

# Initialize a matrix to store Jaccard similarities
overlap_matrix <- matrix(0, nrow = n, ncol = n)

row = 0
# Calculate pairwise Jaccard similarities
for (set in distinct(all.sets, gs)$gs) {
  row <- row + 1
  col = 0 
  for (set2 in distinct(all.sets, gs)$gs) {
    col <- col + 1
    overlap_coef_value <- overlap_coef(all.sets$gene[all.sets$gs == set], all.sets$gene[all.sets$gs == set2])
    overlap_matrix[row, col] <- overlap_coef_value
    overlap_matrix[col, row] <- overlap_coef_value 
  }
}

rownames(overlap_matrix) <- distinct(all.sets, gs)$gs
colnames(overlap_matrix) <- distinct(all.sets, gs)$gs
# the gobp smaller sets have an overlap of 1 with inflammatory response
# 

# GSEA -------------------------------------------------------------------------
run_gsea <- function(results_list, gene_sets, sorting_var){
  results_list <- results_list[order(results_list[[sorting_var]], decreasing = TRUE), ]
  ranked_gene_list <- setNames(results_list[[sorting_var]], results_list$ENSG.ID)
  
  # Run GSEA on the combined gene sets
  set.seed(100)
  gsea_results <- GSEA(
    ranked_gene_list, 
    TERM2GENE = gene_sets, 
    minGSSize = 0, 
    maxGSSize = 2000, 
    pvalueCutoff = 0.05,
    eps = 0
  )
  
  return(gsea_results)
}


## Alcoholic hepatitis ---------------------------------------------------------
gsea_AH <- run_gsea(AH, all.sets, "stat")
p_AH <- dotplot(gsea_AH, 
                x = "NES", 
                size= "GeneRatio",
                showCategory = 13) + ggtitle("Alcoholic hepatitis") 

p_AH
ggsave("figures/05_AH_dotplot.png", 
       p_AH,
       h = 9,
       w = 10)

## IgA nephropathy -------------------------------------------------------------
gsea_igan <- run_gsea(igan, all.sets, "stat")
p_igan <- dotplot(gsea_igan, x = "NES", size= "GeneRatio",showCategory = 13) + ggtitle("IgA Nephropathy")
p_igan
ggsave("figures/05_igan_dotplot.png", 
       p_igan,
       h = 9,
       w = 10)

## Multiple sclerosis ----------------------------------------------------------
gsea_MS <- run_gsea(MS, all.sets, "stat")
p_MS <- dotplot(gsea_MS, x = "NES", size= "GeneRatio",showCategory = 13) + ggtitle("Multiple sclerosis") 
p_MS
ggsave("figures/05_MS_dotplot.png", 
       p_MS,
       h = 9,
       w = 10)

## UC Hansen -------------------------------------------------------------------
UC.hansen <- UC.hansen %>% drop_na(stat)
gsea_UC.hansen <- run_gsea(UC.hansen, all.sets, "stat")
p_hansen <- dotplot(gsea_UC.hansen, x = "NES", size= "GeneRatio",showCategory = 13) + ggtitle("UC Hansen sorting value") 
p_hansen

gsea_UC.hansen.lfc <- run_gsea(UC.hansen, all.sets, "logFC")
dotplot(gsea_UC.hansen.lfc, x = "NES", size= "GeneRatio",showCategory = 13) + ggtitle("UC Hansen LFC") 
# super similar

ggsave("figures/05_UC_hansen_dotplot.png", 
       p_hansen,
       h = 9,
       w = 10)

## UC Andersen -----------------------------------------------------------------
UC.andersen <- UC.andersen %>% 
  group_by(ENSG.ID) %>% 
  filter (P.Value == min(P.Value))

gsea_UC.andersen <- run_gsea(UC.andersen, all.sets, "stat")
p_andersen <- dotplot(gsea_UC.andersen, x = "NES", size= "GeneRatio",showCategory = 13) + ggtitle("UC Andersen stat value") 
p_andersen

gsea_UC.andersen.lfc <- run_gsea(UC.andersen, all.sets, "logFC")
dotplot(gsea_UC.andersen.lfc, x = "NES", size= "GeneRatio",showCategory = 13) + ggtitle("UC Andersen LFC") 

ggsave("figures/05_UC_andersen_dotplot.png", 
       p_andersen,
       h = 9,
       w = 10)

## Obesity test ----------------------------------------------------------------
t2d <- t2d %>% left_join(Final_Annotation_List)
t2d.ob <- run_gsea(drop_na(t2d, sort.t2d.ob), all.sets, "sort.t2d.ob")
dotplot(t2d.ob, x = "NES", size= "GeneRatio",showCategory = 13) + ggtitle("T2d vs Ob") 


t2d.lean <- run_gsea(drop_na(t2d, sort.t2d.lean), all.sets, "sort.t2d.lean")
dotplot(t2d.lean, x = "NES", size= "GeneRatio",showCategory = 13) + ggtitle("T2d vs lean") 

ob.lean <- run_gsea(drop_na(t2d, sort.ob.lean), all.sets, "sort.ob.lean")
dotplot(ob.lean, x = "NES", size= "GeneRatio",showCategory = 13) + ggtitle("ob vs lean") 


## Table with all results ------------------------------------------------------
gsea_AH.res <- gsea_AH@result %>% mutate(dataset = "AH")
gsea_igan.res <- gsea_igan@result %>% mutate(dataset = "igan")
gsea_MS.res <- gsea_MS@result %>% mutate(dataset = "MS")
gsea_UC.hansen.res <- gsea_UC.hansen@result %>% mutate(dataset = "UC.hansen")
gsea_UC.andersen.res <- gsea_UC.andersen@result %>% mutate(dataset = "UC.andersen")

results <- rbind(gsea_AH.res,
                 gsea_igan.res,
                 gsea_MS.res,
                 gsea_UC.hansen.res,
                 gsea_UC.andersen.res)

# making contrast names nice for paper
results <- results %>%
  mutate(dataset_lab = case_when(dataset == "AH" ~ "AH (Massey et al.)",
                                 dataset == "igan" ~ "IgAN (Park et al.)",
                                 dataset == "MS" ~ "MS (Elkjaer et al.)",
                                 dataset == "UC.hansen" ~ "UC (Schniers et al.)",
                                 dataset == "UC.andersen" ~ "UC (Bennike et al.)"),
         Description = if_else(Description == "top100", "Inflammation signature", Description))

#write_tsv(results, "data/05_GSEA_results.tsv")

results %>%
  ggplot(aes(x = reorder(dataset_lab,desc(-log10(p.adjust))), y = reorder(Description,-log10(p.adjust)))) + 
  scale_size_continuous(range = c(1, 15)) +
  geom_point(aes(fill = NES, size = -log10(p.adjust)), shape = 21, colour = "black",alpha = 0.7) +
  geom_point(aes(x = dataset_lab,
                 y = reorder(Description,NES)),
             filter(results, p.adjust<0.05, !(dataset %in% c("AD.keller", "SS", "OA", "RA"))),
             size = .7) +
  scale_fill_gradient2(midpoint = 0, low = "blue",mid = "white", high = "red4", space = "Lab") +
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  ylab(NULL) + xlab(NULL) +
  labs(fill = "Normalized Enrichment Score (NES)",
    size = expression(-log[10]("adj. p-value"))) +
  theme_bw()
  
ggsave("figures/05_GSEA_results_summary.png",
       w = 11, h = 7)


results %>% 
  filter(NES>0) %>%
  group_by(dataset) %>%
  summarise(NES =max(NES)) %>%
  left_join(results)

results %>% 
  filter(NES>0) %>%
  group_by(dataset) %>%
  summarise(p.adjust = min(p.adjust)) %>%
  left_join(results)


# Combined dotplot
results <- results %>%
  dplyr::rowwise() %>%
  mutate(
    log_p = -log10(pvalue), # Calculate -log10(p)
    log_p_adj = -log10(p.adjust),
    enriched_genes = length(unlist(strsplit(core_enrichment, "/"))), # Count enriched genes
    geneRatio = enriched_genes / setSize # Calculate geneRatio
  ) %>%
  ungroup() #%>%
  #filter(dataset != "AD.keller")


results %>%
  #filter(p.adjust < 0.05, !(dataset %in% c("AD.keller", "SS"))) %>%
  #group_by(dataset) %>% # Group by dataset for faceting
  #slice_min(order_by = p.adjust, n = 5) %>% # Select the top 5 significant gene sets (smallest p.adjust)
  #ungroup() %>%
  filter(!(dataset %in% c("AD.keller", "SS"))) %>%
  ggplot(aes(x = NES, y = reorder_within(ID, NES, dataset))) +
  geom_point(aes(size = geneRatio, color = log_p_adj)) + # Size and color by -log10(p)
  scale_size_continuous(range = c(1, 8)) + # Adjust dot size range
  #scale_color_gradientn(
  #  colors = c("gold", "yellowgreen","forestgreen","darkgreen"), # Add intermediate colors
  #  #values = c(0, 0.33, 0.67, 1), # Map colors to gradient stops (normalized to [0,1])
  #  name = expression(-log[10](p.adj)) # Label for legend
  #) +
  scale_color_viridis_c(
    direction = -1,
    option = "inferno", # Or "viridis", "magma", etc.
    #trans = "sqrt",    # Use a square root transformation to emphasize smaller values
    name = expression(-log[10](p.adj))
  ) +
  #scale_color_gradientn(
  #  colors = c("gold", "forestgreen", "darkgreen", "black"),
  #  #values = scales::rescale(c(0, 2, 5, 10)), # Compress lower values into finer steps
  #  name = expression(-log[10](p.adj))
  #)+
  facet_wrap(~ dataset, nrow = 5, scales = "free_y") +
  scale_y_reordered() +
  labs(
    #title = "GSEA Results Across Datasets",
    x = "Normalized Enrichment Score (NES)",
    y = "Gene Set",
    size = "Gene Ratio"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"), # Facet labels
    axis.text.y = element_text(size = 10), # Pathway labels
    axis.text.x = element_text(size = 10), # NES labels
    axis.title = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

ggsave("figures/05_dotplot_combined_allsets.png",
       bg = "white",
       h=15,
       w=22)

# Volcano plots ----------------------------------------------------------------
## IgA neph --------------------------------------------------------------------
igan <- igan %>%
  arrange(desc(stat)) %>%
  mutate(idx = 1:length(stat)) %>%
  mutate(top_2000 = if_else(ENSG.ID %in% filter(all.genes, Position <= 2000)$ENSG.ID, "yes", "no"),
         top_100 = if_else(ENSG.ID %in% filter(all.genes, Position <= 100)$ENSG.ID, "yes", "no"))

ggplot(igan, aes(x = log2FoldChange, y = -log10(pvalue), color = top_2000, alpha = top_2000)) +
#ggplot(igan, aes(x = log2FoldChange, y = -log10(padj), color = top_2000, alpha = top_2000)) +
  geom_point(size = .3) +
  scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
  scale_alpha_manual(values = c("no" = 0.3, "yes" = 1)) +
  labs(
    x = "Log Fold Change",
    y = "-log10(P-Value)",
    title = "IgA nephropathy"
  ) +
  theme_minimal() +
 coord_cartesian(ylim = c(0,20))

## UC Andersen ------------------------------------------------------------------
UC.andersen <- UC.andersen %>%
  arrange(desc(stat)) %>%
  mutate(idx = 1:length(stat)) %>%
  mutate(top_2000 = if_else(ENSG.ID %in% filter(all.genes, Position <= 2000)$ENSG.ID, "yes", "no"),
         top_100 = if_else(ENSG.ID %in% filter(all.genes, Position <= 100)$ENSG.ID, "yes", "no"))

#ggplot(UC.andersen, aes(x = logFC, y = -log10(adj.P.Val), color = top_2000, alpha = top_2000)) +
  ggplot(UC.andersen, aes(x = logFC, y = -log10(P.Value), color = top_2000, alpha = top_2000)) +
  geom_point(size = .3) +
  scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
  scale_alpha_manual(values = c("no" = 0.2, "yes" = 1)) +
  labs(
    x = "Log Fold Change",
    y = "-log10(P-Value)",
    title = "UC Andersen"
  ) +
  theme_minimal()

## Multiple sclerosis ----------------------------------------------------------
MS <- MS %>%
  arrange(desc(stat)) %>%
  mutate(idx = 1:length(stat)) %>%
  mutate(top_2000 = if_else(ENSG.ID %in% filter(all.genes, Position <= 2000)$ENSG.ID, "yes", "no"),
         top_100 = if_else(ENSG.ID %in% filter(all.genes, Position <= 100)$ENSG.ID, "yes", "no"))

ggplot(MS, aes(x = log2FoldChange, y = -log10(pvalue), color = top_2000, alpha = top_2000)) +
  #ggplot(MS, aes(x = log2FoldChange, y = -log10(padj), color = top_2000, alpha = top_2000)) +
  geom_point(size = .3) +
  scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
  scale_alpha_manual(values = c("no" = 0.3, "yes" = 1)) +
  labs(
    x = "Log Fold Change",
    y = "-log10(P-Value)",
    title = "Multiple sclerosis"
  ) +
  theme_minimal()

## UC Hansen -------------------------------------------------------------------
UC.hansen<- UC.hansen %>%
  arrange(desc(stat)) %>%
  mutate(idx = 1:length(stat)) %>%
  mutate(top_2000 = if_else(ENSG.ID %in% filter(all.genes, Position <= 2000)$ENSG.ID, "yes", "no"),
         top_100 = if_else(ENSG.ID %in% filter(all.genes, Position <= 100)$ENSG.ID, "yes", "no"))

# no adj pval
# The “Significant” column contains a “+” if a protein met the selected 
# significance threshold (usually q-value). Additionally, p-values 
# (probability of type I error) and the corresponding q-values (corrected p-value) 
# are provided in the output table. https://link.springer.com/protocol/10.1007/978-1-4939-7493-1_7#Sec16 
# from perseus protocol, I assume p values here are not corrected
UC.hansen$adj_pvalue <- p.adjust(UC.hansen$pvalue, method = "BH")

UC.hansen$significance <- ifelse(
  UC.hansen$adj_pvalue < 0.01, 
  ifelse(UC.hansen$top_2000 == "yes", "significant_inflammatory", "significant"), 
  "not_significant"
)

write_tsv(UC.hansen, "data/05_UC_hansen_processed.tsv")

ggplot(UC.hansen, aes(x = logFC, y = -log10(adj_pvalue), color = significance, alpha = significance)) +
  geom_point(size = .3) +
  scale_color_manual(values = c(
    "significant" = "black", 
    "significant_inflammatory" = "red", 
    "not_significant" = "grey"
  )) +
  scale_alpha_manual(values = c("significant" = .2, 
                                "significant_inflammatory" = 1,
                                "not_significant" = .2)) +
  labs(
    x =  expression(log[2]("fold change")),
    y = expression(-log[10]("adj. p-value")),
    title = "Ulcerative colitis vs. healthy controls (Schniers et al.)"
  ) +
  geom_hline(yintercept = 2, linewidth = .3, linetype = "dashed") +
  coord_cartesian(xlim=c(-2, 2)) +
  plotTheme +
  theme(legend.position = "none")
# weird that some pvals are == 1

ggsave("figures/05_volcano_UCprot_usecase.png", h = 8, w = 10)

## Alcoholic hepatitis ---------------------------------------------------------
AH <- AH %>%
  arrange(desc(stat)) %>%
  mutate(idx = 1:length(stat)) %>%
  mutate(top_2000 = if_else(ENSG.ID %in% filter(all.genes, Position <= 2000)$ENSG.ID, "yes", "no"),
         top_100 = if_else(ENSG.ID %in% filter(all.genes, Position <= 100)$ENSG.ID, "yes", "no"))

ggplot(AH, aes(x = log2FoldChange, y = -log10(pvalue), color = top_2000, alpha = top_2000)) +
  #ggplot(AH, aes(x = log2FoldChange, y = -log10(padj), color = top_2000, alpha = top_2000)) +
  geom_point(size = .3) +
  scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
  scale_alpha_manual(values = c("no" = 0.3, "yes" = 1)) +
  labs(
    x = "Log Fold Change",
    y = "-log10(P-Value)",
    title = "Alcoholic hepatitis"
  ) +
  theme_minimal()

## Plot all --------------------------------------------------------------------
igan <- igan %>%
  mutate(case = "igan")

UC.andersen <- UC.andersen %>%
  mutate(case = "UC.andersen")

AH <- AH %>%
  mutate(case = "AH")

MS <- MS %>%
  mutate(case = "MS")

UC.hansen <- UC.hansen %>%
  mutate(case = "UC.hansen")

test <- rbind(dplyr::select(igan, ENSG.ID, stat, top_2000, top_100, case),
              dplyr::select(UC.andersen, ENSG.ID, stat,top_2000, top_100, case),
              dplyr::select(AH, ENSG.ID, stat, top_2000, top_100, case),
              dplyr::select(MS, ENSG.ID, stat, top_2000, top_100, case),
              dplyr::select(UC.hansen, ENSG.ID, stat, top_2000, top_100, case))

test %>%
  filter(case != "UC.hansen") %>%
  ggplot(aes(x=case, y = stat, color = top_2000, alpha = top_2000, size = top_2000)) +
  geom_jitter() +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c("no" = "grey", "yes" = "red")) +
  scale_alpha_manual(values = c("no" = 0.3, "yes" = 1)) +
  scale_size_manual(values = c("no" = 0.1, "yes" = .3)) +
  labs(
    x = "Dataset",
    y = "stat"
  ) +
  theme_minimal() 

test %>%
  filter(case != "UC.hansen") %>%
  ggplot(aes(x=case, y = stat, color = top_2000, alpha = top_2000, size = top_2000)) +
  geom_jitter() +
  geom_hline(yintercept = 0, linewidth = .3) +
  scale_color_manual(values = c("no" = "black", "yes" = "red")) +
  scale_alpha_manual(values = c("no" = 0.1, "yes" = 1)) +
  scale_size_manual(values = c("no" = 0.1, "yes" = .3)) +
  labs(
    x = "Dataset",
    y = "stat"
  ) +
  theme_minimal() 

test %>%
  filter(case == "UC.hansen") %>%
  ggplot(aes(x=case, y = stat, color = top_2000, alpha = top_2000, size = top_2000)) +
  geom_jitter() +
  geom_hline(yintercept = 0) +
  scale_color_manual(values = c("no" = "black", "yes" = "red")) +
  scale_alpha_manual(values = c("no" = 0.2, "yes" = 1)) +
  scale_size_manual(values = c("no" = 0.1, "yes" = .3)) +
  labs(
    x = "Dataset",
    y = expression(logFC %*% -log(pval))
  ) +
  theme_minimal() 

for(case2 in distinct(test, case)$case){
  print(case2)
  p <- test %>%
    filter(case == case2) %>%
    ggplot(aes(x=case, y = stat, color = top_2000, alpha = top_2000, size = top_2000)) +
    geom_jitter() +
    geom_hline(yintercept = 0) +
    scale_color_manual(values = c("no" = "black", "yes" = "red")) +
    scale_alpha_manual(values = c("no" = 0.1, "yes" = 1)) +
    scale_size_manual(values = c("no" = 0.1, "yes" = .4)) +
    labs(
      x = "Dataset",
      y = "stat"
    ) +
    theme_minimal() 
  
  print(p)
}


## test ------------------------------------------------------------------------
# Load packages ----------------------------------------------------------------
library(DESeq2)
library(GEOquery)

exer <- read.csv("~/Documents/inflammatome_R_project/GSE202295_gene_counts.txt.gz", sep="")

exer <- exer %>% 
  mutate(ENSG.ID = gsub("(.+)(\\.\\d+)", "\\1", Geneid)) 

rownames(exer) <- exer$ENSG.ID

count <- exer %>%
  dplyr::select(- c(gene_name, Geneid, ENSG.ID))

## get series matrix -----------------------------------------------------------
gset <- getGEO("GSE202295")
meta <- pData(gset[[1]])
gset

## clean meta  ----------------------------------------------------------------
meta <- meta %>%
  mutate(
    patient = str_extract(title, "(?<=_)(\\d+)(?=_[basal|post|rec])"),
    sample = str_extract(title, "(?<=\\[)(.*?)(?=\\])"),
    sample = paste0(sample, "Aligned.sortedByCoord.out.bam")
  ) %>%
  rename(timepoint = `timepoint:ch1`)

## arrange samples -------------------------------------------------------------
count <- count[, meta$sample]
all(meta$sample == colnames(count))
# rename samples
colnames(count) <- meta$geo_accession

## create deseq2 object + QC ---------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = meta,
                              design = ~ 1)
dds

# filter low counts 
nrow(dds)
dds <- dds[rowSums(counts(dds)) >= 1,] 
nrow(dds)

count.vst <- vst(dds, blind = FALSE)
count.vst <- assay(count.vst)

length(intersect(rownames(count.vst), top.100))
# only 25 genes
setdiff(top.100, rownames(count.vst))


# add score to meta
meta$score <- colMeans(count.vst[intersect(rownames(count.vst), top.100),])

meta %>% 
  ggplot(aes(timepoint, score)) +
  geom_boxplot()
  
meta %>% 
  ggplot(aes(timepoint, score)) +
  geom_boxplot() +
  facet_wrap(~`diagnosis:ch1`)

meta %>% 
  ggplot(aes(x = timepoint, y = score, group = patient)) +
  geom_point() +  # Add points for the scores
  geom_line() +   # Connect the dots with lines for each patient
  theme_minimal() +  # Optional: for a cleaner plot theme
  labs(x = "Timepoint", y = "Score", title = "Patient Timeline of Scores") +
  theme(legend.position = "none") + # Optional: removes the legend 
  facet_wrap(~patient)

meta %>% 
  ggplot(aes(x = timepoint, y = score, group = patient)) +
  geom_point() +  # Add points for the scores
  geom_line() +   # Connect the dots with lines for each patient
  theme_minimal() +  # Optional: for a cleaner plot theme
  labs(x = "Timepoint", y = "Score", title = "Patient Timeline of Scores") +
  theme(legend.position = "none") + # Optional: removes the legend 
  facet_wrap(~`diagnosis:ch1`)
