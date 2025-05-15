# Clear workspace --------------------------------------------------------------
rm(list = ls())

# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(readr)
library(biomaRt)
library(readxl)
library(VennDiagram)

# Define functions -------------------------------------------------------------
source(file = "scripts/99_project_functions.R")

# Read data --------------------------------------------------------------------
all.genes <- read_tsv("data/04_rank_agg_list.tsv") 
#all.sets <- read_tsv("data/05_inflammationGeneSets.tsv")
all.sets <- read_tsv("data/06_inflammationGeneSets.tsv")
mus.inf.raw <- read_excel("data/raw/msb201224-sup-0001.xls") %>%
  dplyr::rename(gene.symbol = "Gene Symbol")

# Map mouse inflammatome to human genes ----------------------------------------

#human <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "www", version = 113)
#mouse <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = "www", version = 113)

human <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = 113)
mouse <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", version = 113)

df <- mouse@attributes

mus.inf.ensembl <- getBM(attributes = c("ensembl_gene_id" , "mgi_symbol"),
                   filters = "mgi_symbol",
                   values = mus.inf.raw$gene.symbol,
                   mart = mouse)
##### test #############################
mus.inf.ensembl <- getBM(attributes = c("ensembl_gene_id" , "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_orthology_type", "hsapiens_homolog_orthology_confidence"),
                         filters = "mgi_symbol",
                         values = mus.inf.raw$gene.symbol,
                         mart = mouse)

mus.inf.ensembl %>% count(ensembl_gene_id) %>% filter(n>1)

mus.inf.ensembl <- mus.inf.ensembl %>%
  dplyr::rename(ENSG.ID = hsapiens_homolog_ensembl_gene)

length(unique(mus.inf.ensembl$ENSG.ID))
# There are 2150 unique human homologs
length(intersect(mus.inf.ensembl$ENSG.ID, Final_Annotation_List$ENSG.ID))
dim(filter(mus.inf.ensembl, ENSG.ID %in% Final_Annotation_List$ENSG.ID))
# Of those 2148 match our annotation
length(intersect(mus.inf.ensembl$ENSG.ID, filter(all.genes, Position <= 2000)$ENSG.ID))
# 695 are part of the inflammatome 

#######################################################

length(intersect(mus.inf.ensembl$mgi_symbol, mus.inf.raw$gene.symbol))
setdiff(mus.inf.raw$gene.symbol, mus.inf.ensembl$mgi_symbol)

mus.inf.ensembl %>% count(ensembl_gene_id) %>% filter(n>1)

homologs <- getHomologs(ensembl_gene_ids = mus.inf.ensembl$ensembl_gene_id,
                        species_from = 'mmusculus',
                        species_to = 'hsapiens')

homologs <- homologs %>%
  dplyr::rename(ENSG.ID = hsapiens_homolog_ensembl_gene)

homologs %>% count(ENSG.ID) %>% filter(n>1)
homologs %>% count(ensembl_gene_id) %>% filter(n>1)

length(unique(homologs$ENSG.ID))
# There are 2150 unique human homologs
length(intersect(homologs$ENSG.ID, Final_Annotation_List$ENSG.ID))
dim(filter(homologs, ENSG.ID %in% Final_Annotation_List$ENSG.ID))
# Of those 2148 match our annotation
length(intersect(homologs$ENSG.ID, filter(all.genes, Position <= 2000)$ENSG.ID))
# 695 are part of the inflammatome 

########## tests
length(unique(homologs$ENSG.ID))
length(intersect(unique(homologs$ENSG.ID), Final_Annotation_List$ENSG.ID))


mus.inf <- intersect(homologs$ENSG.ID, Final_Annotation_List$ENSG.ID)

# Overlap with GO inf. res. and mouse inf. -------------------------------------
## mouse inflammatome ----------------------------------------------------------
# fraction of human genes that are part of inflammatome
2000/20009
# number of mouse inflammatome genes that would be part of human inf. by chance
2148*(2000/20009) # 214.7034
695/215 # 3.2 times more overlap than by chance

# GO ---------------------------------------------------------------------------
inf <- filter(all.genes, Position <= 2000)$ENSG.ID
go.inf <- filter(all.sets, gs == "GOBP_INFLAMMATORY_RESPONSE")$gene
unique(go.inf)

length(go.inf) * (2000/20009)
length(intersect(go.inf, inf)) / (length(go.inf) * (2000/20009)) 
# 4 times higher than by chance

# Venn diagram -----------------------------------------------------------------
signature <- filter(all.genes, Position <= 100)$ENSG.ID

venn.diagram(list(inf,
                  go.inf,
                  mus.inf),
             category.names = c("Human Inflammatome","GO:BP Inflammatory Response", "Mouse Inflammatome"),
             filename = "figures/06_overlap_top2000_GO_mouse.png",
             imagetype = "png",
             output = T,
             fill = c("red", "blue", "green"),
             margin = .1)

venn.diagram(list(signature,
                  go.inf,
                  mus.inf),
             category.names = c("Human Inflammation Signature","GO:BP Inflammatory Response", "Mouse Inflammatome"),
             filename = "figures/06_overlap_top100_GO_mouse.png",
             imagetype = "png",
             output = T,
             fill = c("red", "blue", "green"),
             margin = .1)

# Create network annotation dataframe ------------------------------------------
network_anno <- data.frame(id = inf)
network_anno <- network_anno %>%
  mutate(gold_std = if_else(id %in% union$ENSG.ID, "yes", "no"),
         GO = if_else(id %in% go.inf, "yes", "no"),
         high_conf = if_else(id %in% signature, "yes", "no"),
         msigdb = if_else(id %in% msigdb.markers$ENSG.ID, "yes", "no"),
         text_mined = if_else(id %in% Inflamm_Top100_marker$ENSG.ID, "yes", "no"),
         mouse_inf = if_else(id %in% mus.inf, "yes", "no")) %>%
  mutate(any = if_else(rowSums(cbind(gold_std == "yes", 
                                     GO == "yes", 
                                     high_conf == "yes", 
                                     msigdb == "yes", 
                                     text_mined == "yes",
                                     mouse_inf == "yes")) > 0, "yes", "no"),
         label_color = if_else(rowSums(cbind( 
                                             GO == "yes", 
                                             high_conf == "yes", 
                                             mouse_inf == "yes")) > 0, "yes", "no"))
write_tsv(network_anno,
          "data/06_network_anno_inflammation_genes.tsv")

# long format
long <- network_anno %>%
  pivot_longer(cols = 2:8, names_to = "set", values_to = "presence") %>%
  filter(presence == "yes")

write_tsv(long,
          "data/06_network_anno_inflammation_genes_long.tsv")

#################################################################################

# Check CDC25B interactors
cdc25b <- read_tsv("~/Downloads/string_protein_annotations.tsv")

all.genes %>%
  filter(Gene.name %in% cdc25b$`#node`, Position <= 2000) %>%
  dplyr::select(Gene.name)

# CDC25B   
# FOXM1    
# CCNB1    
# CCNA2    
# CDC25C   
# PLK1     
# CDK1     
# CDK5R1   
# MAPK11   

