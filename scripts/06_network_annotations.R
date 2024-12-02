# Clear workspace --------------------------------------------------------------
rm(list = ls())

# Load packages ----------------------------------------------------------------
library(tidyverse)
library(readr)
library(biomaRt)

# Load annotations and top 100 markers ------------------------------------------

Final_Annotation_List <- read_tsv("data/Final_Annotation_List_Biomart.tsv") %>%
  filter(Gene.type == "protein_coding")

## text-mined markers ----------------------------------------------------------
Inflamm_Top100_marker <- read_tsv("data/Inflammation_markers_Top100.tsv")

length(intersect(Inflamm_Top100_marker$ENSG.ID, Final_Annotation_List$ENSG.ID))

STRING_ENSP_to_ENSG <- read_tsv( "data/human.aliases.filtered.txt",
                                 col_names = c("SPLIT", "ENSG.ID", "data.type")) %>% 
  separate("SPLIT", c("tax", "ENSP.ID")) %>%
  dplyr::select(-c(tax, data.type)) 

## msigdb inflammation markers -------------------------------------------------
library(msigdbr)

h <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, ensembl_gene, gene_symbol)

infl.markers <- h[h$gs_name == "HALLMARK_INFLAMMATORY_RESPONSE",]
msigdb.markers = infl.markers[c(2,3)]
colnames(msigdb.markers)=c("ENSG.ID","Gene.name")
# check if any are non-coding
msigdb.markers %>% 
  left_join(Final_Annotation_List) %>%
  dplyr::count(Gene.type)
# there are 22 with NA in gene.type
msigdb.markers %>% 
  left_join(Final_Annotation_List) %>%
  filter(is.na(Gene.type)) %>% print(n=22)
# if you look at gene.name they are all in annotation as protein coding,
# but the ENSG.IDs don't match

# keep only markers present in our list of protein coding genes
msigdb.markers <- msigdb.markers %>% 
  filter(ENSG.ID %in% Final_Annotation_List$ENSG.ID)

## Text-mined and msigdb union -------------------------------------------------
union <- full_join(msigdb.markers, Inflamm_Top100_marker)

## GO terms --------------------------------------------------------------------
go.markers <- read_tsv("data/markers.go.tsv")
go.markers.pos <- go.markers %>%
  filter(gene_biotype =="protein_coding", 
         GO.name != "negative regulation of inflammatory response", 
         ENSG.ID %in% Final_Annotation_List$ENSG.ID)

go.markers.neg <- go.markers %>%
  filter(gene_biotype =="protein_coding", 
         GO.name == "negative regulation of inflammatory response", 
         ENSG.ID %in% Final_Annotation_List$ENSG.ID)

go.unique <- data.frame(ENSG.ID = unique(go.markers$ENSG.ID)) %>%
  filter(ENSG.ID %in% Final_Annotation_List$ENSG.ID)

## load mouse inflammation signature -------------------------------------------
# downgrade dbplyr from 2.5.0 to avoid issues with biomart
# (devtools::install_version("dbplyr", version = "2.3.4"))

library(readxl)
musGenes <- read_excel("data/raw/msb201224-sup-0001.xls") %>%
  dplyr::rename(gene.symbol = "Gene Symbol")

# archive version of ensembl to avoid server issues
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mouse <-  useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

mus.inf <- getLDS(attributes = c("mgi_symbol"),
                 filters = "mgi_symbol", 
                 values = musGenes$gene.symbol , 
                 mart = mouse,
                 attributesL = c("hgnc_symbol", "ensembl_gene_id"), 
                 martL = human, 
                 uniqueRows=T)


#mus.inf <- getLDS(attributes = c("entrezgene_id"),
#                  filters = "entrezgene_id", 
#                  values = musGenes$`Entrez ID` , 
#                  mart = mouse,
#                  attributesL = c("hgnc_symbol", "ensembl_gene_id"), 
#                  martL = human, 
#                  uniqueRows=T)

musGenes %>% count(`Entrez ID`) %>% filter(n>1)
musGenes %>% count(Genbank) %>% filter(n>1)
sum(is.na(musGenes$Genbank))
sum(musGenes$Genbank == 0)
sum(musGenes$Genbank == "")
# seems like genbank ids are complete
# but no way to map them with biomart and apparently less reliable 
# https://support.bioconductor.org/p/117211/

sum(mus.inf$Gene.stable.ID == "")
sum(mus.inf$HGNC.symbol == "")
# 58 mouse gene symbols don't match human gene symbols, but all match
# at least 1 ensembl id

# multiple homologs in human
mus.inf %>% count(MGI.symbol) %>% filter(n>1) 

# multiple homologs in mouse
symbol.dups <- mus.inf %>% count(HGNC.symbol) %>% filter(n>1)
id.dups <- mus.inf %>% count(Gene.stable.ID) %>% filter(n>1)

# is it same genes that have dulicated symbol and ensg.id?
id.dups$Gene.stable.ID %in% mus.inf$Gene.stable.ID[mus.inf$HGNC.symbol %in% symbol.dups$HGNC.symbol]
symbol.dups$HGNC.symbol %in% mus.inf$HGNC.symbol[mus.inf$Gene.stable.ID %in% id.dups$Gene.stable.ID]

sum(is.na(Final_Annotation_List$Gene.name)) # in annotation we have 606 genes
# without gene symbol
sum(is.na(Final_Annotation_List$ENSG.ID)) # but ensgs are complete


length(unique(musGenes$gene.symbol)) 
musGenes %>% count(gene.symbol) %>% filter(n>1)
length(unique(mus.inf$MGI.symbol))
# 236 genes are not mapped, maybe no homolog in human?

mus.inf %>% count(MGI.symbol) %>% filter(n>1)
mus.inf %>% count(Gene.stable.ID) %>% filter(n>1)

length(unique(mus.inf$Gene.stable.ID))
# not all are unique
mus.inf <- unique(mus.inf$Gene.stable.ID)

length(intersect(mus.inf, Final_Annotation_List$ENSG.ID))

# save mouse inflammatome
mus.inf.df <- data.frame("ENSG.ID" = intersect(mus.inf, Final_Annotation_List$ENSG.ID))
write_tsv(mus.inf.df, "data/tmp/mouse_inflammatome.tsv")

## load inflammatome -----------------------------------------------------------
inf <- read_tsv("data/04_rank_agg_list.tsv")


# overlap with mouse inflammatome
length(intersect(inf$ENSG.ID[1:2000], mus.inf))
length(intersect(inf$ENSG.ID[1:100], mus.inf))

# overlap with union markers
length(intersect(inf$ENSG.ID[1:2000], union$ENSG.ID))
length(intersect(inf$ENSG.ID[1:100], union$ENSG.ID))

length(intersect(inf$ENSG.ID[1:2000], go.unique$ENSG.ID))
length(intersect(inf$ENSG.ID[1:100], go.unique$ENSG.ID))

# how many inflammatome genes have a mouse homologue?
inf.human.to.mus <- getLDS(attributes = c("ensembl_gene_id"),
                           filters = "ensembl_gene_id", 
                           values = inf$ENSG.ID[1:2000], 
                           mart = human,
                           attributesL = c("mgi_symbol"), 
                           martL = mouse, 
                           uniqueRows=T)

inf.human.to.mus %>% count(MGI.symbol) %>% filter(n>1)
inf.human.to.mus %>% count(Gene.stable.ID) %>% filter(n>1)
# it's a bit confusing because multiple matches in both directions

# Create annotation dataframe --------------------------------------------------
known_inflammation <- data.frame(id = inf$ENSG.ID[1:2000])
known_inflammation <- known_inflammation %>%
  mutate(gold_std = if_else(id %in% union$ENSG.ID, "yes", "no"),
         GO = if_else(id %in% go.unique$ENSG.ID, "yes", "no"),
         high_conf = if_else(id %in% inf$ENSG.ID[1:100], "yes", "no"),
         msigdb = if_else(id %in% msigdb.markers$ENSG.ID, "yes", "no"),
         text_mined = if_else(id %in% Inflamm_Top100_marker$ENSG.ID, "yes", "no"),
         mouse_inf = if_else(id %in% mus.inf, "yes", "no")) %>%
  mutate(any = if_else(rowSums(cbind(gold_std == "yes", 
                                     GO == "yes", 
                                     high_conf == "yes", 
                                     msigdb == "yes", 
                                     text_mined == "yes",
                                     mouse_inf == "yes")) > 0, "yes", "no"))
write_tsv(known_inflammation,
          "data/tmp/network_anno_inflammation_genes.tsv")

# long format
long <- known_inflammation %>%
  pivot_longer(cols = 2:8, names_to = "set", values_to = "presence") %>%
  filter(presence == "yes")

write_tsv(long,
          "data/tmp/network_anno_inflammation_genes_long.tsv")

# Venn diagram -----------------------------------------------------------------
library(VennDiagram)

venn.plot <- venn.diagram(list(inf$ENSG.ID[1:2000],
                               go.unique$ENSG.ID,
                               mus.inf),
                          category.names = c("Human General Inflammatome","GO:BP Inflammatory Response", "Mouse Inflammatome"),
                          filename = "figures/overlap_top2000_GO_mouse.png",
                          imagetype = "png",
                          output = T,
                          fill = c("red", "blue", "green"),
                          margin = .1)

grid.draw(venn.plot)

venn.diagram(list(inf$ENSG.ID[1:100],
                  go.unique$ENSG.ID,
                  mus.inf),
             category.names = c("High-confidence Inflammatome","GO:BP Inflammatory Response", "Mouse Inflammatome"),
             filename = "figures/overlap_top100_GO_mouse.png",
             imagetype = "png",
             output = T,
             fill = c("red", "blue", "green"),
             margin = .1)



df <- read_tsv("data/")

# Cluster data -----------------------------------------------------------------
# cluster column, total number of genes, and then number of genes from each set in each cluster
nodes <- read_csv("data/node_table.csv") %>%
  rename(cluster = "__mclCluster")

nodes %>% 
  group_by(cluster) %>%
  count() 

test <- nodes %>%
  dplyr::select(cluster, `query term`) %>%
  rename(id="query term") %>%
  right_join(long)

clusters_df <- test %>%
  count(cluster, set) %>%
  pivot_wider(names_from = set, values_from = n) %>%
  right_join(nodes %>% 
               group_by(cluster) %>%
               count())

clusters_df

clusters_df <- clusters_df %>%
  mutate(GO_p = GO/n,
         gold_std_p = gold_std/n,
         high_conf_p = high_conf/n,
         mouse_inf_p = mouse_inf/n)

print(clusters_df %>% 
  filter(n>=5) %>%
  arrange(desc(GO_p)), n=48)

print(clusters_df %>% 
        filter(n>=5) %>%
        arrange(desc(high_conf_p)), n=48)

# Drug targets -----------------------------------------------------------------
# checked in pharos that GPX1 is Tbio
nodes$`target::development level`[nodes$`display name`=="GPX1"] = "Tbio"
nodes %>%
  count(`target::development level`)

nodes %>%
  filter(high_conf=="yes") %>%
  count(`target::development level`)

nodes %>%
  filter(high_conf=="yes",
         `target::development level`=="Tdark") %>%
  dplyr::select("display name")

nodes %>%
  filter(`display name` %in% c("PARP14", "GBP1", "IFITM3", "SAMHD1", "TRIM22", 
                               "ISG20", "CDC25B", "HAPLN3", "TYMP", "C1orf162", 
                               "CYTH4", "JAML", "DOCK2", "TAGAP", "IKZF3")) %>%
  dplyr::select(`display name`, `target::development level`)

# Annotation for ALD network ---------------------------------------------------
ALD <- read_csv("data/significant_dysregulated_ALD.csv")

ALD <- ALD %>%
  rename(inflammation = "ANCOVA significant inflammation",
         steatosis = "ANCOVA significant steatosis",
         fibrosis = "ANCOVA significant fibrosis")

ALD_anno <- ALD %>% dplyr::select(`Gene name`, inflammation, steatosis, fibrosis) %>%
  pivot_longer(cols = 2:4, names_to = "set", values_to = "presence") %>%
  drop_na()

write_csv(ALD_anno, "data/ALD_sig_anno.csv")
