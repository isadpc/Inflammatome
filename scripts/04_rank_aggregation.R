# Clear workspace --------------------------------------------------------------
rm(list = ls())

# Load libraries ---------------------------------------------------------------
library(readr)
library(tidyverse)
library(plotly)
library(RobustRankAggreg)

# Define functions -------------------------------------------------------------
source(file = "scripts/99_project_functions.R")

# Contrast selection -----------------------------------------------------------
selected.contrasts <- read_tsv("data/03_selected_contrasts.tsv")

selected.contrasts <- filter(selected.contrasts, 
                             include == "yes")$Contrast

# Rank agreggation -------------------------------------------------------------
set_pool <- list()

for(d in selected.contrasts){
  print(d)
  set_df <- read_tsv(paste("data/", d, ".tsv", sep = "")) %>% 
    dplyr::filter(Gene.type == "protein_coding")
  try({set <- set_df %>% 
    arrange(desc(stat))})
  try({set <- set_df %>% 
    arrange(desc(t))})
  set_pool <- c(set_pool, list(set$ENSG.ID))
}

# Assign normalized ranks (between 0 and 1); all NAs receive rank "1" 
Rank_sets <- RobustRankAggreg::rankMatrix(set_pool) 

Rank_results <- BIRRA_return_all(Rank_sets) # increased number of iterations

entities <- rownames(Rank_sets)

rankedEntities <- entities[order(Rank_results$result)]

head(entities, 10)
head(rankedEntities, 10)

gene.rank <- as.data.frame(rankedEntities)
colnames(gene.rank) <- "ENSG.ID"
gene.rank <- rownames_to_column(gene.rank, var = "Position")

gene.rank <- left_join(gene.rank, Final_Annotation_List, by = "ENSG.ID")
gene.rank <- left_join(gene.rank, STRING_ENSP_to_ENSG, by = "ENSG.ID")
gene.rank$Position <- as.numeric(gene.rank$Position)

# Benchmark rank aggregated list -----------------------------------------------
# count marker occurences
gene.rank.markers <- count_markers_rank(gene.rank, union)

# Create benchmarking plot 
segment.benchmark <- data.frame(x1 = 0,
                                y1 = 0,
                                x2 = dim(filter(Final_Annotation_List, Gene.type == "protein_coding"))[1],
                                y2 = dim(union)[1])

rank_plot <- ggplot(gene.rank.markers, aes(x = Position, y = enr)) +
  geom_point(alpha = 0.3, size=1) + 
  geom_segment(aes(x = x1,
                   y = y1,
                   xend = x2,
                   yend = y2), 
               data = segment.benchmark,
               linewidth = 1, linetype = "dashed") +
  labs(title = "Gold Set Enrichment in Tissue Rank Aggregated List", 
       y = "Occurence of inflammation markers",
       x = "DEGs sorted by rank") +
  theme_classic() +
  guides(size = "none", linetype = "none", colour = "none")

rank_plot

# Zoom in and make interactive plot to see where changes in slope occur
zoom <- rank_plot +
  coord_cartesian(clip = "on",
                  xlim = c(0, 3000),
                  ylim = c(0, gene.rank.markers[gene.rank.markers$Position == 3000,]$enr))

zoom_int <- ggplotly(zoom)
zoom_int

library(htmlwidgets)
saveWidget(zoom_int, "figures/04_rank_agg_zoom_slopechanges.html")

# Save rank aggregated list ----------------------------------------------------
write_tsv(gene.rank.markers, "data/04_rank_agg_list.tsv")

# Downregulated signature ------------------------------------------------------
set_pool <- list()

for(d in selected.contrasts){
  print(d)
  set_df <- read_tsv(paste("data/", d, ".tsv", sep = "")) %>% 
    dplyr::filter(Gene.type == "protein_coding")
  try({set <- set_df %>% 
    arrange(stat)}) # arrange by increasing stat
  try({set <- set_df %>% 
    arrange(t)})
  set_pool <- c(set_pool, list(set$ENSG.ID))
}

## Assign normalized ranks (between 0 and 1); all NAs receive rank "1" 
Rank_sets <- RobustRankAggreg::rankMatrix(set_pool) 

Rank_results <- BIRRA_return_all(Rank_sets) # increased number of iterations

entities <- rownames(Rank_sets)

rankedEntities <- entities[order(Rank_results$result)]

downreg.rank <- as.data.frame(rankedEntities)
colnames(downreg.rank) <- "ENSG.ID"
downreg.rank <- rownames_to_column(downreg.rank, var = "Position")

downreg.rank <- left_join(downreg.rank, Final_Annotation_List, by = "ENSG.ID")
downreg.rank <- left_join(downreg.rank, STRING_ENSP_to_ENSG, by = "ENSG.ID")
downreg.rank$Position <- as.numeric(downreg.rank$Position)

# count marker occurences
downreg.rank.markers <- count_markers_rank(downreg.rank, union)

ggplot(downreg.rank.markers, aes(x = Position, y = enr)) +
  geom_point(alpha = 0.3, size=1) + 
  geom_segment(aes(x = x1,
                   y = y1,
                   xend = x2,
                   yend = y2), 
               data = segment.benchmark,
               linewidth = 1, linetype = "dashed") +
  labs(title = "Gold Set Enrichment in Tissue Rank Aggregated List", 
       y = "Occurence of inflammation markers",
       x = "DEGs sorted by rank") +
  theme_classic() +
  guides(size = "none", linetype = "none", colour = "none")

# save data
write_tsv(downreg.rank.markers, "data/04_rank_agg_downregulated.tsv")
