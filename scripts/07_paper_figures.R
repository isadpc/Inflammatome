# Clear workspace --------------------------------------------------------------
rm(list = ls())

# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(cowplot)
library(readr)

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

# 1. Benchmarking for contrast inclusion (supplementary) -----------------------
rank.contrasts <- read_tsv("data/03_sorted_contrasts.tsv")
selected.contrasts <- read_tsv("data/03_selected_contrasts.tsv")

segment.proportion <- data.frame(x1 = 0,
                                 y1 = 0,
                                 x2 = 1,
                                 y2 = 1)




# 2. Benchmarking rank agg list ------------------------------------------------
rank <- read_tsv("data/04_rank_agg_list.tsv")

## Comparison to contrasts -----------------------------------------------------
rank.p <- rank.contrasts %>%
  filter(Contrast %in% filter(selected.contrasts, include == "yes")$Contrast) %>%
  ggplot(aes(x = p_ROC, y = p_markers)) +
  geom_path(aes(group = Contrast, color = "Ranked DEG lists"), linewidth = .5, alpha = .6) +
  geom_path(aes(idx_ROC/max(idx_ROC), enr/length(union$ENSG.ID), color = "Rank aggregated list"),
            linewidth=.5, data = rank) +
  geom_segment(aes(x = x1,
                   y = y1,
                   xend = x2,
                   yend = y2), 
               data = segment.proportion, linewidth = .5) +
  labs(y = "Cumulative proportion of gold standard genes",
       x = "Cumulative proportion of other genes") +
  coord_cartesian(clip = "off") +
  scale_color_manual(values = c("#d7191c", "gray"), 
                     name = "Gene list") +
  plotTheme +
  theme(legend.position="inside", 
        legend.position.inside = c(.8, .2)) +
  guides(color = guide_legend(override.aes = list(linewidth = 3, alpha = 1)))

## Inflammatome cutoffs --------------------------------------------------------
segment.benchmark <- data.frame(x1 = 0,
                                y1 = 0,
                                x2 = dim(filter(Final_Annotation_List,
                                                Gene.type == "protein_coding"))[1],
                                y2 = dim(union)[1])

rank.zoom <- rank %>%
  filter(Position <= 3000) %>%
  ggplot(aes(Position, enr)) +
  geom_path(linewidth=.5, color = "#d7191c") +
  geom_segment(aes(x = x1,
                   y = y1,
                   xend = x2,
                   yend = y2), 
               data = segment.benchmark, linewidth = .5) +
  labs(y = "Cumulative count of gold standard genes",
       x = "Gene rank") +
  coord_cartesian(clip = "on",
                  xlim = c(0, 3000),
                  ylim = c(0, rank[rank$Position == 3000,]$enr)) +
  geom_vline(xintercept = 2000, linewidth = .5, linetype = "dashed") +
  geom_vline(xintercept = 100, linewidth = .5, linetype = "dashed") +
  plotTheme 

## Merge in one figure ---------------------------------------------------------
rank.fig <- plot_grid(rank.p, rank.zoom, labels = c('a', 'b'), label_size = 12)

ggsave("figures/07_rank_agg_benchmark.png",
       rank.fig,
       h=7,
       w=18)

# 3. GSEA results --------------------------------------------------------------


# 4. Volcano plot --------------------------------------------------------------


# 5. Correlation to scores -----------------------------------------------------