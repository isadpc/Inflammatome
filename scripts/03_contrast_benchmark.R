# Clear workspace --------------------------------------------------------------
rm(list = ls())

# Load libraries ---------------------------------------------------------------
library(readr)
library(tidyverse)

# Define functions -------------------------------------------------------------
source(file = "scripts/99_project_functions.R")

# Set theme --------------------------------------------------------------------
plotTheme <-  theme_classic()+  
  theme(plot.title = element_text(size = 20),  
        strip.text = element_text(size=18),  
        axis.text = element_text(size=15),  
        axis.title = element_text(size=18),  
        legend.position="right",  
        legend.text=element_text(size=8),  
        legend.title=element_text(size=10)) 

plotTheme <- theme_minimal() + 
  theme(
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

# Set palette ------------------------------------------------------------------
pal <- colorRampPalette(c("#d7191c","#fdae61", "#80cdc1","#018571"))
scales::show_col(pal(4))
scales::show_col(pal(8))

myPalette <- c(
  "#1F77B4", # Blue
  "#FF7F0E", # Orange
  "#2CA02C", # Green
  "#D62728", # Red
  "#9467BD", # Purple
  "#8C564B", # Brown
  "#E377C2", # Pink
  "#7F7F7F", # Grey
  "#BCBD22", # Yellow-green
  "#17BECF", # Teal
  "#9B59B6", # Lavender
  "#F39C12", # Yellow
  "#3498DB", # Light Blue
  "#E74C3C", # Dark Red
  "#2ECC71"  # Light Green
)

# Benchmark --------------------------------------------------------------------

## Get filepaths ---------------------------------------------------------------
files.tissue <- list.files(path = "data",
                           pattern = "01_GSE|02_GSE",
                           full.names = T)


## Get overlapping genes  ------------------------------------------------------
# this is only necessary if you want to work with intersect genes between all
# contrasts, for example for plotting purposes
ENSG.IDs <- Final_Annotation_List$ENSG.ID

for (f in files.tissue){
  print(f)
  data <- read_tsv(f)
  ENSG.IDs <- ENSG.IDs[ENSG.IDs %in% data$ENSG.ID]
}

## Apply function to count markers to all datasets -----------------------------
ranked.df <- bind_rows(lapply(files.tissue, 
                              count_markers, 
                              gold_set = union, 
                              coding_only = T,
                              intersect_only = F,
                              intersect_ids = ENSG.IDs))

# Add end row to "pick random genes until end of annotation"
ranked.df <- ranked.df %>%
  group_by(Contrast) %>%
  group_modify(~ add_row(.x, 
                         idx = dim(filter(Final_Annotation_List, Gene.type == "protein_coding"))[1],
                         idx_ROC = dim(filter(Final_Annotation_List, Gene.type == "protein_coding"))[1],
                         marker_count = dim(union)[1],
                         Contrast = .x$Contrast[1])) 

# create proportion of genes covered columns (0 to 1)
ranked.df <- ranked.df %>%
  group_by(Contrast) %>%
  mutate(p_annotation = idx/dim(Final_Annotation_List)[1],
         p_ROC = idx_ROC/dim(Final_Annotation_List)[1],
         p_markers = marker_count/dim(union)[1])

## AUC threshold selection -----------------------------------------------------
# Initialize empty df to store rank agg. lists benchmarked against gold set
select.AUC.df <- data.frame(
  Position = integer(),
  ENSG.ID = character(),
  Gene.type = character(),
  Gene.name = character(),
  Chromosome = character(),
  ENSP.ID = character(),
  stat = double(),
  enr = integer(),
  idx_ROC = integer(),
  auc_threshold = double(),
  num_contrasts = integer()
)

# produce rank agg. lists using different AUC thresholds for contrast inclusion
#for (i in seq(0, 1, 0.1)){
for (i in c(0.5, 0.55, 0.6, 0.65, 0.7, 0.75)){ # at 0.8 no contrasts are included and going over .65 there are too few
  print(i)
  # select contrasts
  selected.contrasts <- select_contrast(ranked.df, i)
  print(sum(selected.contrasts$include == "yes" ))
  selected.contrasts <- filter(selected.contrasts, include == "yes")
  selected.contrasts.vec <- selected.contrasts %>% filter(include == "yes") %>% distinct(Contrast)
  # rank aggregation
  rank <- rank_aggregation(selected.contrasts.vec$Contrast, "protein_coding")
  
  # count marker occurences
  rank.markers <- count_markers_rank(rank, union)
  rank.markers$auc_threshold <- i
  rank.markers$num_contrasts <- length(selected.contrasts.vec$Contrast)
  rank.markers$num_diseases <- length(unique(selected.contrasts$case))
  rank.markers$num_tissues <- length(unique(selected.contrasts$tissue))
  select.AUC.df <- rbind(select.AUC.df, rank.markers)
}

select.AUC.df <- select.AUC.df %>%
  group_by(auc_threshold) %>%
  group_modify(~ add_row(.x, 
                         Position = dim(filter(Final_Annotation_List, Gene.type == "protein_coding"))[1],
                         idx_ROC = dim(filter(Final_Annotation_List, Gene.type == "protein_coding"))[1],
                         enr = dim(union)[1],
                         auc_threshold = .x$auc_threshold[1]))


# create proportion of genes covered columns (0 to 1)
select.AUC.df <- select.AUC.df %>%
  group_by(auc_threshold) %>%
  mutate(p_annotation = Position/dim(Final_Annotation_List)[1],
         p_markers = enr/dim(union)[1],
         p_ROC = idx_ROC/dim(Final_Annotation_List)[1])

# Calculate AUC of resulting rank agg. lists
select.AUC.df <- select.AUC.df %>%
  group_by(auc_threshold) %>%
  nest() %>%
  mutate(AUC = map(data, calculate_auc_rank)) %>%
  unnest(cols = c(data, AUC)) %>%
  ungroup() %>% 
  group_by(auc_threshold) %>%
  nest() %>%
  mutate(AUC_p = map(data, calculate_auc_p)) %>%
  unnest(cols = c(data, AUC_p)) %>%
  ungroup()

write_tsv(select.AUC.df, "data/03_AUC_threshold_selection.tsv")

# compare AUCs and diseases and tissues included
summary_table <- distinct(select.AUC.df, AUC, AUC_ROC, AUC_p, AUC_ROC_p, auc_threshold, num_contrasts, num_diseases, num_tissues) %>% drop_na()
summary_table
write_tsv(summary_table, "data/03_AUC_threshold_selection_summary.tsv")

# AUC very similar, look at gene content
ensg.55 <- select.AUC.df$ENSG.ID[select.AUC.df$auc_threshold == 0.55]
ensg.5 <- select.AUC.df$ENSG.ID[select.AUC.df$auc_threshold == 0.5]
ensg.6 <- select.AUC.df$ENSG.ID[select.AUC.df$auc_threshold == 0.6]
ensg.65 <- select.AUC.df$ENSG.ID[select.AUC.df$auc_threshold == 0.65]

length(intersect(ensg.55[1:2000], ensg.6[1:2000]))
length(intersect(ensg.55[1:2000], ensg.5[1:2000]))
length(intersect(ensg.5[1:2000], ensg.6[1:2000]))
length(intersect(ensg.65[1:2000], ensg.6[1:2000]))

# benchmark rank agg. lists
#segment.benchmark <- data.frame(x1 = 0,
#                                y1 = 0,
#                                x2 = dim(filter(Final_Annotation_List, Gene.type == "protein_coding"))[1],
#                                y2 = dim(union)[1])

segment.benchmark <- data.frame(x1 = 0,
                                y1 = 0,
                                x2 = 1,
                                y2 = 1)

#select.AUC.df <- read_tsv("data/03_AUC_threshold_selection.tsv")

select.AUC.df %>%
  ggplot(aes(x = p_ROC, y = p_markers)) +
  geom_path(aes(color=as.factor(auc_threshold))) +
  geom_segment(aes(x = x1,
                   y = y1,
                   xend = x2,
                   yend = y2), 
               data = segment.benchmark,
               linewidth = .5) +
  labs(y = "Proportion of gold standard genes",
       x = "Normalized gene ranks") +
  scale_color_manual(name = "AUC threshold",
                     values = myPalette) +
  plotTheme 

ggsave("figures/03_AUC_threshold_selection.png",
       h = 7,
       w = 10)

select.AUC.df %>%
  ggplot(aes(x = p_ROC, y = p_markers)) +
  geom_path(aes(color=as.factor(auc_threshold))) +
  geom_segment(aes(x = x1,
                   y = y1,
                   xend = x2,
                   yend = y2), 
               data = segment.benchmark,
               linewidth = .5) +
  labs(#title = "Gold Set Enrichment in Tissue Rank Aggregated Lists", 
    #subtitle = paste0("Evaluating AUC threshold for contrast inclusion"), 
    y = "Proportion of gold standard genes",
    x = "Normalized gene ranks") +
  scale_color_manual(name = "AUC threshold",
                     values = myPalette) +
  plotTheme +
  coord_cartesian(clip = "on",
                  xlim = c(0, select.AUC.df$p_ROC[select.AUC.df$Position == 4000][1]),
                  ylim = c(0, select.AUC.df$p_markers[select.AUC.df$Position == 4000][1]))

ggsave("figures/03_AUC_threshold_selection_zoom.png",
       h = 7,
       w = 10)

## Benchmark plots -------------------------------------------------------------
# create dataframe with selected contrasts
# 0.6 yields highest AUC of rank agg list + includes diversity of 
# diseases and tissues
selected.contrasts <- select_contrast(ranked.df, .6) 

selected.contrasts %>%
  filter(include=="yes") %>%
  dplyr::count(case, tissue) %>% 
  arrange(desc(n))

selected.contrasts %>%
  filter(include=="yes") %>%
  summarise(n_datasets = n_distinct(dataset),
            n_contrasts = n_distinct(Contrast),
            n_diseases = n_distinct(case),
            n_tissues = n_distinct(tissue))

# save dataframe of contrast selection
write_tsv(selected.contrasts,
          "data/03_selected_contrasts.tsv")

tissue.contrasts <- filter(selected.contrasts,
                           include == "yes")$Contrast

# create palette with discarded contrasts in gray
tissue.pal <- sapply(selected.contrasts$Contrast,
                     assign_color,
                     df = tissue.contrasts,
                     pal = pal(length(tissue.contrasts)))

scales::show_col(tissue.pal)

# diagonal with 0.5 AUC 
segment.benchmark <- data.frame(x1 = 0,
                                y1 = 0,
                                x2 = dim(filter(Final_Annotation_List, Gene.type == "protein_coding"))[1],
                                y2 = dim(union)[1])

### Plot all contrasts together ------------------------------------------------
tissue.p <- ranked.df %>%
  ggplot(aes(x = idx_ROC, y = marker_count)) +
  geom_path(aes(color = Contrast)) +
  geom_segment(aes(x = x1,
                   y = y1,
                   xend = x2,
                   yend = y2), 
               data = segment.benchmark,
               linewidth = .5) +
  labs(y = "Number of union inflammation markers",
       x = "non-gold std DEGs sorted by stat value") +
  coord_cartesian(clip = "off") +
  guides(size = "none", linetype = "none", colour = "none")+
  scale_color_manual(values = tissue.pal) +
  plotTheme 

tissue.p

ggsave("figures/03_benchmark_allContrasts.png",
       tissue.p,
       height = 7, width = 10)

write_tsv(ranked.df, gzfile("data/03_sorted_contrasts.tsv.gz"))
# ranked.df <- read.delim(gzfile("data/03_sorted_contrasts.tsv.gz"))


### Plot each dataset separately -----------------------------------------------
for (set in selected.contrasts$dataset){
  p <- ranked.df %>%
    filter(str_detect(Contrast, set)) %>%
    ggplot(aes(x = idx_ROC, y = marker_count)) +
    geom_path(aes(color = Contrast)) +
    geom_segment(aes(x = x1,
                     y = y1,
                     xend = x2,
                     yend = y2), 
                 data = segment.benchmark,
                 linewidth = .5) +
    labs(y = "Number of gold standard genes",
         x = "Non-gold standard DEGs") +
    coord_cartesian(clip = "off") +
    scale_color_manual(values = tissue.pal) +
    plotTheme 
  print(p)
  
  ggsave(paste("figures/03_benchmark_", set, ".png", sep=""), 
         p,
         width = 10,
         height = 7)
}

### Plot per tissue ------------------------------------------------------------
selected.contrasts %>% count(tissue)

# max contrasts per tissue is 15; need a discrete palette of 15
selected.contrasts %>% 
  filter(include == "yes") %>% 
  dplyr::count(tissue)

ranked.df <- ranked.df %>%
  left_join(selected.contrasts)

pal = pal(15)

segment.benchmark <- data.frame(x1 = 0,
                                y1 = 0,
                                x2 = 1,
                                y2 = 1)

ranked.df %>%
  filter(tissue == "skin") %>%
  ggplot(aes(x = p_annotation, y = p_markers)) +
  #geom_path(aes(color = Contrast, linetype = include)) +
  geom_path(aes(color = interaction(Contrast, include), 
                linetype = interaction(Contrast, include))) +
  geom_segment(aes(x = x1,
                   y = y1,
                   xend = x2,
                   yend = y2), 
               data = segment.benchmark,
               linewidth = .5) +
  labs(y = "Proportion of gold standard genes",
       x = "Genes ranked by decreasing test statistic") +
  #coord_cartesian(clip = "off") +
  scale_color_manual(values = myPalette) +
  scale_linetype_manual(values = c("solid", "dashed")) +  # Adjust linetypes if needed
  #guides(
  #  color = guide_legend(order = 1),      # Legend for color (Contrast)
  #  linetype = guide_legend(order = 2)    # Legend for linetype (include)
  #) +
  plotTheme 

ranked.df %>%
  filter(tissue == "skin") %>%
  ggplot(aes(x = p_annotation, y = p_markers)) +
  # Map color to Contrast and linetype to include
  geom_path(aes(color = Contrast, linetype = include)) +
  # Add benchmark segments
  geom_segment(aes(x = x1,
                   y = y1,
                   xend = x2,
                   yend = y2), 
               data = segment.benchmark,
               linewidth = .5) +
  # Labels
  labs(y = "Proportion of gold standard genes",
       x = "Genes ranked by decreasing test statistic") +
  # Manual color palette
  scale_color_manual(values = myPalette) +  # Assign colors to the contrasts
  scale_linetype_manual(values = c("yes" = "solid", "no" = "dashed")) +  # Assign linetypes
  # Adjust legend
  guides(
    color = guide_legend(order = 1),      # Color legend for Contrast
    linetype = guide_legend(order = 2)    # Linetype legend for include
  ) +
  # Optional: Remove the title from the legend
  theme(legend.title = element_blank()) +  
  plotTheme

ranked.df %>%
  filter(tissue == "skin") %>%
  ggplot(aes(x = p_annotation, y = p_markers)) +
  # Combine Contrast and include in the mapping
  geom_path(aes(color = interaction(Contrast, include), linetype = include)) +  
  # Add benchmark segments
  geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), 
               data = segment.benchmark, linewidth = .5) +
  # Labels
  labs(y = "Proportion of gold standard genes",
       x = "Genes ranked by decreasing test statistic") +
  # Custom color palette for Contrast
  scale_color_manual(values = myPalette) +  
  # Set linetypes based on 'include' (solid for "yes", dashed for "no")
  scale_linetype_manual(values = c("yes" = "solid", "no" = "dashed")) + 
  # Adjust legend to show both color and linetype together
  guides(color = guide_legend(title = "Contrast", order = 1),
         linetype = guide_legend(title = "Include", order = 2)) +
  # Optional: Remove legend title if desired
  theme(legend.title = element_blank()) +  
  plotTheme



 for (tiss in selected.contrasts$tissue){
  p <- ranked.df %>%
    filter(tissue == tiss) %>%
    ggplot(aes(x = p_ROC, y = p_markers)) +
    geom_path(aes(color = Contrast)) +
    geom_segment(aes(x = x1,
                     y = y1,
                     xend = x2,
                     yend = y2), 
                 data = segment.benchmark,
                 linewidth = 2) +
    labs(y = "Number of union inflammation markers",
         x = "DEGs sorted by stat value") +
    coord_cartesian(clip = "off") +
    scale_color_manual(values = tissue.pal) +
    plotTheme
  
  print(p)
  
  ggsave(paste("figures/03_benchmark_", tiss, ".png", sep=""), 
         p,
         width = 12,
         height = 8)
 }

### Plot example ---------------------------------------------------------------
ranked.df %>%
  filter(Contrast %in% c("02_GSE166925_CD_CTLlarge", "02_GSE175759_LupusNeph_ctl", "02_GSE142530_AH_CTL", "02_GSE231693_IPF_CTL")) %>%
  ggplot(aes(x = p_ROC, y = p_markers)) +
  #geom_path(aes(color = Contrast, linetype = include)) +
  geom_path(aes(color = Contrast, 
                linetype = include)) +
  geom_segment(aes(x = x1,
                   y = y1,
                   xend = x2,
                   yend = y2), 
               data = segment.benchmark,
               linewidth = .5) +
  labs(y = "Proportion of gold standard genes",
       x = "Genes ranked by decreasing test statistic") +
  #coord_cartesian(clip = "off") +
  scale_color_manual(values = myPalette) +
  scale_linetype_manual(values = c(yes="solid", no="dotted")) +  # Adjust linetypes if needed
  #guides(
  #  color = guide_legend(order = 1),      # Legend for color (Contrast)
  #  linetype = guide_legend(order = 2)    # Legend for linetype (include)
  #) +
  plotTheme

ggsave("figures/03_benchmark_example.png",
       h=7, w=10)
