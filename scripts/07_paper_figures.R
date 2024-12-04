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

# Define color palette ---------------------------------------------------------
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


# 1. Benchmarking for contrast inclusion (supplementary) -----------------------
rank.contrasts <- read_tsv("data/03_sorted_contrasts.tsv")
selected.contrasts <- read_tsv("data/03_selected_contrasts.tsv")

segment.proportion <- data.frame(x1 = 0,
                                 y1 = 0,
                                 x2 = 1,
                                 y2 = 1)

selected.contrasts %>% count(tissue)
# max contrasts per tissue is 15; need a discrete palette of 15
selected.contrasts %>% 
  filter(include == "yes") %>% 
  dplyr::count(tissue)

## Make new contrast names -----------------------------------------------------
selected.contrasts %>% count(dataset)
selected.contrasts %>% distinct(case) %>% print(n=30)
selected.contrasts %>% distinct(control) %>% print(n=30)

selected.contrasts <- selected.contrasts %>%
  mutate(new_case = case_when(case == "AD"       ~ "Aderm",
                              case == "Contact"  ~ "Cderm",
                              case == "Alz"      ~ "AD",
                              case == "Cirr"     ~ "AC" , # alcohol-related cirrhosis
                              case == "AC"       ~ "LUAD", # lung adenocarcinoma
                              case == "DiabNeph" ~ "DN",
                              case == "LupusNeph"~ "LN",
                              case == "Acontact" ~ "ACderm",
                              T ~ case),
         new_control = case_when(control == "ctl"           ~ "Ctl",
                                 control == "marginpaired"  ~ "paired",
                                 control == "postpaired"    ~ "paired (post-op)",
                                 control == "CTL"           ~ "Ctl",
                                 control == "Old"           ~ "old age Ctls",
                                 control == "Young"         ~ "young age Ctls",
                                 control == "Ctlsalmon"     ~ "Ctl",
                                 control == "pairedsalmon"  ~ "paired",
                                 control == "MSinactive"    ~ "Inactive MS (unpaired)",
                                 control == "CTLlarge"      ~ "Ctl (large intestine)",
                                 control == "CTLsmall"      ~ "Ctl (small intestine)",
                                 control == "post_treatment"~ "paired",
                                 control == "ctlearly"      ~ "Ctl (early onset)",
                                 control == "ctllate"       ~ "Ctl (late onset)",
                                 control == "ILD_CTL"       ~ "Ctl",
                                 control == "MHL"           ~ "MHL",
                                 control == "MHO"          ~ "MHO",
                                 T ~ control),
         new_tissue = case_when(tissue == "CNS" ~ "Central nervous system",
                                tissue == "skin" ~ "Skin",
                                tissue == "intestine" ~ "Intestine",
                                tissue == "synovium"  ~ "Synovium",
                                tissue == "liver" ~ "Liver",
                                tissue == "lung" ~ "Lung",
                                tissue == "kidney" ~ "Kidney",
                                tissue == "salivary_glands" ~ "Salivary gland",
                                tissue == "fat" ~ "Adipose tissue",
                                tissue == "muscle" ~ "Muscle"),
         new_contrast = if_else(str_detect("paired", new_control),
                                paste(new_case, new_control, dataset, sep = " "),
                                paste(new_case, "vs.", new_control, dataset, sep = " ")))

## Plot ------------------------------------------------------------------------
rank.contrasts <- rank.contrasts %>%
  left_join(selected.contrasts)

selected.contrasts <- selected.contrasts %>% 
  group_by(tissue) %>%
  arrange(dataset)

plot_list <- lapply(unique(rank.contrasts$new_tissue), function(f){
  
  linetype_val <- if_else(filter(selected.contrasts, new_tissue == f)$include == "yes", "solid", "12")
  names(linetype_val) <- filter(selected.contrasts, new_tissue == f)$new_contrast
  
  rank.contrasts %>%
    filter(new_tissue == f) %>%
    ggplot(aes(x = p_ROC, y = p_markers)) +
    geom_path(aes(color = fct_reorder(new_contrast, dataset), linetype = fct_reorder(new_contrast, dataset)),
              linewidth = .5) +  
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), 
                 data = segment.proportion, linewidth = .5) +
    labs(title = f,
         y = "Cumulative proportion of gold standard genes",
         x = "Cumulative proportion of other genes",
         color = "Contrast",
         linetype = "Contrast") +
    scale_color_manual(values = myPalette) +  
    scale_linetype_manual(values = linetype_val) + 
    guides(color = guide_legend(override.aes = list(linewidth = .9))) +
    plotTheme
})
       
plot_grid(plotlist = plot_list,
          nrow = 5,
          labels = c("a", "b", "c", "d", "e",
                     "f", "g", "h", "i", "j"))

ggsave("figures/07_benchmark_per_tissue.png",
       h=35, w=20)

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