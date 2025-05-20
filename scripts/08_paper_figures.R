# Clear workspace --------------------------------------------------------------
rm(list = ls())

# Load libraries ---------------------------------------------------------------
library(tidyverse)
library(cowplot)
library(readr)
library(ComplexUpset)

# Define functions -------------------------------------------------------------
source(file = "scripts/99_project_functions.R")

# Define plot theme ------------------------------------------------------------
plotTheme <- theme_minimal() + theme(
  plot.title = element_text(size = 22, face = "bold", hjust = 0.5), # Smaller title size, centered
  strip.text = element_text(size = 24), # Size for facet labels
  axis.text = element_text(size = 22), # Axis tick labels
  axis.title = element_text(size = 24), 
  legend.position = "right", 
  legend.text = element_text(size = 22), # Legend text size
  legend.title = element_text(size = 24), # Bold legend title
  #panel.grid.major = element_line(color = "grey80", size = 0.5), # Subtle grid lines
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), # No minor grid lines
  plot.margin = unit(c(1.25, 1.25, 1.25, 1.25), "cm"), # Adjust margins for a clean layout
  axis.ticks = element_line(linewidth = 0.5), # Thinner axis ticks
  axis.line = element_line(linewidth = 0.5) # Subtle axis lines
)

plotTheme2 <- theme_bw() + theme(
  plot.title = element_text(size = 22, face = "bold", hjust = 0.5), # Smaller title size, centered
  strip.text = element_text(size = 24), # Size for facet labels
  axis.text = element_text(size = 22), # Axis tick labels
  axis.title = element_text(size = 24), 
  legend.position = "right", 
  legend.text = element_text(size = 16), # Legend text size
  legend.title = element_text(size = 24), # Bold legend title
  #panel.grid.major = element_line(color = "grey80", size = 0.5), # Subtle grid lines
  #panel.grid.major = element_blank(),
  #panel.grid.minor = element_blank(), # No minor grid lines
  plot.margin = unit(c(1.25, 1.25, 1.25, 1.25), "cm") # Adjust margins for a clean layout
  #axis.ticks = element_line(linewidth = 0.5), # Thinner axis ticks
  #axis.line = element_line(linewidth = 0.5) # Subtle axis lines
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
#rank.contrasts <- read_tsv(gzfile("data/03_sorted_contrasts.tsv.gz"))
# reading with vroom fails
rank.contrasts <- read.delim(gzfile("data/03_sorted_contrasts.tsv.gz"))
selected.contrasts <- read_tsv("data/03_selected_contrasts.tsv")

segment.proportion <- data.frame(x1 = 0,
                                 y1 = 0,
                                 x2 = 1,
                                 y2 = 1)

selected.contrasts %>% dplyr::count(tissue)
# max contrasts per tissue is 15; need a discrete palette of 15
selected.contrasts %>% 
  filter(include == "yes") %>% 
  dplyr::count(tissue)

selected.contrasts %>% 
  filter(include == "yes") %>% 
  dplyr::count(case, dataset, tissue) %>% print(n=33)

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
    coord_fixed(ratio = 1) +
    scale_color_manual(values = myPalette) +  
    scale_linetype_manual(values = linetype_val) + 
    guides(color = guide_legend(override.aes = list(linewidth = .9))) +
    plotTheme +
    theme(legend.justification = c(0,.5))
})
       
plot_grid(plotlist = plot_list,
          nrow = 5,
          labels = c("a", "b", "c", "d", "e",
                     "f", "g", "h", "i", "j"),
          label_size = 24,
          align = "v")

ggsave("figures/08_benchmark_per_tissue.png",
       h=40, w=30)


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
  labs(y = "Cumulative proportion of\ngold standard genes",
       x = "Cumulative proportion of other genes") +
  coord_fixed(ratio = 1) +
  scale_color_manual(values = c("#d7191c", "gray"), 
                     name = "Gene list") +
  plotTheme +
  theme(legend.position="inside", 
        legend.position.inside = c(.8, .2)) +
  guides(color = guide_legend(override.aes = list(linewidth = 3, alpha = 1))) +
  theme(
    axis.title.x = element_text(margin = margin(t = 20)), 
    axis.title.y = element_text(margin = margin(r = 20))  
  )

rank.p

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
  labs(y = "Cumulative count of\ngold standard genes",
       x = "Gene rank") +
  coord_cartesian(clip = "on",
                  xlim = c(0, 3000),
                  ylim = c(0, rank[rank$Position == 3000,]$enr)) +
  geom_vline(xintercept = 2000, linewidth = .5, linetype = "dashed") +
  geom_vline(xintercept = 100, linewidth = .5, linetype = "dashed") +
  plotTheme +
  theme(
    axis.title.x = element_text(margin = margin(t = 20)), 
    axis.title.y = element_text(margin = margin(r = 20))  
  )

## Merge in one figure ---------------------------------------------------------
rank.fig <- plot_grid(rank.p, rank.zoom,
                      labels = c('a', 'b'),
                      rel_widths = c(1, 1.25), 
                      rel_heights = c(1, 1),
                      label_size = 24)

rank.fig

ggsave("figures/08_rank_agg_benchmark.png",
       rank.fig,
       h=7,
       w=18)

# 3. GSEA results --------------------------------------------------------------
GSEA.res <- read_tsv("data/06_GSEA_results.tsv")

GSEA.res %>% distinct(Description)

GSEA.res <- GSEA.res %>%
  mutate(Description = case_when(
    Description == "REACTOME_INTERLEUKIN_1_SIGNALING" ~ "REACTOME: IL-1 signaling",                                   
    Description == "GOBP_ACUTE_INFLAMMATORY_RESPONSE" ~ "GOBP: Acute inflammatory response",                                   
    Description == "REACTOME_INTERFERON_SIGNALING" ~ "REACTOME: IFN signaling",                                   
    Description == "WP_INFLAMMATORY_RESPONSE_PATHWAY" ~ "WP: Inflammatory response pathway",                                 
    Description == "Inflammation signature" ~ "Inflammation signature",                                            
    Description == "MSigDB hallmark inflammatory response" ~ "MSigDB: Hallmark inflammatory response",                             
    Description == "GOBP_INFLAMMATORY_RESPONSE" ~ "GOBP: Inflammatory response",                                         
    Description == "REACTOME_INFLAMMASOMES" ~ "REACTOME: Inflammasomes",                                             
    Description == "WP_CELLS_AND_MOLECULES_INVOLVED_IN_LOCAL_ACUTE_INFLAMMATORY_RESPONSE" ~ "WP: Cells and molecules involved in\nlocal acute inflammatory response",
    Description == "GOBP_CYTOKINE_PRODUCTION_INVOLVED_IN_INFLAMMATORY_RESPONSE" ~ "GOBP: Cytokine production in\ninflammatory response",        
    Description == "REACTOME_TNF_SIGNALING" ~ "REACTOME: TNF signaling",                                            
    Description == "GOBP_CHRONIC_INFLAMMATORY_RESPONSE" ~ "GOBP: Chronic inflammatory response",                                 
    Description == "WP_CYTOKINES_AND_INFLAMMATORY_RESPONSE" ~ "WP: Cytokines and inflammatory response"))


### UpSet plot -----------------------------------------------------------------
all.sets <- read_tsv("data/06_inflammationGeneSets.tsv")

all.sets <- all.sets %>% 
  left_join(distinct(dplyr::select(GSEA.res, ID, Description)), 
            by = join_by(gs == ID)) %>%
  dplyr::select(gene, Description)

# Convert TERM2GENE to a binary presence/absence matrix
binary_matrix <- all.sets %>%
  pivot_wider(names_from = Description, values_from = Description, values_fn = length, values_fill = 0) %>%
  mutate(across(where(is.numeric), ~ ifelse(. > 0, 1, 0))) %>%
  column_to_rownames("gene")

ComplexUpset::upset(
  as.data.frame(binary_matrix),
  intersect = colnames(binary_matrix), # Use gene sets as intersections
  width_ratio = 0.2,
  min_size = 2)

ggsave("figures/08_UpSet_geneSets.png", w = 15, h=11)

### Dotplot --------------------------------------------------------------------

gsea_plot <- GSEA.res %>%
  ggplot(aes(x = reorder(dataset_lab,-log10(p.adjust), decreasing = T), y = reorder(Description,-log10(p.adjust)))) + 
  scale_size_continuous(range = c(1, 15)) +
  geom_point(aes(fill = NES, size = -log10(p.adjust)), shape = 21, colour = "black",alpha = 0.7) +
  geom_point(aes(x = dataset_lab,
                 y = reorder(Description,NES)),
             filter(GSEA.res, p.adjust<0.05, !(dataset %in% c("AD.keller", "SS", "OA", "RA"))),
             size = .7) +
  scale_fill_gradient2(midpoint = 0, low = "blue",mid = "white", high = "red4",
                       space = "Lab"
                       ) +
  scale_x_discrete(guide = guide_axis(angle = 40)) +
  ylab(NULL) + xlab(NULL) +
  labs(fill = "NES",
       size = expression(-log[10]("adj. p-value"))) +
  plotTheme2 +
  theme(
    legend.title = element_text(margin = margin(b = 20)) # Add space below title
  ) 

gsea_plot

ggsave("figures/08_GSEA_results.png",
       gsea_plot,
       w = 13, h = 9)

# 4. Volcano plot --------------------------------------------------------------
UC.hansen <- read_tsv("data/06_UC_hansen_processed.tsv")

sig_up <- dim(UC.hansen %>% filter(adj_pvalue < 0.01, logFC > 0))[1]
sig_up_inf <- dim(UC.hansen %>% filter(adj_pvalue < 0.01, logFC > 0, top_2000 == "yes"))[1]
sig_up_inf/sig_up

all <- dim(UC.hansen)[1]
all_inf <- dim(UC.hansen %>% filter(top_2000 == "yes"))[1]
all_inf/all

volcano <- ggplot(UC.hansen, aes(x = logFC, y = -log10(adj_pvalue), color = significance, alpha = significance)) +
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
    title = "Ulcerative colitis vs. healthy controls\n(Schniers et al.)" #(Schniers et al.)
  ) +
  geom_hline(yintercept = 2, linewidth = .3, linetype = "dashed") +
  coord_cartesian(xlim=c(-2, 2)) +
  plotTheme +
  theme(legend.position = "none") +
  coord_fixed(ratio = .6) 

volcano

ggsave("figures/08_volcano_UCprot_usecase.png",
       volcano,
       h = 8, w = 10)

# combined with GSEA
p <- plot_grid(gsea_plot, volcano,
          nrow = 1, 
          rel_widths = c(2, 1), 
          labels = c('a', 'b'),
          label_size = 24)
p

ggsave("figures/08_GSEA_volcano.png",
       p,
       h=12,
       w=20)

#p2 <- plot_grid(gsea_plot, volcano,
#               nrow = 2, 
#               rel_widths = c(1, 1), 
#               rel_heights = c(1.25, 1),
#               labels = c('a', 'b'),
#               label_size = 24)
#p2
#
#ggsave("figures/07_GSEA_volcano_vertical.png",
#       p2,
#       h=20,
#       w=13)

# 5. Correlation to scores -----------------------------------------------------
scores <- read.table("data/05_score_UC_andersen.tsv", sep = "\t")

correlation <- cor(scores$expr.UC.100.m, scores$score)

ggplot(scores, aes(y = expr.UC.100.m, x = score)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red", se = FALSE) + 
  #labs(x = paste(type,"inflammatome"), x = "Colon inflammation grade score", title = "") +
  labs(y = "Inflammation signature-based score", x = "Colon inflammation grade score", title = "") +
  theme_classic() +
  # Add correlation coefficient as text
  annotate("text", y = max(scores$expr.UC.100.m) + 0.3 , x = max(scores$score), 
           label = paste("Pearson r =", round(correlation, 2)), 
           hjust = 1.5, vjust = 3, color = "red", size = 3.5)

ggsave("figures/08_scatter_score_UseCase.png", device = "png", width = 5, height = 3)



