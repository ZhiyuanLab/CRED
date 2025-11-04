# '=== Header Start ==='
# Title:       Section5_Plot_Unique_Interaction
# Author:      Wanqi Li
# Date:        20240526
# Purpose:     Plot figure for unique interactions mentioned in Section 6
# '=== Header End ==='

# config --------
setwd("E:/wq/zhiyuanlab/20230717-CRED/CRED_project/code/")
data_folder <- "../data/Section5&6_Isoform/"
image_folder <- "../data/Section5&6_Isoform/image/"

# import --------
library('ggplot2')
library('ggpubr')
library('ggvenn')

# color palette --------
colors_GO <- colorRampPalette(c("#d53e4f", "#fc8d59", "#fee08b", "#99d594", "#3288bd"))


# useful function --------
Plot_Unique_Ratio_of_Focus_TF <- function(data_path, focus_gene,
                                          flag_save_plot = FALSE, path_plot, type,
                                          width = 4, height = 4) {
  df <- read.csv(data_path, 
                 stringsAsFactors = FALSE,
                 header = T, encoding = "UTF-8")
  p <- ggplot(df, aes(x = reorder(isoform_name, unique_ratio), y = unique_ratio)) +
    geom_bar(stat = 'identity', width = 0.3) +
    ylim(0,1) +
    labs(title = paste("Unique ratio of ", focus_gene, sep = ""))+
    theme(
      axis.text = element_text(size = 14),
      # axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
      axis.title = element_text(size = 15),
      axis.title.x = element_blank(),
      plot.title = element_text(size = 18),
      # legend.position = c(0.15, 0.2),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      panel.border = element_rect(linewidth = 2, fill = NA),
      panel.background = element_blank(),
    )
  if (flag_save_plot) {
    ggsave(paste0(path_plot, type), p, 
           width = width, height = height)
  }
  return(p)
}

# ----
# 1 regulatory difference of isoforms --------
# 1.1 unique ratio of different isoforms -------
# Fig. S5C
focus_gene <- "HES6"
Plot_Unique_Ratio_of_Focus_TF(data_path = paste0(data_folder, "unique_ratio_", focus_gene, ".csv"),
                              focus_gene = focus_gene,
                              flag_save_plot = TRUE,
                              path_plot = paste0(image_folder, "unique_ratio_", focus_gene), type = ".pdf")

# 1.2 target overlap and difference --------
# get target list, focus on HMGA2-1 4 5
regulon_focus <- read.csv(paste0(data_folder, "regulon_HMGA2.csv"), 
                          stringsAsFactors = FALSE, 
                          header = T, encoding = "UTF-8")
target_list_1 <- regulon_focus[regulon_focus$tf == "TFORF2849-HMGA2", "target"]
target_list_4 <- regulon_focus[regulon_focus$tf == "TFORF2852-HMGA2", "target"]
target_list_5 <- regulon_focus[regulon_focus$tf == "TFORF2853-HMGA2", "target"]
target_list <- list(`HMGA2-1` = target_list_1,
                    `HMGA2-4` = target_list_4,
                    `HMGA2-5` = target_list_5)
p <- ggvenn(target_list, show_elements = FALSE,fill_color = c("#d53e4f", "#fee08b", "#3288bd"),
            stroke_color = "grey30",
            label_sep = "\n", stroke_size = 1.5, set_name_size = 5,
            text_size = 3)+
  theme(
    axis.title = element_blank(),
    plot.title = element_text(size = 18),
    plot.subtitle = element_text(size = 15),
    # legend.position = c(0.15, 0.2),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    panel.background = element_blank()
  )
p
ggsave(paste0(image_folder, "HMGA2_target_overlap", ".png"), p, 
       width = 3, height = 3)

# get overlap target
target_overlap <- intersect(target_list_1, intersect(target_list_4, target_list_5))
regulon_overlap <- regulon_focus[regulon_focus$target %in% target_overlap,]
regulon_overlap <- regulon_overlap[regulon_overlap$tf %in% c("TFORF2849-HMGA2",
                                                             "TFORF2852-HMGA2",
                                                             "TFORF2853-HMGA2"),]
regulon_overlap$ESC_mor <- as.factor(regulon_overlap$ESC_mor)
p <- ggplot() + 
  geom_tile(data = regulon_overlap,
            colour="white",
            aes(x = .data[["target"]],
                y = .data[["tf"]],
                fill = .data[["ESC_mor"]])) +
  
  scale_y_discrete(labels = c("TFORF2849-HMGA2" = "HMGA2-1",
                              "TFORF2852-HMGA2" = "HMGA2-4",
                              "TFORF2853-HMGA2" = "HMGA2-5")) +
  labs(title = "Mode of regulation of overlap targets",
       fill = "MOR") +
  theme(
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
    axis.text.y = element_text(hjust = 0.5),
    axis.title = element_blank(),
    plot.title = element_text(size = 18),
    plot.subtitle = element_text(size = 15),
    # legend.position = c(0.15, 0.2),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    panel.border = element_rect(linewidth = 2, fill = NA),
    panel.background = element_blank()
  )
p
ggsave(paste0(image_folder, "HMGA2_target_overlap_mor", ".pdf"), p, 
       width = 9, height = 3)

# ----
# 2 unique ratio distribution of gene set --------
# Fig. 5E
unique_interaction_focus_all <- read.csv(paste0(data_folder, "unique_interaction_focus_GO_all.csv"), 
                                         stringsAsFactors = FALSE,
                                         header = T, encoding = "UTF-8")

p <- ggplot(unique_interaction_focus_all, aes(x = GO, y = unique_interaction_ratio_mean)) +
  geom_boxplot(aes(color = GO), linewidth = 1) +
  scale_x_discrete(position = "top", expand = expansion(add = c(0.8,1.8))) +
  stat_compare_means(method = "anova", aes(label = "p.format"),
                     label.y.npc = 0, label.x.npc = 1, vjust = -1.8, hjust = 0,
                     size = 6) +
  # scale_colour_discrete(palette = "Spectral") +
  scale_color_manual(values = colors_GO(length(unique(unique_interaction_focus_all$GO)))) +
  coord_flip() +
  theme(
    axis.text = element_text(size = 14),
    # axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5),
    # axis.text.y = element_text(hjust = 0.5),
    axis.title = element_text(size = 15),
    plot.title = element_blank(),
    plot.subtitle = element_text(size = 15),
    legend.position = "none",
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    panel.border = element_rect(linewidth = 2, fill = NA),
    panel.background = element_blank()
  )
ggsave(paste0(image_folder, "unique_interaction_ratio_of_GO", ".png"), p, 
       width = 8, height = 4)





