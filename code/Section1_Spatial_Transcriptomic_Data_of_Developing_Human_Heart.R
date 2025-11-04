# '=== Header Start ==='
# Title:       Section1_Spatial_Transcriptomic_Data_of_Developing_Human_Heart
# Author:      Wanqi Li
# Date:        20240320
# Purpose:     Analysis and figure for Spatial_Transcriptomic mentioned in Section 1
# Source Data: Asp, M. et al. A Spatiotemporal Organ-Wide Gene Expression and 
#              Cell Atlas of the Developing Human Heart. 
#              Cell 179, 1647-1660.e19 (2019).
# '=== Header End ==='

# config --------
setwd("E:/wq/zhiyuanlab/20230717-CRED/CRED_project/code/")
data_folder <- "../data/Dataset_Spatial_Transcriptomic_Data_of_Human_Developing_Heart/"
image_folder <- "../data/Dataset_Spatial_Transcriptomic_Data_of_Human_Developing_Heart/image/"

# import --------
library('viper')
library('stringr')
library('ggplot2')
library('ggpubr')
library('gridExtra')

# color palette --------
colors_heart_region = c('#312E82', '#AF4499', '#1B793D', '#BABABA')
names(colors_heart_region) = c("Ventricle", "Outflow tract / large vessels", "Atrium", "Other")
colors_expression = "Spectral"
colors_activity = "Spectral"

# useful function --------
Calculate_Activity <- function(ExpMat.logcpm, meta.info, regulon, method = "viper", minsize = 4) {
  # generate ExpressionSet ----
  # create phenotype data
  # meta data should include sample name
  meta.info$sample.name <- rownames(meta.info)
  # rownames of meta data
  meta.description <- data.frame(labelDescription = colnames(meta.info))
  # create AnnotatedDataFrame object
  phenoData <- new("AnnotatedDataFrame", data = meta.info, 
                   varMetadata = meta.description)
  # creats ExpressionSet object
  ExpSet <- ExpressionSet(as.matrix(ExpMat.logcpm),
                          phenoData = phenoData)
  
  if (method == "viper") {
    viper.all <- viper(ExpSet, regulon, method = "ttest", minsize = minsize,
                       eset.filter = FALSE, cores = 1, verbose = TRUE)
    activity <- viper.all@assayData[["exprs"]]
  }
  return(activity)
}

Plot_Exp_Act_by_region <- function(ExpMat_logcpm, ActMat, meta_info, 
                                   focus_gene, major_group = c("Ventricle", "Outflow tract / large vessels"),
                                   flag_save_plot = FALSE, path_plot, type = ".png") {
  focus_df <- merge(meta_info, t(ExpMat_logcpm)[, focus_gene, drop = FALSE], by = "row.names")
  rownames(focus_df) <- focus_df$Row.names
  focus_df <- focus_df[,-1]
  p1 <- ggplot(focus_df, aes(x = new_x, y = new_y)) +
    geom_point(aes(color = .data[[focus_gene]]), size = 2.5) +
    scale_colour_distiller(palette = colors_expression) +
    labs(
      x = "",
      y = "",
      title = paste("Expression of ", focus_gene, sep = ""),
      color = "Log2 CPM") +
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 18),
      # legend.position = c(0.15, 0.2),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      panel.border = element_blank(),
      panel.background = element_blank(),
    )
  
  # boxplot only show major groups
  focus_df_major_group <- focus_df[focus_df$region %in% major_group,]
  p2 <- ggplot(focus_df_major_group, aes(x = region, y = .data[[focus_gene]], fill = region)) +
    geom_boxplot() +
    scale_fill_manual(values = colors_heart_region,
                      labels = c("Ventricle", "Outflow tract /\nlarge vessels", "Atrium", "Other")) +
    scale_x_discrete(labels = c("Ventricle", "Outflow tract /\nlarge vessels", "Atrium", "Other")) +
    stat_compare_means(
      comparisons = list(c("Ventricle", "Outflow tract / large vessels")),
      method = "t.test",
      label = "p.signif"
    ) +
    labs(
      y = "Expression",
      title = "",
      color = "region")+
    theme(
      axis.line = element_line(linewidth = 1),
      axis.text = element_text(size = 14),
      axis.text.x = element_blank(),
      axis.title = element_text(size = 15),
      axis.title.x = element_blank(),
      plot.title = element_text(size = 18),
      # legend.position = c(0.15, 0.2),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      # panel.border = element_rect(size = 2)
      panel.background = element_blank(),
    )
  
  focus_df <- merge(meta_info, t(ActMat)[, focus_gene, drop=FALSE], by = "row.names")
  rownames(focus_df) <- focus_df$Row.names
  focus_df <- focus_df[,-1]
  p3 <- ggplot(focus_df, aes(x = new_x, y = new_y)) +
    geom_point(aes(color = .data[[focus_gene]]), size = 2.5) +
    scale_colour_distiller(palette = colors_activity) +
    labs(
      x = "",
      y = "",
      title = paste("Activity of ", focus_gene, sep = ""),
      color = "Activity") +
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 18),
      # legend.position = c(0.15, 0.2),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      panel.border = element_blank(),
      panel.background = element_blank(),
    )
  
  # boxplot only show major groups
  focus_df_major_group <- focus_df[focus_df$region %in% major_group,]
  p4 <- ggplot(focus_df_major_group, aes(x = region, y = .data[[focus_gene]], fill = region)) +
    geom_boxplot() +
    scale_fill_manual(values = colors_heart_region,
                      labels = c("Ventricle", "Outflow tract /\nlarge vessels", "Atrium", "Other")) +
    scale_x_discrete(labels = c("Ventricle", "Outflow tract /\nlarge vessels", "Atrium", "Other")) +
    stat_compare_means(
      comparisons = list(c("Ventricle", "Outflow tract / large vessels")),
      method = "t.test",
      label = "p.signif"
    ) +
    labs(
      y = "Activity",
      title = "",
      color = "region")+
    theme(
      axis.line = element_line(linewidth = 1),
      axis.text = element_text(size = 14),
      axis.text.x = element_blank(),
      axis.title = element_text(size = 15),
      axis.title.x = element_blank(),
      plot.title = element_text(size = 18),
      # legend.position = c(0.15, 0.2),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      # panel.border = element_rect(size = 2)
      panel.background = element_blank(),
    )
  
  p <- grid.arrange(p1, p2, p3, p4, ncol = 4, widths=c(2,2.5,2,2.5))
  
  if (flag_save_plot) {
    ggsave(paste0(path_plot, type), p, 
           width = 16, height = 4)
  }
  return(p)
}

Plot_Scatter_XY <- function(ExpMat, meta_info, name_x, name_y, name_color, lab_color, 
                            name_title,
                            flag_save_plot = FALSE, path_plot, type = ".png") {
  focus_df <- merge(meta_info, t(ExpMat)[, c(name_x, name_y), drop = FALSE], by = "row.names")
  rownames(focus_df) <- focus_df$Row.names
  focus_df <- focus_df[,-1]
  
  p <- ggplot(focus_df, aes(x = .data[[name_x]], y = .data[[name_y]])) +
    geom_point(size = 3, aes(colour = .data[[name_color]]))+
    theme_test()+
    labs(
      x = name_x,
      y = name_y,
      title = eval(name_title),
      color = eval(lab_color)) +
    # scale_colour_distiller(palette = "Spectral") +
    theme(
      axis.line = element_blank(),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 15),
      plot.title = element_text(size = 18),
      # legend.position = c(0.15, 0.2),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      panel.border = element_rect(linewidth = 2)
    )
  
  if (lab_color == "region") {
    p <- p +
      scale_color_manual(values = colors_heart_region,
                         labels = c("Ventricle", "Outflow tract /\nlarge vessels", "Atrium", "Other"))
  } else {
    cat("WARNING: unknown color lab!")
  }
  
  if (flag_save_plot) {
    ggsave(paste0(path_plot, type), p, 
           width = 5, height = 4)
  }
  return(p)
}

# ----
# 0 read data --------
meta_info = read.table(paste0(data_folder, "meta_data.tsv"), 
                       sep='\t', header=T, row.names=1)
dim(meta_info) # [1] 3111   14
ExpMat <- read.csv(paste0(data_folder, "imputated_gene_mtx_16_marker100.tsv"), 
                   stringsAsFactors = FALSE, sep='\t', 
                   row.names = 1, header = T, encoding = "UTF-8")
colnames(ExpMat) <- toupper(colnames(ExpMat))
ExpMat <- as.data.frame(t(ExpMat))
dim(ExpMat) # [1] 15323   243
# regulon
regulon_TD_Viper <- readRDS(file = "../data/regulon/regulon_TD_Viper.rds")

# ----
# 1 normalize data --------
ExpMat_cpm <- apply(ExpMat,2,function(x){
  x/sum(x) * 10^6
})
ExpMat_logcpm <- log2(ExpMat_cpm + 1)

# ----
# 2 plot by position --------
# use data of section 16
focus_section <- '16'
sample_focus <- meta_info[str_subset(rownames(meta_info), paste0('^', focus_section, 'x')), ]
nrow(sample_focus) # [1] 243 spots
sample_focus$res.0.8 = factor(sample_focus$res.0.8, levels=c(0:9))
sample_focus$region = c()
sample_focus[sample_focus$res.0.8 %in% c(0,1,2,3), "region"] = "Ventricle"
sample_focus[sample_focus$res.0.8 %in% c(4), "region"] = "Atrium"
sample_focus[sample_focus$res.0.8 %in% c(5), "region"] = "Outflow tract / large vessels"
sample_focus[sample_focus$res.0.8 %in% c(6,7,8,9), "region"] = "Other"
sample_focus$region<-factor(sample_focus$region,
                            levels =c("Ventricle", "Outflow tract / large vessels", "Atrium", "Other"))
# Fig. 1D
p_region <- ggplot(sample_focus, aes(x = new_x, y = new_y)) +
  geom_point(aes(color = region), size = 2.5) +
  scale_color_manual(values = colors_heart_region,
                     labels = c("Ventricle", "Outflow tract /\nlarge vessels", "Atrium", "Other")) +
  theme_void()
ggsave(paste0(image_folder, 'section_16_region.pdf'), p_region, width=3, height=3)
ggsave(paste0(image_folder, 'section_16_region.png'), p_region, width=3, height=3)

# ----
# 3 calculate activity --------
ActMat <- Calculate_Activity(ExpMat_logcpm, sample_focus, regulon_TD_Viper, method = "viper", minsize = 4)
write.csv(ActMat, file = paste0(data_folder, "TFactivity_section_16_243_spot_1370_TF_viper_TD_ttest.csv"),
          na = 'NA', fileEncoding="UTF-8")
ActMat <- read.csv(paste0(data_folder, "TFactivity_section_16_243_spot_1370_TF_viper_TD_ttest.csv"), stringsAsFactors = FALSE,
                   row.names = 1, header = T, encoding = "UTF-8", check.names = FALSE)

# ----
# 4 plot expression / activity by position --------
# Fig. 1E
focus_gene <- "MESP1"
focus_gene <- "MEOX2"
Plot_Exp_Act_by_region(ExpMat_logcpm, ActMat, sample_focus, 
                       focus_gene, major_group = c("Ventricle", "Outflow tract / large vessels"),
                       flag_save_plot = TRUE, 
                       path_plot = paste0(image_folder, "section_16_", focus_gene, "_byregion_compare"), 
                       type = ".png")

# ----
# 5 plot scatter of focus gene pair --------
# Fig. 1F
focus_gene1 <- "MESP1"
focus_gene2 <- "MEOX2"
# expression
Plot_Scatter_XY(ExpMat_logcpm, sample_focus, focus_gene1, focus_gene2, name_color = "region", lab_color = "region", 
                name_title = "Expression (log2 CPM)",
                flag_save_plot = TRUE, 
                path_plot = paste0(image_folder, "section_16_", focus_gene1, "_", focus_gene2, "_exp"), 
                type = ".pdf")
# activity
Plot_Scatter_XY(ActMat, sample_focus, focus_gene1, focus_gene2, name_color = "region", lab_color = "region", 
                name_title = "Activity",
                flag_save_plot = TRUE, 
                path_plot = paste0(image_folder, "section_16_", focus_gene1, "_", focus_gene2, "_act"), 
                type = ".png")
