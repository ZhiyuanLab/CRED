# '=== Header Start ==='
# Title:       Section2_Atlas_of_Human_Blood_Cells
# Author:      Wanqi Li
# Date:        20240322
# Purpose:     analysis and figure for Human Blood Cells mentioned in Section 2
# Source Data: Xie, X. et al. Single-cell transcriptomic landscape of human 
#              blood cells. National Science Review 8, nwaa180 (2021).
# '=== Header End ==='

# config --------
setwd("D:/wq/zhiyuanlab/20230717-CRED/CRED_project/code/")
data_folder <- "../data/Dataset_Atlas_of_Human_Blood_Cells/"
image_folder <- "../data/Dataset_Atlas_of_Human_Blood_Cells/image/"

# import --------
library('viper')
library('ggplot2')
library('ggpubr')
library('ggExtra')

# color palette --------
colors_hemato_celltype <- colorRampPalette(c("#d53e4f", "#fc8d59", "#fee08b", "#99d594", "#3288bd"))

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

Plot_Count_by_Celltype <- function(meta_info, name_group, lab_color, name_y, max_y = NA,
                                      name_title, name_subtitle,
                                      flag_save_plot = FALSE, path_plot, type = ".png",
                                      width = 5, height = 4) {
  meta_info_count <- aggregate(meta_info[,name_group], by = list(meta_info[,name_group]), FUN = length)
  colnames(meta_info_count) <- c(name_group, name_y)
  meta_info_count[, name_group] <- factor(meta_info_count[, name_group], levels = levels(meta_info[, name_group]))
  if (is.na(max_y)) {
    max_y = max(meta_info_count[, name_y])
  }
  p <- ggplot(meta_info_count, aes(x = .data[[name_group]], y = .data[[name_y]], fill = .data[[name_group]])) +
    geom_col() +
    geom_text(aes(label = .data[[name_y]]), size = 3, vjust = -0.5) +
    ylim(0, max_y) +
    scale_fill_manual(values = colors_hemato_celltype(nrow(meta_info_count))) +
    labs(
      x = "cell clusters",
      y = name_y,
      title = eval(name_title),
      subtitle = name_subtitle,
      color = eval(lab_color)) +
    # scale_colour_distiller(palette = "Spectral") +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 15),
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      plot.title = element_text(size = 18),
      plot.subtitle = element_text(size = 15),
      # legend.position = c(0.15, 0.2),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      panel.border = element_rect(linewidth = 2, fill = NA),
      panel.background = element_blank()
    )
  
  if (flag_save_plot) {
    ggsave(paste0(path_plot, type), p, 
           width = width, height = height)
  }
  return(p)
}

Plot_Correlation_Bar_by_Celltype <- function(ExpMat, meta_info, name_x, name_y, name_group, lab_color,
                                             name_title, name_subtitle,
                                             name_cor_method = "pearson", min_y = NA, max_y = NA,
                                             flag_save_plot = FALSE, path_plot, type = ".png",
                                             width = 5, height = 4) {
  clusters <- unique(meta_info[,name_group])
  cor_by_cluster <- data.frame("cluster" = c(),
                               "cor" = c())
  for (group in clusters) {
    ExpMat_tmp <- ExpMat[c(name_x, name_y), rownames(meta_info[meta_info[,name_group] == group, ])]
    ExpMat_tmp <- as.data.frame(t(ExpMat_tmp))
    corr <- cor(ExpMat_tmp[,c(name_x, name_y)], method = name_cor_method)[name_x, name_y]
    cor_by_cluster <- rbind(cor_by_cluster, c(group, corr))
  }
  colnames(cor_by_cluster) <- c("cluster", "cor")
  cor_by_cluster$cluster<-factor(cor_by_cluster$cluster,
                                 levels = levels(meta_info[, name_group]))
  cor_by_cluster$cor <- as.numeric(cor_by_cluster$cor)
  
  if (is.na(min_y)) {
    min_y = min(cor_by_cluster$cor)
  }
  if (is.na(max_y)) {
    max_y = max(cor_by_cluster$cor)
  }
  
  p <- ggplot(cor_by_cluster, aes(x = cluster, y = cor, fill = cluster)) +
    geom_col()+
    # scale_fill_brewer(palette = "Spectral") +
    scale_fill_manual(values = colors_hemato_celltype(nrow(cor_by_cluster))) +
    ylim(min_y, max_y) +
    theme_test()+
    labs(
      x = "cell clusters",
      y = "correlation",
      title = eval(name_title),
      subtitle = name_subtitle,
      color = eval(lab_color)) +
    # scale_colour_distiller(palette = "Spectral") +
    theme(
      axis.text = element_text(size = 14),
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      axis.title = element_text(size = 15),
      plot.title = element_text(size = 18),
      plot.subtitle = element_text(size = 15),
      # legend.position = c(0.15, 0.2),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      panel.border = element_rect(linewidth = 2)
    )
  
  if (flag_save_plot) {
    ggsave(paste0(path_plot, type), p, 
           width = width, height = height)
  }
  return(p)
}

Plot_Boxplot_XY <- function(Exp_tsne_Mat, name_x, name_y,
                            max_x = NA, min_x = NA, flag_ttest = FALSE, 
                            name_title, 
                            flag_save_plot = FALSE, path_plot, type = ".png") {
  info_df <- data.frame("Expression" = c(Exp_tsne_Mat[, name_x],
                                         Exp_tsne_Mat[, name_y]),
                        "gene" = c(rep(name_x, nrow(Exp_tsne_Mat)),
                                   rep(name_y, nrow(Exp_tsne_Mat))))
  p.tsne <- ggplot(info_df, aes(x = gene, y = Expression)) +
    geom_boxplot(size = 1.5, aes(colour = gene)) +
    scale_color_manual(values = c("#cf181d", "#2e4785")) +
    theme_test() +
    ylim(0, 15) +
    labs(
      y = "Expression (Log2 CPM)",
      title = eval(name_title),
      color = "gene") +
    # scale_colour_distiller(palette = "Spectral") +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 15),
      plot.title = element_text(size = 18),
      plot.subtitle = element_text(size = 15),
      # legend.position = c(0.15, 0.2),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      panel.border = element_rect(linewidth = 2)
    )
  
  if (!is.na(max_x) | !is.na(min_x)) {
    if (is.na(max_x)) {
      max_x <- max(max(Exp_tsne_Mat[,name_x]),
                   max(Exp_tsne_Mat[,name_y]))
    }
    if (is.na(min_x)) {
      min_x <- min(min(Exp_tsne_Mat[,name_x]),
                   min(Exp_tsne_Mat[,name_y]))
    }
    p.tsne <- p.tsne + 
      xlim(min_x, max_x) +
      ylim(min_x, max_x)
  }
  if (flag_ttest) {
    p.tsne <- p.tsne +
      stat_compare_means(
        comparisons = list(c(name_x, name_y)),
        method = "t.test",
        label = "p.signif",
        tip.length = 0.03,
        bracket.size = 1
      ) 
  }
  if (flag_save_plot) {
    ggsave(paste0(path_plot, type), p.tsne, 
           width = 5, height = 4)
  }
  return(p.tsne)
}


Plot_Scatter_XY <- function(Exp_tsne_Mat, name_x, name_y, name_color, lab_color, alpha = 1,
                            max_x = NA, min_x = NA, flag_cor = FALSE, name_cor_method = "pearson",
                            name_title, flag_margin = FALSE, 
                            flag_save_plot = FALSE, path_plot, type = ".png") {
  p.tsne <- ggplot(Exp_tsne_Mat, aes(x = .data[[name_x]], y = .data[[name_y]])) +
    geom_point(size = 3, aes(colour = .data[[name_color]], alpha = alpha)) +
    theme_test()+
    labs(
      x = name_x,
      y = name_y,
      title = eval(name_title),
      color = eval(lab_color)) +
    # scale_colour_distiller(palette = "Spectral") +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 15),
      plot.title = element_text(size = 18),
      plot.subtitle = element_text(size = 15),
      # legend.position = c(0.15, 0.2),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      panel.border = element_rect(linewidth = 2)
    )
  
  if (!is.na(max_x) | !is.na(min_x)) {
    if (is.na(max_x)) {
      max_x <- max(max(Exp_tsne_Mat[,name_x]),
                   max(Exp_tsne_Mat[,name_y]))
    }
    if (is.na(min_x)) {
      min_x <- min(min(Exp_tsne_Mat[,name_x]),
                   min(Exp_tsne_Mat[,name_y]))
    }
    p.tsne <- p.tsne + 
      xlim(min_x, max_x) +
      ylim(min_x, max_x)
  }
  
  if (flag_cor) {
    corr <- cor.test(Exp_tsne_Mat[, name_x], Exp_tsne_Mat[, name_y], method = name_cor_method)
    name_subtitle <- paste0(name_cor_method, " r = ", round(corr$estimate,3), ", p = ", round(corr$p.value,3))
    p.tsne <- p.tsne +
      labs(subtitle = name_subtitle)
  }
  
  if (lab_color == "Expression") {
    p.tsne <- p.tsne +
      scale_color_gradient(low = "gray", high = "firebrick",
                           limits = c(0,3.5))
  } else if (lab_color == "Activity") {
    p.tsne <- p.tsne +
      scale_color_gradient2(low = "navyblue", mid = "gray", high = "firebrick",
                            limits = c(-2.6, 2.6), midpoint = 0)
  } else if (lab_color == "PseudoTime"){
    p.tsne <- p.tsne +
      scale_color_gradientn(colours = colors_pseudotime(7),
                            limits = c(0, max_color))
  } else if (lab_color == "celltype") {
    p.tsne <- p.tsne +
      scale_color_manual(values = colors_hemato_celltype(length(levels(Exp_tsne_Mat[, name_color]))))
  } else {
    cat("WARNING: unknown color lab!")
  }
  
  if (flag_margin) {
    p.tsne <- ggMarginal(
      p.tsne,
      type = 'density',
      margins = 'both',
      size = 5,
      groupColour = TRUE,
      groupFill = TRUE
    )
  }
  
  if (flag_save_plot) {
    ggsave(paste0(path_plot, type), p.tsne, 
           width = 5, height = 4)
  }
  return(p.tsne)
}

# ----
# 0 read data --------
ExpMat <- read.csv(paste0(data_folder, "ABC_umi_matrix_7551_cells.csv"), stringsAsFactors = FALSE,
                   header = T, encoding = "UTF-8")
ExpMat <- as.data.frame(t(ExpMat))
dim(ExpMat) # [1] 19813  7551
meta_info <- read.csv(paste0(data_folder, "meta_info.csv"), stringsAsFactors = FALSE,
                      header = T, row.names = 1, encoding = "UTF-8")
rownames(meta_info) <- meta_info$cell_id
meta_info$celltype <- factor(meta_info$celltype,
                             levels =c("HSC", "MPP", "CMP", "MEP", "LMPP", "MLP", "BNK", "GMP",
                                       "proB", "preB", "immB", "regB", "naiB", "memB", "plasma",
                                       "ery", 
                                       "CLP", "NKP", "toxiNK", "kineNK",
                                       "CD4T", "CD8T",
                                       "hMDP", "cMOP", "preM", "claM", "interM", "nonM",
                                       "proN", "myeN", "metaN", "matureN"))
meta_info$celltype2 <- factor(meta_info$celltype2,
                             levels =c("HSC", "MPP", "CMP", "CLP", "MEP", "GMP", "BNK",
                                       "Erythrocyte", "Monocyte", "Granulocyte", "T", "B", "NK"))
meta_info$lineage <- factor(meta_info$lineage,
                                 levels =c("HSPC", "B", "Erythrocyte", "NK", "T", "Monocyte", "Neutrophil"))
dim(meta_info) # [1] 7551    10
# dropout samples
focus_gene1 <- "GATA1"
focus_gene2 <- "SPI1"
ExpMat_dropout <- ExpMat[, ExpMat[focus_gene1, ] == 0 &
                           ExpMat[focus_gene2, ] == 0]
dim(ExpMat_dropout) # [1] 19813  3196
meta_info_dropout <- meta_info[colnames(ExpMat_dropout),]

# regulon
regulon_TD_known_Viper <- readRDS(file = "../data/regulon/regulon_TD_known_Viper.rds")

# ----
# 1 normalize data --------
ExpMat_cpm <- apply(ExpMat,2,function(x){
  x/sum(x) * 10^6
})
ExpMat_logcpm <- log2(ExpMat_cpm + 1)
ExpMat_dropout_logcpm <- ExpMat_logcpm[, ExpMat_logcpm[focus_gene1, ] == 0 &
                                         ExpMat_logcpm[focus_gene2, ] == 0]

# ----
# 2 calculate activity --------
# all samples
ActMat <- Calculate_Activity(ExpMat_logcpm, meta_info, regulon_TD_known_Viper, method = "viper", minsize = 4)
dim(ActMat) # [1] 1388 7551
write.csv(ActMat, file = paste0(data_folder,  "TFactivity_7551_sample_1388_TF_viper_TD_ttest.csv"),
          na = 'NA', fileEncoding="UTF-8")
ActMat <- read.csv(paste0(data_folder, "TFactivity_7551_sample_1388_TF_viper_TD_ttest.csv"), stringsAsFactors = FALSE,
                   row.names = 1, header = T, encoding = "UTF-8", check.names = FALSE)
# dropout samples
ActMat_dropout <- Calculate_Activity(ExpMat_dropout_logcpm, meta_info_dropout, regulon_TD_known_Viper, method = "viper", minsize = 4)
dim(ActMat_dropout) # [1] 1388 3196
write.csv(ActMat_dropout, file = paste0(data_folder,  "TFactivity_3196_sample_1388_TF_viper_TD_ttest.csv"),
          na = 'NA', fileEncoding="UTF-8")
ActMat_dropout <- read.csv(paste0(data_folder, "TFactivity_3196_sample_1388_TF_viper_TD_ttest.csv"), stringsAsFactors = FALSE,
                   row.names = 1, header = T, encoding = "UTF-8", check.names = FALSE)

# ----
# 3 plot count of celltype --------
# Supp. Fig. 3A
# all samples
cell_number <- nrow(meta_info)
name_group <- "celltype2"
Plot_Count_by_Celltype(meta_info, name_group = name_group, lab_color = name_group, name_y = "count", max_y = 1850,
                       name_title = paste0("Counts of different ", name_group), 
                       name_subtitle = paste0(cell_number, " samples"), 
                       flag_save_plot = TRUE, path_plot = paste0(image_folder, "counts_by_", name_group, "_", cell_number, "samples"), 
                       type = ".pdf",
                       width = 9, height = 5)

# dropout samples
focus_gene1 <- "GATA1"
focus_gene2 <- "SPI1"
cell_number <- nrow(meta_info_dropout)
name_group <- "celltype"
Plot_Count_by_Celltype(meta_info_dropout, name_group = name_group, lab_color = name_group, name_y = "count", max_y = 400,
                       name_title = paste("Counts of different ", name_group, sep = ""), 
                       name_subtitle = paste0(cell_number, " samples, ", focus_gene1, "=0&", focus_gene2, "=0"), 
                       flag_save_plot = TRUE, path_plot = paste0(image_folder, "counts_by_", name_group, "_", cell_number, "samples"), 
                       type = ".pdf",
                       width = 12, height = 4)
# ----
# 4 plot correlation by celltype --------
focus_gene1 <- "GATA1"
focus_gene2 <- "SPI1"
cell_number <- nrow(meta_info)
name_group <- "celltype"
name_cor_method <- "spearman"
# activity
Plot_Correlation_Bar_by_Celltype(ActMat, meta_info, name_x = focus_gene1, name_y = focus_gene2, name_group = name_group, lab_color = name_group,
                                 name_title = paste0("Activity correlation of cell clusters, ", cell_number, " samples"),
                                 name_subtitle = paste0(focus_gene1, " VS ", focus_gene2, ", ", name_cor_method, " correlation"),
                                 name_cor_method = name_cor_method, # min_y = -1, max_y = 1,
                                 flag_save_plot = TRUE, 
                                 path_plot = paste0(image_folder, "act_cor_by_", name_group, "_", 
                                                    focus_gene1, "_", focus_gene2, "_", name_cor_method, "_", 
                                                    cell_number, "samples"), 
                                 type = ".png",
                                 width = 12, height = 4)


cell_number <- nrow(meta_info_dropout)
# activity
Plot_Correlation_Bar_by_Celltype(ActMat, meta_info_dropout, name_x = focus_gene1, name_y = focus_gene2, name_group = name_group, lab_color = name_group,
                                 name_title = paste0("Activity correlation of cell clusters, ", cell_number, " samples"),
                                 name_subtitle = paste0(focus_gene1, " VS ", focus_gene2, ", ", name_cor_method, " correlation"),
                                 name_cor_method = name_cor_method, # min_y = -1, max_y = 1,
                                 flag_save_plot = TRUE, 
                                 path_plot = paste0(image_folder, "act_cor_by_", name_group, "_", 
                                                    focus_gene1, "_", focus_gene2, "_", name_cor_method, "_", 
                                                    cell_number, "samples"), 
                                 type = ".png",
                                 width = 12, height = 4)

# ----
# 5 plot boxplot of focue gene --------
# expression
# Fig. 2F
focus_gene1 <- "GATA1"
focus_gene2 <- "SPI1"
focus_df_Exp <- merge(meta_info,
                      as.data.frame(t(ExpMat_logcpm[c(focus_gene1, focus_gene2), ])),
                      by = "row.names")
rownames(focus_df_Exp) <- focus_df_Exp[,1]
focus_df_Exp <- focus_df_Exp[, -1]
focus_celltype <- "ery"
Plot_Boxplot_XY(focus_df_Exp[focus_df_Exp$celltype == focus_celltype, ], 
                name_x = focus_gene1, name_y = focus_gene2,
                # max_x = 3.5, min_x = 0,
                # max_x = 4000,
                flag_ttest = TRUE,
                name_title = paste0("Expression of focus gene pair"),
                flag_save_plot = TRUE, path_plot = paste0(image_folder, "boxplot_exp_cpm_", focus_gene1, "_", focus_gene2, "_", focus_celltype), type = ".pdf")

# ----
# 6 plot scatter of focue gene --------
# activity
# Fig. 2G
focus_gene1 <- "GATA1"
focus_gene2 <- "SPI1"
focus_df_Act <- merge(meta_info,
                      as.data.frame(t(ActMat[c(focus_gene1, focus_gene2), ])),
                      by = "row.names")
rownames(focus_df_Act) <- focus_df_Act[,1]
focus_df_Act <- focus_df_Act[, -1]
focus_celltype <- "ery"
Plot_Scatter_XY(focus_df_Act[focus_df_Act$celltype == focus_celltype, ], 
                name_x = focus_gene1, name_y = focus_gene2, name_color = "celltype", lab_color = "celltype", alpha = 0.6,
                # max_x = 3.5, min_x = 0,
                flag_cor = TRUE,
                name_title = paste0("Activity of gene ", focus_gene1, " and ", focus_gene2),
                flag_save_plot = TRUE, path_plot = paste0(image_folder, "scatter_act_", focus_gene1, "_", focus_gene2, "_", focus_celltype), type = ".pdf")


