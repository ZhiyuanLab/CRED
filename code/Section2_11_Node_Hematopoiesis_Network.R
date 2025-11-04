# '=== Header Start ==='
# Title:       Section2_11_Node_Hematopoiesis_Network
# Author:      Wanqi Li
# Date:        20240321
# Purpose:     analysis and figure for 11-Node Hematopoiesis Network mentioned in Section 2
# '=== Header End ==='

# config --------
setwd("D:/wq/zhiyuanlab/20230717-CRED/CRED_project/code/")
data_folder <- "../data/Simulation_11_Node_Hematopoiesis_Network/"
image_folder <- "../data/Simulation_11_Node_Hematopoiesis_Network/image/"

# import --------
library('viper')
library('ggplot2')
library('Seurat')

# color palette --------
colors_pseudotime <- colorRampPalette(c("#00007F", "blue", "#007FFF", 
                                        "cyan", "#7FFF7F", "yellow", 
                                        "#FF7F00", "red", "#7F0000"))
colors_hemato_celltype = c("#d53e4f", "#fc8d59", "#fee08b", "#99d594", "#3288bd")
names(colors_hemato_celltype) = c("HSC", "Mega.", "Ery.", "Mono.", "Gran.")
levels_hemato_celltype <- names(colors_hemato_celltype)

# useful function --------
Format_Regulon_DoRothEA2Viper <- function (dorothea_regulon) {
  # change regulon from DoRothEA format to Viper format
  dorothea.TF <- unique(dorothea_regulon$tf)
  regulon.Viper <- list()
  for (TF in dorothea.TF) {
    interactionset <- dorothea_regulon[
      dorothea_regulon$tf == TF,]
    geneset <- interactionset$target
    regulon.Viper[[TF]]$tfmode = interactionset$mor
    names(regulon.Viper[[TF]]$tfmode) = geneset
    regulon.Viper[[TF]]$likelihood = rep(1, length(geneset))
  }
  return(regulon.Viper)
}

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

Plot_Dimentional_Reduction <- function(Exp_tsne_Mat, name_title = "t-SNE visualization based on gene expression", 
                                       name_color = "PseudoTime", lab_color = "PseudoTime", max_color = NA,
                                       flag_save_plot = FALSE, path_reduction_plot, type = ".png") {
  # draw dimentional reduction plot
  p.tsne <- ggplot(Exp_tsne_Mat, aes(x = tSNE_1, y = tSNE_2, text = cell)) +
    geom_point(size = 3, aes(colour = .data[[name_color]]))+
    theme_test()+
    labs(
      x = 'tSNE-1',
      y = 'tSNE-2',
      title = eval(name_title),
      color = eval(lab_color)) +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 15),
      plot.title = element_text(size = 18),
      # legend.position = c(0.15, 0.2),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      panel.border = element_rect(linewidth = 2)
    )
  
  if (lab_color == "Expression") {
    p.tsne <- p.tsne +
      scale_color_gradient(low = "gray", high = "firebrick",
                           limits = c(0,3.5))
  } else if (lab_color == "Activity") {
    p.tsne <- p.tsne +
      scale_color_gradient2(low = "navyblue", mid = "gray", high = "firebrick",
                            limits = c(-2.6, 2.6), midpoint = 0)
  } else if (lab_color == "PseudoTime"){
    if (is.na(max_color)) {
      max_color <- max(Exp_tsne_Mat[,name_color])
    }
    p.tsne <- p.tsne +
      scale_color_gradientn(colours = colors_pseudotime(7),
                            limits = c(0, max_color))
  } else if (lab_color == "celltype") {
    p.tsne <- p.tsne +
      scale_color_manual(values = colors_hemato_celltype,
                         # labels = c("HSC", "Gran.", "Mono.", "Ery.", "Mega.")
                         )
  } else {
    cat("WARNING: unknown color lab!")
  }
  
  if (flag_save_plot) {
    ggsave(paste0(path_reduction_plot, type), p.tsne, 
           width = 5, height = 4)
  }
  return(p.tsne)
}

Plot_Count_by_Celltype <- function(meta.info.celltype_count, name_x, name_y, max_y = NA,
                                   name_group, lab_color, name_title,
                                   flag_save_plot = FALSE, path_plot, type = ".png",
                                   width = 5, height = 4) {
  if (is.na(max_y)) {
    max_y = max(meta.info.celltype_count[,name_y])
  }
  p <- ggplot(meta.info.celltype_count, aes(x = .data[[name_group]], y = .data[[name_y]], fill = .data[[name_group]])) +
    geom_col() +
    scale_fill_manual(values = colors_hemato_celltype,
                       # labels = c("HSC", "Gran.", "Mono.", "Ery.", "Mega.")
                      ) +
    geom_text(aes(label = .data[[name_y]]), size = 3, vjust = -0.5) +
    ylim(0, max_y) +
    labs(
      x = name_x,
      y = name_y,
      title = eval(name_title),
      color = eval(lab_color)) +
    # scale_colour_distiller(palette = "Spectral") +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 15),
      # axis.text.x = element_text(angle = 90, vjust = 0.5),
      plot.title = element_text(size = 18),
      plot.subtitle = element_text(size = 15),
      legend.position = "none",
      # legend.position = c(0.15, 0.2),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      panel.border = element_rect(linewidth = 2, fill = NA),
      panel.background = element_blank()
    )
  
  if (flag_save_plot) {
    ggsave(paste(path_plot, type, sep = ""), p, 
           width = width, height = height)
  }
  return(p)
}

Plot_Scatter_XY <- function(Exp_tsne_Mat, name_x, name_y, name_color, lab_color, alpha = 1,
                            max_x = NA, min_x = NA, flag_cor = FALSE, name_cor_method = "pearson",
                            name_title,
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
      scale_color_manual(values = colors_hemato_celltype,
                         # labels = c("HSC", "Gran.", "Mono.", "Ery.", "Mega.")
                         )
  } else {
    cat("WARNING: unknown color lab!")
  }
  
  if (flag_save_plot) {
    ggsave(paste0(path_plot, type), p.tsne, 
           width = 5, height = 4)
  }
  return(p.tsne)
}

Plot_Correlation_Bar_by_Celltype <- function(Exp_tsne_Mat_cluster, name_x, name_y, name_group, lab_color,
                                             name_title, name_subtitle,
                                             name_cor_method = "pearson", min_y = NA, max_y = NA,
                                             flag_save_plot = FALSE, path_plot, type = ".png") {
  clusters <- unique(Exp_tsne_Mat_cluster[,name_group])
  cor_by_cluster <- data.frame("cluster" = c(),
                               "cor" = c())
  for (group in clusters) {
    Exp_tsne_Mat_cluster.tmp <- Exp_tsne_Mat_cluster[Exp_tsne_Mat_cluster[,name_group] == group,]
    corr <- cor(Exp_tsne_Mat_cluster.tmp[,c(name_x, name_y)], method = name_cor_method)[name_x, name_y]
    cor_by_cluster <- rbind(cor_by_cluster, c(group, corr))
  }
  colnames(cor_by_cluster) <- c("cluster", "cor")
  cor_by_cluster$cluster<-factor(cor_by_cluster$cluster,
                                     levels = levels(Exp_tsne_Mat_cluster[, name_group]))
  cor_by_cluster$cor <- as.numeric(cor_by_cluster$cor)
  
  if (is.na(min_y)) {
    min_y = min(cor_by_cluster$cor)
  }
  if (is.na(max_y)) {
    max_y = max(cor_by_cluster$cor)
  }
  
  p.tsne <- ggplot(cor_by_cluster, aes(x = cluster, y = cor, fill = cluster)) +
    geom_col()+
    # scale_fill_brewer(palette = "Spectral") +
    scale_fill_manual(values = colors_hemato_celltype) +
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
      axis.title = element_text(size = 15),
      plot.title = element_text(size = 18),
      plot.subtitle = element_text(size = 15),
      # legend.position = c(0.15, 0.2),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      panel.border = element_rect(linewidth = 2)
    )
  
  if (flag_save_plot) {
    ggsave(paste(path_plot, type, sep = ""), p.tsne, 
           width = 5, height = 4)
  }
  return(p.tsne)
}

Plot_Boxplot_XY <- function(Exp_tsne_Mat, focus_gene1, focus_gene2, name_y, name_color, lab_color, alpha = 1,
                            max_x = NA, min_x = NA, flag_cor = FALSE, name_cor_method = "pearson",
                            name_title, width = 5, height = 4,
                            flag_save_plot = FALSE, path_plot, type = ".png") {
  focus_df <- data.frame(
    "gene" = c(rep(focus_gene1, nrow(Exp_tsne_Mat)), rep(focus_gene2, nrow(Exp_tsne_Mat))),
    name_y = c(Exp_tsne_Mat[, focus_gene1], Exp_tsne_Mat[, focus_gene2])
  )
  focus_df$celltype <- rep(Exp_tsne_Mat$celltype, 2)
  colnames(focus_df) <- c("gene", name_y, "celltype")
  p <- ggplot(focus_df, aes(x = gene, y = .data[[name_y]])) +
    geom_boxplot(size = 1.5, aes(colour = gene)) +
    scale_color_manual(values = c("#cf181d", "#2e4785")) +
    # scale_x_discrete(labels = c(focus_gene1, focus_gene2)) +
    # stat_compare_means(
    #   comparisons = list(c("Ventricle", "Outflow tract / large vessels")),
    #   method = "t.test",
    #   label = "p.signif"
    # ) +
    labs(
      y = name_y,
      title = "",
      color = lab_color)+
    theme(
      axis.line = element_line(linewidth = 1),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 15),
      axis.title.x = element_blank(),
      plot.title = element_text(size = 18),
      # legend.position = c(0.15, 0.2),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      # panel.border = element_rect(size = 2)
      panel.background = element_blank(),
    )
  if (flag_save_plot) {
    ggsave(paste0(path_plot, type), p, 
           width = width, height = height)
  }
  return(p)
}

# ----
# 0 read data --------
ExpMat <- read.csv(paste0(data_folder, "ExpressionData.csv"), 
                   stringsAsFactors = FALSE,
                   row.names = 1, header = T, encoding = "UTF-8")
dim(ExpMat) # [1]   11 8000
# meta info with annotation and prediction
meta_info <- read.csv(paste0(data_folder, "PseudoTime.csv"),
                      stringsAsFactors = FALSE,
                      row.names = 1, header = T, encoding = "UTF-8")
dim(meta_info) # [1] 8000    3
# regulon
regulon <- read.csv(paste0(data_folder, "regulon_known.csv"),
                    stringsAsFactors = FALSE,
                    header = T, encoding = "UTF-8")
regulon_Viper <- Format_Regulon_DoRothEA2Viper(regulon)

# ----
# 1 calculate activity --------
ActMat <- Calculate_Activity(ExpMat, meta_info, regulon_Viper, minsize = 1)
dim(ActMat) # [1]   11 8000
write.csv(ActMat, file = paste0(data_folder, "Activity_known.csv"),
          na = 'NA', fileEncoding="UTF-8")
ActMat <- read.csv(paste0(data_folder, "Activity_known.csv"), stringsAsFactors = FALSE,
                   row.names = 1, header = T, encoding = "UTF-8", check.names = FALSE)
ActMat <- as.data.frame(t(ActMat))

# ----
# 2 do dimentional reduction and clustering --------
serobj <- CreateSeuratObject(counts = ExpMat,
                             meta.data = meta_info,
                             min.cells = 0,
                             min.features = 0 + 1)
serobj <- NormalizeData(serobj,
                        normalization.method = "LogNormalize",
                        scale.factor = 10000 # mean
                        )
serobj <- ScaleData(serobj,
                    features = rownames(ExpMat))
serobj <- RunPCA(serobj,
                 features = rownames(ExpMat),
                 verbose = FALSE, npcs = 10)
ElbowPlot(serobj, ndims = 10)
serobj <- FindNeighbors(serobj, dims = 1:10)
serobj <- FindClusters(serobj, resolution = 0.09)#1.6
serobj <- RunTSNE(serobj, dims = 1:10, check_duplicates = FALSE)
# DimPlot(object = serobj,reduction = "tsne",
#         group.by = "seurat_clusters",
#         # shape.by = "panel_id",
#         pt.size = 3)
# save reduction result
tsne_df <- serobj@reductions[["tsne"]]@cell.embeddings
write.csv(tsne_df, paste0(data_folder, "tsne.csv"), row.names=TRUE)
# save cluster result
meta_info_cluster <- merge(meta_info, 
                           serobj@meta.data[,"seurat_clusters", drop = FALSE],
                           by = 'row.names')
rownames(meta_info_cluster) <- meta_info_cluster[,1]
meta_info_cluster <- meta_info_cluster[,-1]
meta_info_cluster$celltype <- c()
meta_info_cluster[meta_info_cluster$seurat_clusters == 2, "celltype"] = "HSC"
meta_info_cluster[meta_info_cluster$seurat_clusters == 0, "celltype"] = "Gran."
meta_info_cluster[meta_info_cluster$seurat_clusters == 1, "celltype"] = "Mono."
meta_info_cluster[meta_info_cluster$seurat_clusters == 3, "celltype"] = "Ery."
meta_info_cluster[meta_info_cluster$seurat_clusters == 4, "celltype"] = "Mega."
meta_info_cluster$celltype <- factor(meta_info_cluster$celltype,
                                     levels = levels_hemato_celltype)
write.csv(meta_info_cluster, paste0(data_folder, "meta.info.csv"), row.names=TRUE)
meta_info_cluster <- read.csv(paste0(data_folder, "meta.info.csv"), 
                    stringsAsFactors = FALSE,
                    row.names = 1, header = T, encoding = "UTF-8")

# count cell number of each celltype
celltype_count <- aggregate(meta_info_cluster[,"celltype"], by = list(meta_info_cluster[,"celltype"]), FUN = length)
colnames(celltype_count) <- c("celltype", "count")
celltype_count$celltype <- factor(celltype_count$celltype,
                                  levels = levels_hemato_celltype)

# Supp. Fig. 2B
Plot_Count_by_Celltype(celltype_count, name_x = "celltype", name_y = "count", max_y = 3900, 
                       name_group = "celltype", lab_color = "celltype",
                       name_title = "Counts of different celltype",
                       flag_save_plot = TRUE, path_plot = paste0(image_folder, "celltype_count"), type = ".pdf")

# ----
# 3 plot dimentional reduction --------
tsne_df <- read.csv(paste0(data_folder, "tsne.csv"), 
                    stringsAsFactors = FALSE,
                    row.names = 1, header = T, encoding = "UTF-8")
colnames(tsne_df) <- c("tSNE_1", "tSNE_2")
# merge expression, position and annotation
Exp_tsne_info_Mat <- merge(meta_info_cluster, tsne_df,
                           by = 'row.names')
rownames(Exp_tsne_info_Mat) <- Exp_tsne_info_Mat[,1]
colnames(Exp_tsne_info_Mat)[1] <- "cell"
Exp_tsne_info_Mat <- merge(Exp_tsne_info_Mat, t(ExpMat),
                           by = 'row.names')
rownames(Exp_tsne_info_Mat) <- Exp_tsne_info_Mat[,1]
Exp_tsne_info_Mat <- Exp_tsne_info_Mat[, -1]

# plot reduction
# Fig. 2B
Plot_Dimentional_Reduction(Exp_tsne_info_Mat, name_title = "t-SNE visualization\nbased on gene expression", 
                           name_color = "Time", lab_color = "PseudoTime", max_color = 800,
                           flag_save_plot = TRUE, path_reduction_plot = paste0(image_folder, "tsne"), type = ".png")
# Fig. 2C
Plot_Dimentional_Reduction(Exp_tsne_info_Mat, name_title = "t-SNE visualization\nbased on gene expression", 
                           name_color = "celltype", lab_color = "celltype",
                           flag_save_plot = TRUE, path_reduction_plot = paste0(image_folder, "tsne_celltype"), type = ".png")
# plot expression on reduction
# Supp. Fig. 2A
focus_gene <- "Gata1" # change this
Plot_Dimentional_Reduction(Exp_tsne_info_Mat, name_title = paste0("Gene expression of ",focus_gene, "\nwith t-SNE embedding"), 
                           name_color = focus_gene, lab_color = "Expression",
                           flag_save_plot = TRUE, path_reduction_plot = paste0(image_folder, "tsne_exp_", focus_gene), type = ".png")

# merge activity and position
Act_tsne_info_Mat <- merge(meta_info_cluster, tsne_df,
                      by = 'row.names')
rownames(Act_tsne_info_Mat) <- Act_tsne_info_Mat[,1]
colnames(Act_tsne_info_Mat)[1] <- "cell"
Act_tsne_info_Mat <- merge(Act_tsne_info_Mat, ActMat,
                           by = 'row.names')
rownames(Act_tsne_info_Mat) <- Act_tsne_info_Mat[,1]
Act_tsne_info_Mat <- Act_tsne_info_Mat[, -1]

# plot activity on reduction
# Supp. Fig. 2C
focus_gene <- "Gata1" # change this
Plot_Dimentional_Reduction(Act_tsne_info_Mat, name_title = paste0("Gene activity of ",focus_gene, "\nwith t-SNE embedding"), 
                           name_color = focus_gene, lab_color = "Activity",
                           flag_save_plot = TRUE, path_reduction_plot = paste0(image_folder, "tsne_act_", focus_gene), type = ".pdf")

# ----
# 4 plot scatter of X and Y --------
# Supp. Fig. 2D
focus_gene1 <- "Gata1"
focus_gene2 <- "Pu1"
# expression
Plot_Scatter_XY(Exp_tsne_info_Mat, name_x = focus_gene1, name_y = focus_gene2, name_color = "celltype", lab_color = "celltype", 
                max_x = 3.5, min_x = 0,
                name_title = paste0("Expression of gene ", focus_gene1, " and ", focus_gene2),
                flag_save_plot = TRUE, path_plot = paste0(image_folder, "scatter_exp_", focus_gene1, "_", focus_gene2), type = ".pdf")
# activity
Plot_Scatter_XY(Act_tsne_info_Mat, name_x = focus_gene1, name_y = focus_gene2, name_color = "celltype", lab_color = "celltype", 
                max_x = 2.5, min_x = -2.5,
                name_title = paste0("Activity of gene ", focus_gene1, " and ", focus_gene2),
                flag_save_plot = TRUE, path_plot = paste0(image_folder, "scatter_act_", focus_gene1, "_", focus_gene2), type = ".pdf")

# focus on celltype
# Fig. 2DE, Supp. Fig. 2E
focus_celltype <- "HSC"
# cor = cor.test(Act_tsne_info_Mat[Act_tsne_info_Mat$celltype == focus_celltype,"Gata1"], 
#                Act_tsne_info_Mat[Act_tsne_info_Mat$celltype == focus_celltype,"Pu1"], 
#                method = "pearson")
Plot_Scatter_XY(Exp_tsne_info_Mat[Exp_tsne_info_Mat$celltype == focus_celltype, ], 
                name_x = focus_gene1, name_y = focus_gene2, name_color = "celltype", lab_color = "celltype", alpha = 0.6,
                # max_x = 3.5, min_x = 0,
                flag_cor = TRUE,
                name_title = paste0("Expression of gene ", focus_gene1, " and ", focus_gene2),
                flag_save_plot = TRUE, path_plot = paste0(image_folder, "scatter_exp_", focus_gene1, "_", focus_gene2, "_", focus_celltype), type = ".pdf")
Plot_Scatter_XY(Act_tsne_info_Mat[Act_tsne_info_Mat$celltype == focus_celltype, ], 
                name_x = focus_gene1, name_y = focus_gene2, name_color = "celltype", lab_color = "celltype", alpha = 0.6,
                # max_x = 2.5, min_x = -2.5, 
                flag_cor = TRUE,
                name_title = paste0("Activity of gene ", focus_gene1, " and ", focus_gene2),
                flag_save_plot = TRUE, path_plot = paste0(image_folder, "scatter_act_", focus_gene1, "_", focus_gene2, "_", focus_celltype), type = ".pdf")
# Fig. 2DE, Supp. Fig. 2E
name_y <- "Expression"
Plot_Boxplot_XY(Exp_tsne_info_Mat[Exp_tsne_info_Mat$celltype == focus_celltype, ], 
                focus_gene1, focus_gene2, name_y, name_color = "gene", lab_color = "gene", alpha = 1,
                max_x = NA, min_x = NA, flag_cor = FALSE, name_cor_method = "pearson",
                name_title = paste0("Expression of gene ", focus_gene1, " and ", focus_gene2),
                width = 4, height = 4,
                flag_save_plot = TRUE, path_plot = paste0(image_folder, "boxplot_exp_", focus_gene1, "_", focus_gene2, "_", focus_celltype), type = ".pdf")
name_y <- "Activity"
Plot_Boxplot_XY(Act_tsne_info_Mat[Act_tsne_info_Mat$celltype == focus_celltype, ], 
                focus_gene1, focus_gene2, name_y, name_color = "gene", lab_color = "gene", alpha = 1,
                max_x = NA, min_x = NA, flag_cor = FALSE, name_cor_method = "pearson",
                name_title = paste0("Activity of gene ", focus_gene1, " and ", focus_gene2),
                width = 4, height = 4,
                flag_save_plot = TRUE, path_plot = paste0(image_folder, "boxplot_act_", focus_gene1, "_", focus_gene2, "_", focus_celltype), type = ".pdf")

# ----
# 5 plot correlation by cell type --------
focus_gene1 <- "Gata1"
focus_gene2 <- "Pu1"
name_cor_method <- "pearson"
# expression
Plot_Correlation_Bar_by_Celltype(Exp_tsne_info_Mat, name_x = focus_gene1, name_y = focus_gene2, name_group = "celltype", lab_color = "celltype",
                                 name_title = "Expression correlation of cell clusters",
                                 name_subtitle = paste0(focus_gene1, " VS ", focus_gene2, ", ", name_cor_method, " correlation"),
                                 name_cor_method = name_cor_method, min_y = -1, max_y = 1,
                                 flag_save_plot = TRUE, 
                                 path_plot = paste0(image_folder, "correlation_bar_exp_", focus_gene1, "_", focus_gene2, "_", name_cor_method), 
                                 type = ".png")
# activity
Plot_Correlation_Bar_by_Celltype(Act_tsne_info_Mat, name_x = focus_gene1, name_y = focus_gene2, name_group = "celltype", lab_color = "celltype",
                                 name_title = "Activity correlation of cell clusters",
                                 name_subtitle = paste0(focus_gene1, " VS ", focus_gene2, ", ", name_cor_method, " correlation"),
                                 name_cor_method = name_cor_method, min_y = -1, max_y = 1,
                                 flag_save_plot = TRUE, 
                                 path_plot = paste0(image_folder, "correlation_bar_act_", focus_gene1, "_", focus_gene2, "_", name_cor_method), 
                                 type = ".png")
