# '=== Header Start ==='
# Title:       Section6_Stemcell-derived Peri-gastruloids
# Author:      Wanqi Li
# Date:        20240513
# Purpose:     Analysis and figure for Stemcell-derived Peri-gastruloids mentioned in Section 6
# Source Data: Liu, L. et al. Modeling post-implantation stages of human development 
#              into early organogenesis with stem-cell-derived peri-gastruloids. 
#              Cell 186, 3776-3792.e16 (2023).
# '=== Header End ==='

# config --------
setwd("E:/wq/zhiyuanlab/20230717-CRED/CRED_project/code/")
data_folder <- "../data/Dataset_Stemcell-derived_Peri-gastruloids/"
image_folder <- "../data/Dataset_Stemcell-derived_Peri-gastruloids/image/"

# import --------
library('ggplot2')
library('ggsci')
library('viper')
library('reshape2')
library('tidyr')
library('pheatmap')
library('ggtree')
library('aplot')
library('ggExtra')
library('ggpubr')

# color palette --------

# useful function --------
Plot_Dimentional_Reduction <- function(Exp_tsne_Mat, name_title = "t-SNE visualization based on gene expression", 
                                       name_color = "PseudoTime", lab_color = "PseudoTime", max_color = NA,
                                       color_scale = FALSE,
                                       flag_save_plot = FALSE, path_reduction_plot, type = ".png") {
  if (color_scale) {
    Exp_tsne_Mat[,name_color] <- scale(Exp_tsne_Mat[,name_color])
  }
  # draw dimentional reduction plot
  p.tsne <- ggplot(Exp_tsne_Mat, aes(x = UMAP1, y = UMAP2)) +
    geom_point(size = 3, aes(colour = .data[[name_color]]))+
    theme_test()+
    labs(
      x = 'UMAP1',
      y = 'UMAP2',
      title = eval(name_title),
      color = eval(lab_color)) +
    # scale_colour_distiller(palette = "Spectral") +
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
      scale_color_gradient2(low = "navyblue", mid = "gray", high = "firebrick")
  } else if (lab_color == "Activity") {
    p.tsne <- p.tsne +
      scale_color_gradient2(low = "navyblue", mid = "gray", high = "firebrick")
  } else if (lab_color == "Cell Type"){
    p.tsne <- p.tsne +
      scale_color_frontiers()
  } else {
    cat("WARNING: unknown color lab!")
  }
  if (flag_save_plot) {
    ggsave(paste0(path_reduction_plot, type), p.tsne, 
           width = 5, height = 4)
  }
  return(p.tsne)
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

Calculate_Mean_Activity_by_Celltype <- function(act.mat, meta.info, name_group, celltype_list) {
  activity_mean_by_celltype <- data.frame(matrix(ncol = length(celltype_list), nrow = nrow(act.mat)))
  colnames(activity_mean_by_celltype) <- celltype_list
  rownames(activity_mean_by_celltype) <- rownames(act.mat)
  for (celltype in celltype_list) {
    act.mat.tmp <- act.mat[, meta.info[colnames(act.mat), name_group] == celltype]
    act.mean <- rowMeans(act.mat.tmp)
    activity_mean_by_celltype[,celltype] <- act.mean[rownames(activity_mean_by_celltype)]
  }
  return(activity_mean_by_celltype)
}

Save_Pheatmap_Pdf <- function(x, filename, width=7, height=7){
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

Calculate_recall_activity <- function(ActMat_Mean_Lineage, TFome_65, regulon, class) {
  recall_TFome <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(recall_TFome) <- c("Rank of TFA", "Percent Recall", "Used Regulome", "Class")
  TFome_intersect <- intersect(TFome_65$V1,rownames(ActMat_Mean_Lineage))
  
  if (class == "calculated") {
    ActMat_Mean_Lineage_Ect <- ActMat_Mean_Lineage[order(ActMat_Mean_Lineage$Ectoderm,
                                                         decreasing = T),]
    ActMat_Mean_Lineage_Endo <- ActMat_Mean_Lineage[order(ActMat_Mean_Lineage$Endoderm,
                                                          decreasing = T),]
    ActMat_Mean_Lineage_Meso <- ActMat_Mean_Lineage[order(ActMat_Mean_Lineage$Mesoderm,
                                                          decreasing = T),]
    # save result
    for (top in seq(0, 0.2, 0.05)) {
      n <- floor(nrow(ActMat_Mean_Lineage) * top)
      temp <- union(rownames(ActMat_Mean_Lineage_Ect)[1:n],rownames(ActMat_Mean_Lineage_Meso)[1:n])
      temp <- union(temp,rownames(ActMat_Mean_Lineage_Endo)[1:n])
      overlap <- intersect(temp, TFome_intersect)
      overlap_ratio <- length(overlap) / length(TFome_intersect)
      recall_TFome <- rbind(recall_TFome, c(top * 100, overlap_ratio, regulon, class))
    }
  } else {
    set.seed(305)
    for (top in seq(0, 0.2, 0.05)) {
      n <- floor(nrow(ActMat_Mean_Lineage) * top)
      ori <- 0
      round <- 1000
      for(i in 1:round){
        sample <- sample(rownames(ActMat_Mean_Lineage), n, replace = FALSE)
        overlap <- intersect(sample, TFome_intersect)
        temp <- length(overlap)/55
        ori <- ori + temp
      }
      overlap_ratio <- ori/round
      recall_TFome <- rbind(recall_TFome, c(top * 100, overlap_ratio, regulon, class))
    }
  }
  colnames(recall_TFome) <- c("Rank of TFA", "Percent Recall", "Used Regulome", "Class")
  return(recall_TFome)
}

Plot_Correlation_of_Exp_and_Act <-  function(Exp.Mat, act.mat, focus_isoform_list, focus_isoform_list_name, focus_gene, cor_method,
                                             name_title = "", name_subtitle = "",
                                             flag_save_plot = FALSE, path_plot, type,
                                             width = 5, height = 4) {
  cor_df <- data.frame(
    "TF" = focus_isoform_list_name,
    "cor" = c(rep(NA, length(focus_isoform_list)))
  )
  rownames(cor_df) <- focus_isoform_list
  for (isoform in focus_isoform_list) {
    corr = cor(t(Exp.Mat[focus_gene,]), t(act.mat[isoform, colnames(Exp.Mat)]), method = cor_method)
    cor_df[isoform, "cor"] <- corr
  }
  
  p <- ggplot(cor_df, aes(x = TF, y = cor, fill = TF)) +
    geom_col()+
    scale_fill_brewer(palette = "Spectral") +
    # ylim(min_y, max_y) +
    # theme_test()+
    # geom_text(aes(label=.data[[name_y]]),size=3,vjust=0)+
    labs(
      x = "isoform",
      y = "correlation",
      title = eval(name_title),
      subtitle = name_subtitle,
      # color = eval(lab_color)
    ) +
    # scale_colour_distiller(palette = "Spectral") +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 15),
      axis.text.x = element_blank(),
      plot.title = element_text(size = 18),
      plot.subtitle = element_text(size = 15),
      # legend = element_blank(),
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

Plot_Heatmap_of_TF_Pair <- function(df_long, name_value, name_x, name_y,
                                    flag_text = FALSE, name_text,
                                    cluster_rows = FALSE, cluster_cols = FALSE,
                                    name_title = "", lab_color,
                                    flag_save_plot = FALSE, path_plot, type,
                                    width = 5, height = 4) {
  df_long <- df_long[, c(name_x, name_y, name_value)]
  df_wide <- spread(df_long, key = name_y, value = name_value)
  rownames(df_wide) <- df_wide[,1]
  df_wide <- df_wide[,-1]
  if (cluster_rows) {
    byrow <- hclust(dist(df_wide))
    df_wide <- df_wide[byrow$order,] 
    v <- ggtree(byrow,layout = "rectangular",branch.length = "none", size = 0.5)+layout_dendrogram() + # 绘制行聚类树
      geom_aline(linetype="solid", linewidth = 0.5) + # 补齐末端
      labs(title = name_title) +
      theme(
        plot.title = element_text(size = 18, hjust = 0.25)
      )
  }
  if (cluster_cols) {
    bycol <- hclust(dist(t(df_wide)))
    df_wide <- df_wide[,bycol$order]
    h <- ggtree(bycol,layout = "rectangular",branch.length = "none", size = 0.5) # 绘制列聚类树
  }
  # change to long format
  df_wide[, name_x] <- row.names(df_wide)
  df_long <- melt(df_wide,id.vars = name_x, variable.name = name_y, value.name = name_value)   #宽数据变为长数据
  df_long[,1] <- factor(df_long[,1], level = rev(unique(df_long[,1])))
  
  # df_long$num <- rep(c(1:26),6)   #绘图时的纵坐标
  # df_long$x <- rep(c(1:6),each = 26)   #绘图时的横坐标
  
  p <- ggplot() + 
    geom_tile(data = df_long,
              colour="white",
              aes(x = .data[[name_x]],
                  y = .data[[name_y]],
                  fill = .data[[name_value]])) +
    scale_fill_gradient2(low = "navyblue", mid = "white", high = "firebrick", midpoint = 0,
                         # limit = c(-1,1),
                         name = lab_color) + 
    scale_y_discrete(position = c("right")) + 
    labs(
      # title = eval(name_title),
      color = eval(lab_color)
    ) +
    scale_color_gradient2(low = "navyblue", mid = "gray", high = "firebrick") +
    theme(
      axis.text = element_text(size = 14),
      axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5),
      axis.text.y = element_text(hjust = 0.5),
      axis.title = element_blank(),
      plot.title = element_blank(),
      plot.subtitle = element_text(size = 15),
      # legend.position = c(0.15, 0.2),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      panel.border = element_rect(linewidth = 2, fill = NA),
      panel.background = element_blank()
    )
  
  if (flag_text) {
    p <- p +
      # add TF name and activity
      geom_text(data = df_long,
                size = 1.8,
                aes(x = .data[[name_x]],
                    y = .data[[name_y]],
                    label = round(.data[[name_text]],3)),
                vjust = 0) 
  }
  
  if (cluster_cols) {
    p <- p %>% insert_left(h,width = 0.3)
  }
  
  if (cluster_rows) {
    p <- p %>% insert_top(v,height = 0.2)
  }
  
  if (flag_save_plot) {
    ggsave(paste(path_plot, type, sep = ""), p, 
           width = width, height = height)
  }
  return(p)
}

Calculate_Correlation_of_Activity_by_Focus <- function(act.mat, focus_TF_pair_list, focus_cell_list,
                                                       cor_method = "pearson", type = "r") {
  cor_df <- focus_TF_pair_list
  for (i in 1:nrow(focus_TF_pair_list)) {
    for (cell_list_name in names(focus_cell_list)) {
      cell_list <- focus_cell_list[[cell_list_name]]
      act.mat.focus <- act.mat[unlist(focus_TF_pair_list[i,1:2]), cell_list]
      act.mat.focus <- t(act.mat.focus)
      corr <- cor.test(act.mat.focus[, 1], act.mat.focus[, 2], method = cor_method)
      # corr <- cor(t(act.mat.focus), method = cor_method) [1,2]
      if (type == "r") {
        cor_df[i, cell_list_name] <- corr$estimat
      } else if (type == "p") {
        cor_df[i, cell_list_name] <- corr$p.value
      } else {
        cat("ERROR: unknown statistic!")
        return()
      }
    }
  }
  return(cor_df)
}

Plot_Activity_Scatter <- function(act.mat, meta.info, focus_gene1, focus_gene2, name_group,
                                  name_x, name_y,
                                  name_title, flag_scale = FALSE, cor_method = "pearson",
                                  flag_save_fig = FALSE, path_fig) {
  info.focus <- act.mat[c(focus_gene1, focus_gene2),]
  info.focus <- as.data.frame(t(info.focus))
  info.focus <- merge(meta.info, info.focus, by = "row.names")
  rownames(info.focus) <- info.focus[,1]
  info.focus <- info.focus[,-1]
  
  if (flag_scale) {
    info.focus[, focus_gene1] <- scale(info.focus[, focus_gene1])
    info.focus[, focus_gene2] <- scale(info.focus[, focus_gene2])
  }
  corr <- cor.test(info.focus[,focus_gene1], info.focus[,focus_gene2], method = cor_method)
  cat(corr$estimat, corr$p.value)
  
  p <- ggplot(info.focus, aes(x = .data[[focus_gene1]], y = .data[[focus_gene2]])) +
    geom_point(size = 2, aes(color = .data[[name_group]]))+
    # xlim(min_x, max_x) +
    # ylim(min_x, max_x) +
    # theme_test()+
    labs(
        x = name_x,
        y = name_y,
      title = eval(name_title),
      #   color = eval(lab_color)
    ) +
    # scale_colour_distiller(palette = "Spectral") +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 15),
      plot.title = element_text(size = 18),
      plot.subtitle = element_text(size = 15),
      # legend.position = c(0.15, 0.2),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 14),
      panel.border = element_rect(linewidth = 2, fill = NA),
      panel.background = element_blank()
    )
  p <- ggMarginal(
    p,
    type = 'density',
    margins = 'both',
    size = 5,
    groupColour = TRUE,
    groupFill = TRUE
  )
  if (flag_save_fig) {
    ggsave(paste(path_fig, sep = ""), p, 
           width = 5, height = 4)
  }
  return(p)
}

# ----
# 0 read data --------
ExpMat <- read.csv(paste0(data_folder, "ExpMat_lognorm_1930_cells_33538_genes.csv"), stringsAsFactors = FALSE,
                    row.names = 1, header = T, encoding = "UTF-8")
ExpMat <- as.data.frame(t(ExpMat))
dim(ExpMat) # [1] 33538  1930
meta_info <- read.csv(paste0(data_folder, "meta_info_1930_cells.csv"), stringsAsFactors = FALSE,
                      row.names = 1, header = T, encoding = "UTF-8")
dim(meta_info) # [1] 1930    3
embedding_umap <- read.csv(paste0(data_folder, "embedding_umap_1930_cells.csv"), stringsAsFactors = FALSE,
                           row.names = 1, header = T, encoding = "UTF-8")
dim(embedding_umap) # [1] 1930    2

# ----
# 1 plot UMAP --------
Exp_tsne_Mat <- merge(meta_info, embedding_umap, by = "row.names")
rownames(Exp_tsne_Mat) <- Exp_tsne_Mat[,1]
Exp_tsne_Mat <- Exp_tsne_Mat[, -1]
Exp_tsne_Mat <- merge(Exp_tsne_Mat, t(ExpMat), by = "row.names")
rownames(Exp_tsne_Mat) <- Exp_tsne_Mat[,1]
Exp_tsne_Mat <- Exp_tsne_Mat[, -1]
# Fig. 6A
Plot_Dimentional_Reduction(Exp_tsne_Mat, name_title = "UMAP visualization\nbased on gene expression", 
                           name_color = "celltype", lab_color = "Cell Type", max_color = NA,
                           flag_save_plot = TRUE, 
                           path_reduction_plot = paste0(image_folder, "UMAP"), type = ".pdf")


# ----
# 2 calculate activity --------
# 2.1 VIPER --------
# regulon
regulon_TD_Viper <- readRDS(file = "../data/regulon/regulon_TD_Viper.rds")
regulon_E_Viper <- readRDS(file = "../data/regulon/regulon_E_Viper.rds")
regulon_E_i_Viper <- readRDS(file = "../data/regulon/regulon_E_isoform_Viper.rds")

# viper.all.E <- Calculate_Activity(ExpMat, meta_info, regulon_E_Viper, method = "viper", minsize = 2)
# dim(viper.all.E) # [1] 1566   73
# write.csv(viper.all.E, file = paste0(data_folder, "TFactivity.VIPER.regulonE.1566TF_1930cell.csv"),
#           na = 'NA', fileEncoding="UTF-8")
viper.all.E <- read.csv(paste0(data_folder, "TFactivity.VIPER.regulonE.1566TF_1930cell.csv"), stringsAsFactors = FALSE,
                        row.names = 1, header = T, encoding = "UTF-8", check.names = FALSE)

# viper.all.E_i <- Calculate_Activity(ExpMat, meta_info, regulon_E_i_Viper, method = "viper", minsize = 2)
# dim(viper.all.E_i) # [1] 3237   73
# write.csv(viper.all.E_i, file = paste0(data_folder, "TFactivity.VIPER.regulonE_i.3237TF_1930cell.csv"),
#           na = 'NA', fileEncoding="UTF-8")
viper.all.E_i <- read.csv(paste0(data_folder, "TFactivity.VIPER.regulonE_i.3237TF_1930cell.csv"), stringsAsFactors = FALSE,
                          row.names = 1, header = T, encoding = "UTF-8", check.names = FALSE)
# 
# viper.all.TD <- Calculate_Activity(ExpMat, meta_info, regulon_TD_Viper, method = "viper", minsize = 2)
# dim(viper.all.TD) # [1] 1457   73
# write.csv(viper.all.TD, file = paste0(data_folder, "TFactivity.VIPER.regulonTD.1457TF_1930cell.csv"),
#           na = 'NA', fileEncoding="UTF-8")
viper.all.TD <- read.csv(paste0(data_folder, "TFactivity.VIPER.regulonTD.1457TF_1930cell.csv"), stringsAsFactors = FALSE,
                         row.names = 1, header = T, encoding = "UTF-8", check.names = FALSE)

# 2.2 calculate mean activity --------
celltype_list <- unique(meta_info$celltype)
lineage_list <- unique(meta_info$lineage)
ActMat_E_Mean_Lineage <- Calculate_Mean_Activity_by_Celltype(viper.all.E, meta_info, name_group = "lineage", celltype_list = lineage_list)
# write.csv(ActMat_E_Mean_Lineage, file = paste0(data_folder, "TFactivity.VIPER.regulonE.1566TF_1930cell.mean_by_lineage.csv"),
#           na = 'NA', fileEncoding="UTF-8")
ActMat_i_Mean_Lineage <- Calculate_Mean_Activity_by_Celltype(viper.all.E_i, meta_info, name_group = "lineage", celltype_list = lineage_list)
# write.csv(ActMat_i_Mean_Lineage, file = paste0(data_folder, "TFactivity.VIPER.regulonE_i.3237TF_1930cell.mean_by_lineage.csv"),
#           na = 'NA', fileEncoding="UTF-8")
ActMat_TD_Mean_Lineage <- Calculate_Mean_Activity_by_Celltype(viper.all.TD, meta_info, name_group = "lineage", celltype_list = lineage_list)
# write.csv(ActMat_TD_Mean_Lineage, file = paste0(data_folder, "TFactivity.VIPER.regulonTD.1457TF_1930cell.mean_by_lineage.csv"),
#           na = 'NA', fileEncoding="UTF-8")

# ----
# 3 comparison of regulon database --------
# 3.1 TF for gastrula --------
geneset <- c("ID1","PAX3","NEUROG1","PAX2",
             "PAX8","CDX1","FOXA2","SALL1",
             "SOX17","HHEX","FOXA3","KLF3")
index.12TF <- data.frame(gene = geneset, 
                         class = 
                           c(rep("Ectoderm-associated TFs",4),
                             rep("Mesoderm-associated TFs",4),
                             rep("Endoderm-associated TFs",4)))
index.12TF.anno <- index.12TF
rownames(index.12TF.anno) <- index.12TF[,1]
index.12TF.anno <- index.12TF.anno[,2,drop = FALSE]
ActMat_E_Mean_Lineage_12TF <- ActMat_E_Mean_Lineage[
  rownames(ActMat_E_Mean_Lineage)%in%geneset,]
ActMat_E_Mean_Lineage_12TF <- ActMat_E_Mean_Lineage_12TF[geneset,]
# Figure 6B.1
TFA_E_12TF <- pheatmap(ActMat_E_Mean_Lineage_12TF,
                       show_rownames = T,show_colnames = T,
                       angle_col = 45,
                       cluster_cols = F,cluster_rows = F,
                       # border_color = 'NA',
                       scale = 'row',
                       cellwidth = 13,
                       cellheight = 13,
                       annotation_row = index.12TF.anno,
                       color= colorRampPalette(c("navy", "white", "red"))(50),
                       main = 'TFA cross 3 lineage (CRED)',
                       # cutree_rows = 8
                       # display_numbers = T,number_color = 'blue4',
                       # cellwidth = 15,cellheight = 15
                       # annotation_colors = mycolors,
                       # filename = 'spearman distance.pdf'
)
Save_Pheatmap_Pdf(TFA_E_12TF, 
                  paste0(image_folder, "Activity_CRED_12TF_mean_by_lineage.pdf"), width = 10, height = 10)
png(paste0(image_folder, "Activity_CRED_12TF_mean_by_lineage.png"), units="in", width=6, height=6, res = 300)
print(TFA_E_12TF)
dev.off()

ActMat_TD_Mean_Lineage_12TF <- ActMat_TD_Mean_Lineage[
  rownames(ActMat_TD_Mean_Lineage)%in%geneset,]
ActMat_TD_Mean_Lineage_12TF <- ActMat_TD_Mean_Lineage_12TF[geneset,]
rownames(ActMat_TD_Mean_Lineage_12TF)[1] <- "ID1"
# Figure 6B.2
TFA_TD_12TF <- pheatmap(ActMat_TD_Mean_Lineage_12TF,
                      show_rownames = T,show_colnames = T,
                      angle_col = 45,
                      cluster_cols = F,cluster_rows = F,
                      # border_color = 'NA',
                      scale = 'row',
                      cellwidth = 13,
                      cellheight = 13,
                      annotation_row = index.12TF.anno,
                      color= colorRampPalette(c("navy", "white", "red"))(50),
                      main = 'TFA cross 3 lineage (TD)',
                      # cutree_rows = 8
                      # display_numbers = T,number_color = 'blue4',
                      # cellwidth = 15,cellheight = 15
                      # annotation_colors = mycolors,
                      # filename = 'spearman distance.pdf'
)
Save_Pheatmap_Pdf(TFA_TD_12TF, 
                  paste0(image_folder, "Activity_TD_12TF_mean_by_lineage.pdf"), width = 10, height = 10)
png(paste0(image_folder, "Activity_TD_12TF_mean_by_lineage.png"), units="in", width=6, height=6, res = 300)
print(TFA_TD_12TF)
dev.off()

# 3.2 TF for mesoderm differentiation --------
geneset2 <- c("MYCN","FOS","HSF1","TP53","KLF8")
index.5TF <- data.frame(gene = geneset2, 
                        class = 
                          c(rep("Li et al.(2022)",4),
                            "Chu et al.(2016)"))
index.5TF.anno <- index.5TF
rownames(index.5TF.anno) <- index.5TF[,1]
index.5TF.anno <- index.5TF.anno[,2,drop = FALSE]

ActMat_E_Mean_Lineage_5TF <- ActMat_E_Mean_Lineage[
  rownames(ActMat_E_Mean_Lineage)%in%geneset2,]
ActMat_E_Mean_Lineage_5TF <- ActMat_E_Mean_Lineage_5TF[geneset2,]

# Figure 6D.1
TFA_E_5TF <- pheatmap(ActMat_E_Mean_Lineage_5TF,
                  show_rownames = T,show_colnames = T,
                  angle_col = 45,
                  cluster_cols = F,cluster_rows = F,
                  # border_color = 'NA',
                  scale = 'row',
                  cellwidth = 18,
                  cellheight = 18,
                  annotation_row = index.5TF.anno,
                  color= colorRampPalette(c("navy", "white", "red"))(50),
                  main = 'TFA key TF in Endo (CRED)',
                  # cutree_rows = 8
                  # display_numbers = T,number_color = 'blue4',
                  # cellwidth = 15,cellheight = 15
                  # annotation_colors = mycolors,
                  # filename = 'spearman distance.pdf'
)
Save_Pheatmap_Pdf(TFA_E_5TF, 
                  paste0(image_folder, "Activity_E_5TF_mean_by_lineage.pdf"), width = 10, height = 10)
png(paste0(image_folder, "Activity_E_5TF_mean_by_lineage.png"), units="in", width=6, height=6, res = 300)
print(TFA_E_5TF)
dev.off()

ActMat_TD_Mean_Lineage_5TF <- ActMat_TD_Mean_Lineage[
  rownames(ActMat_TD_Mean_Lineage)%in%geneset2,]
ActMat_TD_Mean_Lineage_5TF <- ActMat_TD_Mean_Lineage_5TF[geneset2,
]

# Figure 6D.2
TFA_TD_5TF <- pheatmap(ActMat_TD_Mean_Lineage_5TF,
                     show_rownames = T,show_colnames = T,
                     angle_col = 45,
                     cluster_cols = F,cluster_rows = F,
                     # border_color = 'NA',
                     scale = 'row',
                     cellwidth = 18,
                     cellheight = 18,
                     annotation_row = index.5TF.anno,
                     color= colorRampPalette(c("navy", "white", "red"))(50),
                     main = 'TFA key TF in Endo (TD)',
                     # cutree_rows = 8
                     # display_numbers = T,number_color = 'blue4',
                     # cellwidth = 15,cellheight = 15
                     # annotation_colors = mycolors,
                     # filename = 'spearman distance.pdf'
)
Save_Pheatmap_Pdf(TFA_TD_5TF, 
                  paste0(image_folder, "Activity_TD_5TF_mean_by_lineage.pdf"), width = 10, height = 10)
png(paste0(image_folder, "Activity_TD_5TF_mean_by_lineage.png"), units="in", width=6, height=6, res = 300)
print(TFA_TD_5TF)
dev.off()

# 3.3 TF from TFome --------
# Ng, A. H. M. et al. A comprehensive library of human transcription factors for 
# cell fate engineering. Nat Biotechnol 39, 510–519 (2021).
TFome_65 <- read.csv(paste0(data_folder, 'TFOME_65.CSV'), header = F)

recall_TFome_E <- Calculate_recall_activity(ActMat_E_Mean_Lineage, TFome_65, regulon = "CRED", class = "calculated")
recall_TFome_TD <- Calculate_recall_activity(ActMat_TD_Mean_Lineage, TFome_65, regulon = "TD", class = "calculated")
recall_TFome_E_bg <- Calculate_recall_activity(ActMat_E_Mean_Lineage, TFome_65, regulon = "CRED background", class = "background")
recall_TFome_TD_bg <- Calculate_recall_activity(ActMat_TD_Mean_Lineage, TFome_65, regulon = "TD background", class = "background")

recall_TFome <- rbind(recall_TFome_E, recall_TFome_TD, recall_TFome_E_bg, recall_TFome_TD_bg)
write.csv(recall_TFome, file = paste0(data_folder, "CRED_vs_TD_in_TFome.csv"),
          na = 'NA', fileEncoding="UTF-8")
recall_TFome <- read.csv(paste0(data_folder, "CRED_vs_TD_in_TFome.csv"), stringsAsFactors = FALSE,
                         row.names = 1, header = T, encoding = "UTF-8")
recall_TFome$Rank.of.TFA=factor(recall_TFome$Rank.of.TFA,levels=unique(recall_TFome$Rank.of.TFA))
recall_TFome$Class=factor(recall_TFome$Class,levels=c('calculated','background'))
# Figure 6C
p2 <- ggplot(recall_TFome, aes(Rank.of.TFA,Percent.Recall,
                          group=Used.Regulome,col=Used.Regulome)) + 
  geom_step(aes(linetype=Class),linewidth=1.5)+
  theme_bw()+
  theme(
    legend.title=element_text(size = 15),
    legend.text=element_text(size = 14),
    axis.text.x = element_text(size = 14), 
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 15),
    panel.border = element_rect(linewidth = 2))+
  labs(x="Percentage Rank of TFA",y="Precent Recall")+
  scale_color_nejm()
p2
ggsave(paste0(image_folder, "TFome_recall", ".pdf"), p2, 
       width = 6, height = 4)

# 3.4 examples --------
# FOXP1
focus_gene <- "FOXP1"
focus_isoform_list <- c("TFORF0214-FOXP1", "TFORF0218-FOXP1")
focus_isoform_name_list <- c("FOXP1-2", "FOXP1-6")
act.mat <- viper.all.E_i
meta.info <- meta_info
group <- "lineage"
name_title <- paste0("Activity of ", focus_gene)
path_fig = paste0(image_folder, "Activity_", focus_gene, "_by_", group)
type <- ".png"
width = 5
height = 4
Plot_Activity_Violinplot_bycelltype <- function(act.mat, meta.info, group,
                                             focus_gene, focus_isoform_list, focus_isoform_name_list,
                                             name_title,
                                             flag_save_fig = FALSE, 
                                             path_fig = paste0("./", "Activity_", focus_gene, "_by_", group), type = ".png",
                                             width = 5, height = 4) {
  focus_df <- merge(meta.info[, group, drop = FALSE], t(act.mat)[, focus_isoform_list, drop = FALSE], by = "row.names")
  rownames(focus_df) <- focus_df$Row.names
  focus_df <- focus_df[,-1]
  colnames(focus_df) <- c(group, focus_isoform_name_list)
  df.long <- melt(focus_df, id.vars = group, variable.name = "isoform", value.name = "Activity")
  p <- ggplot(df.long, aes(x = isoform, y = Activity, fill = .data[[group]])) +
    geom_hline(yintercept = 0,
               linewidth = 1, linetype = "dashed", color = "gray") +
    geom_violin(linewidth = 1, position = "dodge", alpha = 0.5) +
    geom_boxplot(linewidth = 1, width = 0.2, position = position_dodge(0.9), notch = T) +
    theme_bw()+
    stat_compare_means(
      comparisons = list(focus_isoform_name_list),
      method = "wilcox.test",
      label = "p.signif"
    ) +
    labs(
      # x = "isoform",
      # y = "Activity",
      title = eval(name_title),
      # color = eval(lab_color)
    ) +
    theme(
      axis.line = element_line(linewidth = 1),
      axis.text.x = element_text(size = 14), 
      axis.text.y = element_text(size = 14),
      axis.title = element_text(size = 15),
      plot.title = element_text(size = 18),
      legend.title=element_text(size = 15),
      legend.text=element_text(size = 14),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_blank()
    )
  if (flag_save_fig) {
    ggsave(paste0(path_fig, type), p, 
           width = 5, height = 4)
  }
  return(p)
}
Plot_Activity_Violinplot_bycelltype(act.mat = viper.all.E_i, meta.info = meta_info, group = "lineage",
                                    focus_gene, focus_isoform_list, focus_isoform_name_list,
                                    name_title = paste0("Activity of ", focus_gene),
                                    flag_save_fig = TRUE, 
                                    path_fig = paste0(image_folder, "Activity_", focus_gene, "_by_", group, "_violinplot"), type = ".png",
                                    width = 5, height = 6)

Plot_Activity_Barplot_bycelltype <- function(act.mat, meta.info, group,
                                             focus_gene, focus_isoform_list, focus_isoform_name_list,
                                             name_title,
                                             flag_save_fig = FALSE, 
                                             path_fig = paste0("./", "Activity_", focus_gene, "_by_", group), type = ".png",
                                             width = 5, height = 4) {
  focus_df <- merge(meta.info[, group, drop = FALSE], t(act.mat)[, focus_isoform_list, drop = FALSE], by = "row.names")
  rownames(focus_df) <- focus_df$Row.names
  focus_df <- focus_df[,-1]
  colnames(focus_df) <- c(group, focus_isoform_name_list)
  df.long <- melt(focus_df, id.vars = group, variable.name = "isoform", value.name = "Activity")
  p <- ggplot(df.long, aes(x = isoform, y = Activity, fill = .data[[group]], group = .data[[group]])) +
    geom_hline(yintercept = 0,
               linewidth = 1, color = "gray") +
    # geom_bar(linewidth = 1,
    #          stat = "identity", position = position_dodge(0.9), width = 0.7) +
    geom_jitter(data = df.long, aes(x = .data[["isoform"]], y = .data[["Activity"]],
                                    fill = .data[[group]], color = .data[[group]]),
                position = position_jitterdodge(0.1), 
                size = 0.5, alpha = 0.1,
                # height = 0.02, width = 0.1
    ) +
    stat_summary(fun = mean, geom = "bar", width = 0.7, position = position_dodge(0.9),
                 fill = NA, linewidth = 1, aes(color = .data[[group]])) +
    stat_summary(fun.data = mean_se, geom = "errorbar", color = "black",
                 width = 0.25, position = position_dodge(0.9)) +
    theme_bw()+
    # stat_compare_means(
    #   comparisons = list(focus_isoform_name_list),
    #   method = "t.test",
    #   label = "p.signif"
    # ) +
    labs(
      # x = "isoform",
      # y = "Activity",
      title = eval(name_title),
      # color = eval(lab_color)
    ) +
    theme(
      axis.line = element_blank(),
      axis.text.x = element_text(size = 14), 
      axis.text.y = element_text(size = 14),
      axis.title = element_text(size = 15),
      plot.title = element_text(size = 18),
      legend.title=element_text(size = 15),
      legend.text=element_text(size = 14),
      panel.border = element_rect(linewidth = 2, fill = NA),
      panel.grid = element_blank(),
      panel.background = element_blank()
    )
  if (flag_save_fig) {
    ggsave(paste0(path_fig, type), p, 
           width = 5, height = 4)
  }
  return(p)
}

Plot_Activity_Barplot_bycelltype(act.mat = viper.all.E_i, meta.info = meta_info, group = "lineage",
                                 focus_gene, focus_isoform_list, focus_isoform_name_list,
                                 name_title = paste0("Activity of ", focus_gene),
                                 flag_save_fig = TRUE, 
                                 path_fig = paste0(image_folder, "Activity_", focus_gene, "_by_", group, "_barplot"), type = ".png",
                                 width = 5, height = 6)


# ----
# 4 Exp and Act of focus gene --------
# Fig. 6E, I, M, Fig. 7
focus_gene <- "POU5F1"
Plot_Dimentional_Reduction(Exp_tsne_Mat, name_title = paste0("Expression of ", focus_gene), 
                           name_color = focus_gene, lab_color = "Expression", max_color = NA, color_scale = TRUE,
                           flag_save_plot = TRUE, 
                           path_reduction_plot = paste0(image_folder, "Expression_", focus_gene), type = ".png")

Act_tsne_Mat <- merge(meta_info, embedding_umap, by = "row.names")
rownames(Act_tsne_Mat) <- Act_tsne_Mat[,1]
Act_tsne_Mat <- Act_tsne_Mat[, -1]
Act_tsne_Mat <- merge(Act_tsne_Mat, t(viper.all.E_i), by = "row.names")
rownames(Act_tsne_Mat) <- Act_tsne_Mat[,1]
Act_tsne_Mat <- Act_tsne_Mat[, -1]
# Fig. 6G, K, N, Fig.S7A, B
focus_gene <- "TFORF0993-KLF8"
focus_gene_name <- "KLF8-4"
focus_gene <- "TFORF0141-POU5F1"
focus_gene_name <- "POU5F1-4"
focus_gene <- "TFORF1486-NANOG"
focus_gene_name <- "NANOG-2"
Plot_Dimentional_Reduction(Act_tsne_Mat, name_title = paste0("Activity of ", focus_gene_name), 
                           name_color = focus_gene, lab_color = "Activity", max_color = NA, color_scale = TRUE,
                           flag_save_plot = TRUE, 
                           path_reduction_plot = paste0(image_folder, "Activity_", focus_gene), type = ".png")

# ----
# 5 comparison of exp and act --------
# Fig. 6J
POU5F1_list <- c("TFORF0141-POU5F1", "TFORF1340-POU5F1", "TFORF1341-POU5F1", "TFORF3477-POU5F1")
focus_isoform_list <- POU5F1_list
focus_isoform_list_name <- c("POU5F1-4", "POU5F1-1", "POU5F1-2", "POU5F1-3")
focus_gene <- "POU5F1"
cor_method <- "spearman"
Plot_Correlation_of_Exp_and_Act(ExpMat, viper.all.E_i, focus_isoform_list, focus_isoform_list_name, focus_gene, cor_method,
                                name_title = "Correlation of expression and activity", name_subtitle = paste0("isoforms of ", focus_gene),
                                flag_save_plot = TRUE, path_plot = paste0(image_folder, "cor_of_exp_and_act_isoform_", focus_gene), 
                                type = ".pdf", width = 5, height = 4)

# ----
# 6 activity heatmap of TF pair --------
name_group <- "lineage"
ActMat <- ActMat_i_Mean_Lineage
focus_gene <- "KLF8"
focus_gene <- "NANOG"
focus_isoform_list <- rownames(ActMat)[grepl(paste0(focus_gene, "$"), rownames(ActMat))]
focus_act_mat <- as.data.frame(t(ActMat[focus_isoform_list,]))
for (col in colnames(focus_act_mat)) {
  focus_act_mat[,col] <- scale(focus_act_mat[,col, drop = F])
}
colnames(focus_act_mat) <- c("KLF8-1", "KLF8-2", "KLF8-3", "KLF8-4")
colnames(focus_act_mat) <- c("NANOG-1", "NANOG-2")
focus_act_mat[,name_group] <- rownames(focus_act_mat)
focus_act_mat.long <- melt(focus_act_mat, id.vars = name_group, variable.name = "isoform", value.name = "activity")
# Fig. 6F, Fig. S7C
Plot_Heatmap_of_TF_Pair(focus_act_mat.long, name_value = "activity", name_x = name_group, name_y = "isoform",
                        name_title = paste0("Mean activity in lineage"), lab_color = "Scaled\nActivity",
                        flag_save_plot = TRUE,
                        path_plot = paste0(image_folder, "Activity_", focus_gene, "_mean_by_", name_group), type = ".pdf",
                        width = 4, height = 2)

# ----
# 7 activity correlation of TF pair --------
EOMES_list <- c("TFORF1158-EOMES", "TFORF1159-EOMES", "TFORF1160-EOMES")
HMGA2_list <- c("TFORF2849-HMGA2", "TFORF2850-HMGA2", "TFORF2851-HMGA2", "TFORF2852-HMGA2", "TFORF2853-HMGA2")

focus_gene1 <- "EOMES"
focus_gene2 <- "HMGA2"
focus_gene_list1 <- EOMES_list
focus_gene_list2 <- HMGA2_list
focus_TF_pair_list <- data.frame(
  "TF1" = sort(c(rep(focus_gene_list1, length(focus_gene_list2)))),
  "TF2" = c(rep(focus_gene_list2, length(focus_gene_list1)))
)

focus_cell_list <- list(
  "Ectoderm" = rownames(meta_info[meta_info$lineage == "Ectoderm",]),
  # "Endoderm" = rownames(meta.info[meta.info$lineage == "Endoderm",]),
  # "Mesoderm" = rownames(meta.info[meta.info$lineage == "Mesoderm",]),
  "Meso&Endo" = rownames(meta_info[meta_info$lineage %in% c("Mesoderm", "Endoderm"),]),
  "whole" = rownames(meta_info)
)

cor_method <- "pearson"
cor_df_focus <- Calculate_Correlation_of_Activity_by_Focus(act.mat = viper.all.E_i, focus_TF_pair_list, focus_cell_list,
                                                           cor_method = cor_method) 
# rename
cor_df_focus$TF1 <- rep(c("EOMES-1", "EOMES-2", "EOMES-3"),each = 5)
cor_df_focus$TF2 <- rep(c("HMGA2-1", "HMGA2-2", "HMGA2-3", "HMGA2-4", "HMGA2-5"), 3)

# Fig.7 C
df_long <- cor_df_focus
name_value <- "whole"
name_x <- "TF1"
name_y <- "TF2"
focus_gene1 <- "EOMES"
focus_gene2 <- "HMGA2"
cor_method <- "pearson"
lab_color <- "correlation"
Plot_Heatmap_of_TF_Pair(cor_df_focus, name_value, name_x, name_y, 
                        cluster_rows = TRUE, cluster_cols = TRUE,
                        name_title = paste0("Activity correlation in ", name_value, " cells (", cor_method, ")"),
                        lab_color = lab_color,
                        flag_save_plot = TRUE, 
                        path_plot = paste0(image_folder, "Activity_", focus_gene1, "_", focus_gene2, "_correlation"), type = ".png",
                        width = 5, height = 4)

# Fig.7 D
focus_gene1 <- "TFORF1159-EOMES"
focus_gene1_name <- "EOMES-2"
focus_gene2 <- "TFORF2849-HMGA2"
focus_gene2_name <- "HMGA2-1"
focus_gene2 <- "TFORF2852-HMGA2"
focus_gene2_name <- "HMGA2-4"
focus_gene2 <- "TFORF2853-HMGA2"
focus_gene2_name <- "HMGA2-5"
name_group <- "lineage"
Plot_Activity_Scatter(act.mat = viper.all.E_i, meta.info = meta_info, focus_gene1 = focus_gene1, focus_gene2 = focus_gene2, name_group = name_group,
                      name_x = focus_gene1_name, name_y = focus_gene2_name,
                      name_title = "Activity of gene pair", flag_scale = FALSE,
                      flag_save_fig = TRUE, path_fig = paste0(image_folder, "Activity_", focus_gene1_name, "_", focus_gene2_name, "_scatter_by_", name_group, ".pdf"))

# Fig.7 D
focus_gene1 <- "EOMES"
focus_gene1_name <- "EOMES-2"
focus_gene2 <- "HMGA2"
focus_gene2_name <- "HMGA2-1"
focus_gene2 <- "TFORF2852-HMGA2"
focus_gene2_name <- "HMGA2-4"
focus_gene2 <- "TFORF2853-HMGA2"
focus_gene2_name <- "HMGA2-5"
name_group <- "lineage"
Plot_Activity_Scatter(act.mat = ExpMat, meta.info = meta_info, focus_gene1 = focus_gene1, focus_gene2 = focus_gene2, name_group = name_group,
                      name_x = focus_gene1, name_y = focus_gene2,
                      name_title = "Expression of gene pair", flag_scale = FALSE,
                      flag_save_fig = TRUE, path_fig = paste0(image_folder, "Expression_", focus_gene1, "_", focus_gene2, "_scatter_by_", name_group, ".png"))


