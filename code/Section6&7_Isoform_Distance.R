# '=== Header Start ==='
# Title:       Section6_Isoform_Distance
# Author:      Wanqi Li
# Date:        20240913
# Purpose:     Analysis and figure for isoform distance mentioned in Section 6
# '=== Header End ==='

# config --------
setwd("E:/wq/zhiyuanlab/20230717-CRED/CRED_project/code/")
data_folder <- "../data/Section5&6_Isoform/"
image_folder <- "../data/Section5&6_Isoform/image/"

# import --------
library('reshape2')
library('ggplot2')
library('ggtree')
library('aplot')

# color palette --------

# useful function --------
Plot_Heatmap_of_TF_Pair_Symmetry <- function(df_long, name_value, name_x, name_y,
                                    flag_text = FALSE, name_text,
                                    cluster_rows = FALSE, cluster_cols = FALSE,
                                    label_node = FALSE, rotate_node = NA, color_type = NA,
                                    name_title = "", lab_color,
                                    flag_save_plot = FALSE, path_plot, type,
                                    width = 5, height = 4) {
  
  # 这种行列一致的热图聚类需要把其中一个方向翻转一下，达到对称的效果
  # 所以需要反向定义横向坐标的次序
  # 翻转树没找到好用的方法，这里是先打印每个节点的编号，然后手动翻转
  df_long <- df_long[, c(name_x, name_y, name_value)]
  df_wide <- spread(df_long, key = name_y, value = name_value)
  rownames(df_wide) <- df_wide[,1]
  df_wide <- df_wide[,-1]
  if (cluster_cols) {
    bycol <- hclust(dist(df_wide))
    df_wide <- df_wide[bycol$order,] 
    v <- ggtree(bycol,layout = "rectangular",branch.length = "none", size = 0.5)+layout_dendrogram() + # 绘制行聚类树
      geom_aline(linetype="solid", linewidth = 0.5) + # 补齐末端
      labs(title = name_title) +
      theme(
        plot.title = element_text(size = 18, hjust = 0.25)
      )
    if (label_node) {
      v <- v +
        geom_text(aes(label = node))
    }
    if (!is.na(rotate_node)) {
      v <- ggtree::rotate(v, rotate_node)
    }
  }
  if (cluster_rows) {
    byrow <- hclust(dist(t(df_wide)))
    df_wide <- df_wide[,byrow$order]
    h <- ggtree(byrow,layout = "rectangular",branch.length = "none", size = 0.5) + # 绘制列聚类树
      geom_aline(linetype="solid", linewidth = 0.5) # 补齐末端
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
    scale_y_discrete(position = c("right")) + 
    labs(
      # title = eval(name_title),
      color = eval(lab_color)
    ) +
    # scale_color_gradient2(low = "navyblue", mid = "gray", high = "firebrick") +
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
  
  if (is.na(color_type) || color_type == "distance") {
    p <- p +
      scale_fill_gradient2(low = "navyblue", mid = "firebrick", high = "white", midpoint = 0,
                           # limit = c(-1,1),
                           name = lab_color) 
  } else if (color_type == "similarity") {
    p <- p +
      scale_fill_gradient2(low = "navyblue", mid = "white", high = "firebrick", midpoint = 0.3,
                           limit = c(0,1),
                           name = lab_color)
  }
  
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
  
  if (cluster_rows) {
    p <- p %>% insert_left(h,width = 0.3)
  }
  
  if (cluster_cols) {
    p <- p %>% insert_top(v,height = 0.2)
  }
  
  if (flag_save_plot) {
    ggsave(paste(path_plot, type, sep = ""), p, 
           width = width, height = height)
  }
  return(p)
}

# ----
# 1 sequence distance --------
# calculate using Matlab, plot with following code
# three level: transcript, CDS, animo acid
# 1.1 transcript --------
focus_gene <- "HMGA2"
dis_mat <- read.csv(paste0(data_folder, "sequence_distance/", "HMGA2transcript.csv"), stringsAsFactors = FALSE,
                    row.names = 1, header = T, encoding = "UTF-8", check.names = F)
dis_mat$isoform1 <- rownames(dis_mat)
dis_mat_focus <- melt(dis_mat, id.vars = "isoform1", variable.name = "isoform2", value.name = "dist")
rotate_node <- 6

focus_gene <- "EOMES"
dis_mat <- read.csv(paste0(data_folder, "sequence_distance/", "EOMEStranscript.csv"), stringsAsFactors = FALSE,
                    row.names = 1, header = T, encoding = "UTF-8", check.names = F)
dis_mat$isoform1 <- rownames(dis_mat)
dis_mat_focus <- melt(dis_mat, id.vars = "isoform1", variable.name = "isoform2", value.name = "dist")
rotate_node <- 4

focus_gene <- "JDP2"
dis_mat <- read.csv(paste0(data_folder, "sequence_distance/", "JDP2transcript.csv"), stringsAsFactors = FALSE,
                    row.names = 1, header = T, encoding = "UTF-8", check.names = F)
dis_mat$isoform1 <- rownames(dis_mat)
dis_mat_focus <- melt(dis_mat, id.vars = "isoform1", variable.name = "isoform2", value.name = "dist")
rotate_node <- 5

# Fig.7 H
df_long <- dis_mat_focus
name_value <- "dist"
name_x <- "isoform1"
name_y <- "isoform2"
lab_color <- "distance"
Plot_Heatmap_of_TF_Pair_Symmetry(df_long, name_value, name_x, name_y, 
                                 cluster_rows = TRUE, cluster_cols = TRUE,
                                 label_node = F, rotate_node = rotate_node,
                                 name_title = paste0("isoform distance of ", focus_gene),
                                 lab_color = lab_color,
                                 flag_save_plot = TRUE, 
                                 path_plot = paste0(image_folder, "Distance_transcript_", focus_gene), type = ".pdf",
                                 width = 5, height = 4)

# 1.2 CDS --------
focus_gene <- "HMGA2"
dis_mat <- read.csv(paste0(data_folder, "sequence_distance/", "ntHMGA2.csv"), stringsAsFactors = FALSE,
                    row.names = 1, header = T, encoding = "UTF-8", check.names = F)
dis_mat$isoform1 <- rownames(dis_mat)
dis_mat_focus <- melt(dis_mat, id.vars = "isoform1", variable.name = "isoform2", value.name = "dist")
rotate_node <- 6


focus_gene <- "EOMES"
dis_mat <- read.csv(paste0(data_folder, "sequence_distance/", "ntEOMES.csv"), stringsAsFactors = FALSE,
                    row.names = 1, header = T, encoding = "UTF-8", check.names = F)
dis_mat$isoform1 <- rownames(dis_mat)
dis_mat_focus <- melt(dis_mat, id.vars = "isoform1", variable.name = "isoform2", value.name = "dist")
rotate_node <- 4


focus_gene <- "JDP2"
dis_mat <- read.csv(paste0(data_folder, "sequence_distance/", "ntJDP2.csv"), stringsAsFactors = FALSE,
                    row.names = 1, header = T, encoding = "UTF-8", check.names = F)
dis_mat$isoform1 <- rownames(dis_mat)
dis_mat_focus <- melt(dis_mat, id.vars = "isoform1", variable.name = "isoform2", value.name = "dist")
rotate_node <- 3


# Fig.7 H
df_long <- dis_mat_focus
name_value <- "dist"
name_x <- "isoform1"
name_y <- "isoform2"
lab_color <- "distance"
Plot_Heatmap_of_TF_Pair_Symmetry(df_long, name_value, name_x, name_y, 
                                 cluster_rows = TRUE, cluster_cols = TRUE,
                                 label_node = F, rotate_node = rotate_node,
                                 name_title = paste0("isoform distance of ", focus_gene),
                                 lab_color = lab_color,
                                 flag_save_plot = TRUE, 
                                 path_plot = paste0(image_folder, "Distance_CDS_", focus_gene), type = ".pdf",
                                 width = 5, height = 4)

# 1.3 amino acid -------
focus_gene <- "HMGA2"
dis_mat <- read.csv(paste0(data_folder, "sequence_distance/", "aaHMGA2.csv"), stringsAsFactors = FALSE,
                    row.names = 1, header = T, encoding = "UTF-8", check.names = F)
dis_mat$isoform1 <- rownames(dis_mat)
dis_mat_focus <- melt(dis_mat, id.vars = "isoform1", variable.name = "isoform2", value.name = "dist")
rotate_node <- 6

focus_gene <- "EOMES"
dis_mat <- read.csv(paste0(data_folder, "sequence_distance/", "aaEOMES.csv"), stringsAsFactors = FALSE,
                    row.names = 1, header = T, encoding = "UTF-8", check.names = F)
dis_mat$isoform1 <- rownames(dis_mat)
dis_mat_focus <- melt(dis_mat, id.vars = "isoform1", variable.name = "isoform2", value.name = "dist")
rotate_node <- 4

focus_gene <- "JDP2"
dis_mat <- read.csv(paste0(data_folder, "sequence_distance/", "aaJDP2.csv"), stringsAsFactors = FALSE,
                    row.names = 1, header = T, encoding = "UTF-8", check.names = F)
dis_mat$isoform1 <- rownames(dis_mat)
dis_mat_focus <- melt(dis_mat, id.vars = "isoform1", variable.name = "isoform2", value.name = "dist")
rotate_node <- 3

# Fig.7 H
df_long <- dis_mat_focus
name_value <- "dist"
name_x <- "isoform1"
name_y <- "isoform2"
lab_color <- "distance"
Plot_Heatmap_of_TF_Pair_Symmetry(df_long, name_value, name_x, name_y, 
                                 cluster_rows = TRUE, cluster_cols = TRUE,
                                 label_node = F, rotate_node = rotate_node,
                                 name_title = paste0("isoform distance of ", focus_gene),
                                 lab_color = lab_color,
                                 flag_save_plot = TRUE, 
                                 path_plot = paste0(image_folder, "Distance_AA_", focus_gene), type = ".png",
                                 width = 5, height = 4)

# --------
# 2 regulatory distance --------
# 2.1 Jaccard matrix --------
focus_gene <- "HMGA2"
dis_mat <- read.csv(paste0(data_folder, "regulatory_distance/", "Jaccard_Mat_HMGA2.csv"), stringsAsFactors = FALSE,
                    row.names = 1, header = T, encoding = "UTF-8", check.names = F)
dis_mat$isoform1 <- rownames(dis_mat)
dis_mat_focus <- melt(dis_mat, id.vars = "isoform1", variable.name = "isoform2", value.name = "dist")
rotate_node <- 6

focus_gene <- "EOMES"
dis_mat <- read.csv(paste0(data_folder, "regulatory_distance/", "Jaccard_Mat_EOMES.csv"), stringsAsFactors = FALSE,
                    row.names = 1, header = T, encoding = "UTF-8", check.names = F)
dis_mat$isoform1 <- rownames(dis_mat)
dis_mat_focus <- melt(dis_mat, id.vars = "isoform1", variable.name = "isoform2", value.name = "dist")
rotate_node <- 4


focus_gene <- "JDP2"
dis_mat <- read.csv(paste0(data_folder, "regulatory_distance/", "Jaccard_Mat_JDP2.csv"), stringsAsFactors = FALSE,
                    row.names = 1, header = T, encoding = "UTF-8", check.names = F)
dis_mat$isoform1 <- rownames(dis_mat)
dis_mat_focus <- melt(dis_mat, id.vars = "isoform1", variable.name = "isoform2", value.name = "dist")
rotate_node <- 3

# Fig.7 H
df_long <- dis_mat_focus
name_value <- "dist"
name_x <- "isoform1"
name_y <- "isoform2"
lab_color <- "similarity"
Plot_Heatmap_of_TF_Pair_Symmetry(df_long, name_value, name_x, name_y, 
                                 cluster_rows = TRUE, cluster_cols = TRUE,
                                 label_node = F, rotate_node = rotate_node, color_type = "similarity",
                                 name_title = paste0("isoform distance of ", focus_gene),
                                 lab_color = lab_color,
                                 flag_save_plot = TRUE, 
                                 path_plot = paste0(image_folder, "Distance_regulon_", focus_gene), type = ".pdf",
                                 width = 5, height = 4)





