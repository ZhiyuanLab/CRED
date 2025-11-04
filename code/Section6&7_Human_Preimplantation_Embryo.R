# '=== Header Start ==='
# Title:       Section6_Human_Preimplantation_Embryo
# Author:      Wanqi Li
# Date:        20240515
# Purpose:     Analysis and figure for human preimplantation embryo mentioned in Section 6
# Source Data: Torre, D. et al. Isoform-resolved transcriptome of the human 
#              preimplantation embryo. Nat Commun 14, 6902 (2023).
# '=== Header End ==='

# config --------
setwd("E:/wq/zhiyuanlab/20230717-CRED/CRED_project/code/")
data_folder <- "../data/Dataset_Human_Preimplantation_Embryo/"
image_folder <- "../data/Dataset_Human_Preimplantation_Embryo/image/"

# import --------
library('ggplot2')
library('reshape2')

# color palette --------

# useful function --------
Plot_Exp_by_Stage <- function(Exp.Mat_i.cpm, focus_gene, flag_multi_tr = FALSE,
                              meta.info_i, name_group, name_y,
                              name_title,
                              flag_save_plot, path_fig, type, 
                              width = 5, height = 4) {
  # draw exp by developmental stage
    if (flag_multi_tr) {
      # merge different transcripts
      exp.focus <- Exp.Mat_i.cpm[focus_gene, ]
      df <- merge(meta.info_i, colSums(exp.focus), by = "row.names")
      rownames(df) <- df[,1]
      df <- df[,-1]
      colnames(df)[ncol(df)] <- focus_gene_name
    }
    else {
      exp.focus <- Exp.Mat_i.cpm[focus_gene, ]
      df <- merge(meta.info_i, exp.focus, by = "row.names")
      rownames(df) <- df[,1]
      df <- df[,-1]
      colnames(df)[ncol(df)] <- focus_gene_name
    }
  
  colnames(df) <- c(name_group, "expression")
  # df.long <- melt(df, id.vars = name_group, variable.name = "isoform", value.name = "expression")
  df_mean <- aggregate(df[, "expression"], list(df[,name_group]), mean)
  colnames(df_mean) <- c(name_group, "expression")
  # df_mean.long <- melt(df_mean, id.vars = name_group, variable.name = "isoform", value.name = "expression")
  
  p <- ggplot() +
    geom_bar(data = df_mean, aes(x = .data[[name_group]], y = .data[["expression"]], color = .data[[name_group]], group = name_group), 
             fill = "white", linewidth = 1.5,
             stat = "identity", position = position_dodge(0.9), width = 0.7) +
    geom_jitter(data = df, aes(x = .data[[name_group]], y = .data[["expression"]],
                                    fill = .data[[name_group]], color = .data[[name_group]]),
                position = position_jitterdodge(0.1), 
                size = 1.5, 
                # height = 0.02, width = 0.1
    ) +
    labs(
      x = name_group,
      y = name_y,
      title = eval(name_title),
      # color = eval(lab_color)
    ) +
    scale_fill_brewer(palette = "Spectral") +
    scale_color_brewer(palette = "Spectral") +
    # scale_fill_manual(values = c("#584B9F", "#A71B4B")) +
    # scale_color_manual(values = c("#584B9F", "#A71B4B")) +
    # scale_color_gradient2(low = "navyblue", mid = "gray", high = "firebrick") +
    theme(
      axis.text = element_text(size = 14),
      axis.text.x = element_blank(),
      axis.title = element_text(size = 15),
      plot.title = element_text(size = 18),
      plot.subtitle = element_text(size = 15),
      # legend.position = c(0.15, 0.2),
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      panel.border = element_rect(linewidth = 2, fill = NA),
      panel.background = element_blank()
    )
  if (flag_save_plot) {
    ggsave(paste0(path_fig, type), p, 
           width = width, height = height)
  }
  return(p)
}

Plot_Exp_by_Stage_multi_gene <- function(Exp.Mat_i.cpm, focus_gene, 
                                         meta.info_i, name_group, name_y,
                                         name_title,
                                         flag_save_plot, path_fig, type, 
                                         width = 5, height = 4) {
  # draw exp by developmental stage
  df <- meta.info_i
  for (focue_TF in names(focus_gene)) {
    focus_gene_name <- unlist(focus_gene[focue_TF])
    exp.focus <- Exp.Mat_i.cpm[focus_gene_name, ]
    if (length(focus_gene_name) > 1) {
      exp.focus <- data.frame(focue_TF = colSums(exp.focus))
    }
    df <- merge(df, exp.focus, by = "row.names")
    rownames(df) <- df[,1]
    df <- df[,-1]
    colnames(df)[ncol(df)] <- focue_TF
  }
  
  df.long <- melt(df, id.vars = name_group, variable.name = "isoform", value.name = "expression")
  df_mean <- aggregate(df[, names(focus_gene)], list(df[,name_group]), mean)
  colnames(df_mean) <- c(name_group, names(focus_gene))
  df_mean.long <- melt(df_mean, id.vars = name_group, variable.name = "isoform", value.name = "expression")
  
  p <- ggplot() +
    geom_bar(data = df_mean.long, aes(x = .data[[name_group]], y = .data[["expression"]], color = .data[["isoform"]], group = isoform), 
             fill = "white", linewidth = 1.5,
             stat = "identity", position = position_dodge(0.9), width = 0.7) +
    geom_jitter(data = df.long, aes(x = .data[[name_group]], y = .data[["expression"]],
                                    fill = .data[["isoform"]], color = .data[["isoform"]]),
                position = position_jitterdodge(0.1), 
                size = 1.5, 
                # height = 0.02, width = 0.1
    ) +
    labs(
      x = name_group,
      y = name_y,
      title = eval(name_title),
      # color = eval(lab_color)
    ) +
    # scale_fill_brewer(palette = "Spectral") +
    # scale_color_brewer(palette = "Spectral") +
    scale_fill_manual(values = c("#584B9F", "#A71B4B")) +
    scale_color_manual(values = c("#584B9F", "#A71B4B")) +
    # scale_color_gradient2(low = "navyblue", mid = "gray", high = "firebrick") +
    theme(
      axis.text = element_text(size = 14),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.title = element_text(size = 15),
      plot.title = element_text(size = 18),
      plot.subtitle = element_text(size = 15),
      # legend.position = c(0.15, 0.2),
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      panel.border = element_rect(linewidth = 2, fill = NA),
      panel.background = element_blank()
    )
  if (flag_save_plot) {
    ggsave(paste0(path_fig, type), p, 
           width = width, height = height)
  }
  return(p)
}

# ----
# 0 read data --------
ExpMat_i_raw <-  read.csv(paste0(data_folder, "human_embryo-transcript_counts.tsv"), stringsAsFactors = FALSE,
                           row.names = 1, header = T, sep = "\t", encoding = "UTF-8")
dim(ExpMat_i_raw) # [1] 213873     73
meta_info_i <- read.csv(paste0(data_folder, "meta.info.transcript.csv"), stringsAsFactors = FALSE,
                        row.names = 1, header = T, encoding = "UTF-8")
meta_info_i$developmental_stage <- factor(meta_info_i$developmental_stage,
                                          levels =c("1C", "2C", "4C", "8C", "morula", "blastocyst"))
# gene name reference
gene_reference <- read.csv(paste0(data_folder, "gene_fpkm_group.xls"), stringsAsFactors = FALSE,
                           header = T, sep = "\t", encoding = "UTF-8")

# ----
# 1 normalize data --------
# change to CPM
ExpMat_i_cpm <- apply(ExpMat_i_raw,2,function(x){
  x/sum(x) * 10^6
})
write.csv(ExpMat_i_cpm, file = paste0(data_folder, "ExpMat.cpm.73_samples.213873_transcripts.csv"),
          na = 'NA', fileEncoding="UTF-8")

# ----
# 2 exp pattern during development --------

# Fig. 6L
focus_gene <- c("ENST00000471529", "ENST00000513407")
focus_gene_name <- "POU5F1-1" # TFORF1340-POU5F1
name_group <- "developmental_stage"
Plot_Exp_by_Stage(ExpMat_i_cpm, focus_gene, flag_multi_tr = TRUE,
                  meta_info_i, name_group, name_y = "Expression (CPM)",
                  name_title = paste("Expression of", focus_gene_name),
                  flag_save_plot = TRUE, path_fig = paste0(image_folder, "expression_", focus_gene_name), type = ".png", 
                  width = 5, height = 4)

# fig. 6H
focus_gene_name <- "KLF8"
focus_gene <- list(
  "KLF8-2" = c("ENST00000468660"),
  "KLF8-3" = c("ENST00000640927")
)
# Fig. 6O
focus_gene_name <- "NANOG"
focus_gene <- list(
  "NANOG-1" = c("ENST00000229307"),
  "NANOG-2" = c("ENST00000526286")
)
name_group <- "developmental_stage"
Plot_Exp_by_Stage_multi_gene(ExpMat_i_cpm, focus_gene,
                             meta_info_i, name_group, name_y = "Expression (CPM)",
                             name_title = paste("Expression of", focus_gene_name),
                             flag_save_plot = TRUE, path_fig = paste0(image_folder, "expression_", focus_gene_name), type = ".pdf", 
                             width = 8, height = 4)


