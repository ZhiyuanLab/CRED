# '=== Header Start ==='
# Title:       Section3_Unipotent_Network
# Author:      Wanqi Li
# Date:        20240322
# Purpose:     analysis and figure for Unipotent Network mentioned in Section 3
# '=== Header End ==='

# config --------
setwd("D:/wq/zhiyuanlab/20230717-CRED/CRED_project/code/")
data_folder <- "../data/Simulation_Unipotent_Network/"
image_folder <- "../data/Simulation_Unipotent_Network/image/"

# import --------
library('viper')
library('Rtsne')
library('ggplot2')

# color palette --------
colors_pseudotime <- colorRampPalette(c("#00007F", "blue", "#007FFF", 
                                        "cyan", "#7FFF7F", "yellow", 
                                        "#FF7F00", "red", "#7F0000"))

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

Do_Dimentional_Reduction <- function(Exp.Mat, 
                                     method = "t-SNE", p = 400, seed = 5,
                                     flag_save_pos = FALSE, path_pos) {
  # do t-SNE and save position
  if (method != "t-SNE") {
    cat("ERROR: unknown method!")
    return()
  }
  set.seed(seed)
  tsne_obj <- Rtsne(t(Exp.Mat),
                    perplexity = p, check_duplicates = FALSE)
  tsne_df <- data.frame(tSNE_1 = tsne_obj$Y[, 1], 
                        tSNE_2 = tsne_obj$Y[, 2],
                        row.names=colnames(Exp.Mat))
  if (flag_save_pos) {
    write.csv(tsne_df, path_pos, row.names=TRUE)
  }
  return(tsne_df)
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
  } else {
    cat("WARNING: unknown color lab!")
  }
  
  if (flag_save_plot) {
    ggsave(paste0(path_reduction_plot, type), p.tsne, 
           width = 5, height = 4)
  }
  return(p.tsne)
}

# ----
# 0 read data --------
ExpMat <- read.csv(paste0(data_folder, "ExpressionData.csv"), 
stringsAsFactors = FALSE,
row.names = 1, header = T, encoding = "UTF-8")
dim(ExpMat) # [1]   15 2000
# meta info with annotation and prediction
meta_info <- read.csv(paste0(data_folder, "PseudoTime.csv"),
                      stringsAsFactors = FALSE,
                      row.names = 1, header = T, encoding = "UTF-8")
dim(meta_info) # [1] 2000    3
# regulon
regulon <- read.csv(paste0(data_folder, "regulon.csv"),
                    stringsAsFactors = FALSE,
                    header = T, encoding = "UTF-8")
regulon_Viper <- Format_Regulon_DoRothEA2Viper(regulon)
regulon_transited <- read.csv(paste0(data_folder, "regulon_transited.csv"),
                              stringsAsFactors = FALSE,
                              header = T, encoding = "UTF-8")
regulon_transited_Viper <- Format_Regulon_DoRothEA2Viper(regulon_transited)

# ----
# 1 calculate activity --------
ActMat <- Calculate_Activity(ExpMat, meta_info, regulon_Viper, minsize = 1)
dim(ActMat) # [1]    2 2000
write.csv(ActMat, file = paste0(data_folder, "Activity.csv"),
          na = 'NA', fileEncoding="UTF-8")
ActMat <- read.csv(paste0(data_folder, "Activity.csv"), stringsAsFactors = FALSE,
                   row.names = 1, header = T, encoding = "UTF-8", check.names = FALSE)
ActMat <- as.data.frame(t(ActMat))

ActMat_transited <- Calculate_Activity(ExpMat, meta_info, regulon_transited_Viper, minsize = 1)
dim(ActMat_transited) # [1]    2 2000
write.csv(ActMat_transited, file = paste0(data_folder, "Activity_transited.csv"),
          na = 'NA', fileEncoding="UTF-8")
ActMat_transited <- read.csv(paste0(data_folder, "Activity_transited.csv"), stringsAsFactors = FALSE,
                             row.names = 1, header = T, encoding = "UTF-8", check.names = FALSE)
ActMat_transited <- as.data.frame(t(ActMat_transited))

# ----
# 2 plot dimentional reduction --------
# do reduction and save position
tsne_df <- Do_Dimentional_Reduction(ExpMat, 
                                    method = "t-SNE", p = 300, seed = 15,
                                    flag_save_pos = TRUE, path_pos = paste0(data_folder, "tsne.csv"))
tsne_df <- read.csv(paste0(data_folder, "tsne.csv"), 
                    stringsAsFactors = FALSE,
                    row.names = 1, header = T, encoding = "UTF-8")
colnames(tsne_df) <- c("tSNE_1", "tSNE_2")

# merge expression and position
Exp_tsne_Mat <- merge(t(ExpMat), tsne_df,
                      by = 'row.names')
rownames(Exp_tsne_Mat) <- Exp_tsne_Mat[,1]
colnames(Exp_tsne_Mat)[1] <- "cell"
Exp_tsne_Mat$PseudoTime <- as.numeric(as.data.frame(matrix(unlist(strsplit(as.character(Exp_tsne_Mat$cell), '_')), ncol=2, byrow=TRUE))[,2])

# get focus time window
para_focus_time_max <- 500
Exp_tsne_Mat_focus <- Exp_tsne_Mat[Exp_tsne_Mat$PseudoTime <= para_focus_time_max,]
dim(Exp_tsne_Mat_focus) # [1] 1265   19

# plot reduction
# Fig. 3C
Plot_Dimentional_Reduction(Exp_tsne_Mat_focus, name_title = "t-SNE visualization\nbased on gene expression", 
                           name_color = "PseudoTime", lab_color = "PseudoTime", max_color = 500,
                           flag_save_plot = TRUE, path_reduction_plot = paste0(image_folder, "tsne_early"), type = ".pdf")

# plot expression on reduction
# Supp. Fig. 4
focus_gene <- "S3" # change this
focus_gene_name <- "S3"
Plot_Dimentional_Reduction(Exp_tsne_Mat_focus, name_title = paste0("Gene expression of ",focus_gene_name, "\nwith t-SNE embedding"), 
                           name_color = focus_gene, lab_color = "Expression",
                           flag_save_plot = TRUE, path_reduction_plot = paste0(image_folder, "tsne_exp_", focus_gene_name), type = ".pdf")

# merge activity and position
Act_tsne_Mat <- merge(ActMat, tsne_df,
                      by = 'row.names')
rownames(Act_tsne_Mat) <- Act_tsne_Mat[,1]
colnames(Act_tsne_Mat)[1] <- "cell"
Act_tsne_Mat$PseudoTime <- as.numeric(as.data.frame(matrix(unlist(strsplit(as.character(Act_tsne_Mat$cell), '_')), ncol=2, byrow=TRUE))[,2])
Act_tsne_Mat_focus <- Act_tsne_Mat[Act_tsne_Mat$PseudoTime <= para_focus_time_max,]
dim(Act_tsne_Mat_focus) # [1] 125   6

# plot activity on reduction
# Fig. 3D
focus_gene <- "S3" # change this
Plot_Dimentional_Reduction(Act_tsne_Mat_focus, name_title = paste0("Gene activity of ",focus_gene, "\nwith t-SNE embedding"), 
                           name_color = focus_gene, lab_color = "Activity",
                           flag_save_plot = TRUE, path_reduction_plot = paste0(image_folder, "tsne_act_", focus_gene), type = ".png")

# merge transited activity and position
Act_tsne_Mat_transited <- merge(ActMat_transited, tsne_df,
                      by = 'row.names')
rownames(Act_tsne_Mat_transited) <- Act_tsne_Mat_transited[,1]
colnames(Act_tsne_Mat_transited)[1] <- "cell"
Act_tsne_Mat_transited$PseudoTime <- as.numeric(as.data.frame(matrix(unlist(strsplit(as.character(Act_tsne_Mat_transited$cell), '_')), ncol=2, byrow=TRUE))[,2])
Act_tsne_Mat_transited_focus <- Act_tsne_Mat_transited[Act_tsne_Mat_transited$PseudoTime <= para_focus_time_max,]
dim(Act_tsne_Mat_transited_focus) # [1] 125   6

# plot activity on reduction
# Fig. 3D
focus_gene <- "S3" # change this
Plot_Dimentional_Reduction(Act_tsne_Mat_transited_focus, name_title = paste0("Gene activity of ",focus_gene, " (transited)\nwith t-SNE embedding"), 
                           name_color = focus_gene, lab_color = "Activity",
                           flag_save_plot = TRUE, path_reduction_plot = paste0(image_folder, "tsne_act_transited_", focus_gene), type = ".pdf")


