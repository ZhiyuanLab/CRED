# '=== Header Start ==='
# Title:       Section6_Comparison_of_Regulome
# Author:      Wanqi Li
# Date:        20240323
# Purpose:     analysis and figure for Comparison of different regulome database mentioned in Section 6
# Source Data: Joung, J. et al. A transcription factor atlas of directed 
#              differentiation. Cell 186, 209-229.e26 (2023).
# '=== Header End ==='

# config --------
setwd("E:/wq/zhiyuanlab/20230717-CRED/CRED_project/code/")
data_folder <- "../data/Dataset_Combinatorial_Differentiation/"
image_folder <- "../data/Dataset_Combinatorial_Differentiation/image/"


# import --------
library('viper')
library('parallel')
library('ggplot2')

# color palette --------
colors_regulon = c('#3b4992', '#ee0000', '#008b45')
names(colors_regulon) = c("regulon_E", "regulon_TD", "regulon_TDE")

# useful function --------
Viper_ttest_Enhanced <- function (eset, pheno, control_group,
                                  regulon, dnull = NULL, pleiotropy = FALSE, nes = TRUE, 
                                  method = c("none", "scale", "rank", "mad", "ttest", "ttest2"), 
                                  bootstraps = 0, 
                                  minsize = 25, adaptive.size = FALSE, eset.filter = TRUE, 
                                  mvws = 1, pleiotropyArgs = list(regulators = 0.05, shadow = 0.05, 
                                                                  targets = 10, penalty = 20, method = "adaptive"), cores = 1, 
                                  verbose = TRUE) 
{
  method <- match.arg(method)
  pdata <- NULL
  if (is(eset, "viperSignature")) {
    dnull <- eset$nullmodel
    eset <- eset$signature
    method = "none"
    if (bootstraps > 0) {
      bootstraps <- 0
      warning("Using a null model, bootstraps iterations are ignored.", 
              call. = FALSE)
    }
  }
  if (pleiotropy & bootstraps > 0) {
    bootstraps <- 0
    warning("Using pleiotropic correction, bootstraps iterations are ignored.", 
            call. = FALSE)
  }
  if (is(eset, "ExpressionSet")) {
    ExpSet <- eset
    pdata <- phenoData(eset)
    eset <- exprs(eset)
  } else {
    warning("test2 method need ExpressionSet input", 
            call. = FALSE)
  }
  if (is.null(nrow(eset))) 
    eset <- matrix(eset, length(eset), 1, dimnames = list(names(eset), 
                                                          NULL))
  if (verbose) 
    message("\nComputing the association scores")
  if (names(regulon[[1]])[1] == "tfmode") 
    regulon <- list(regulon = regulon)
  if (bootstraps > 0) {
    return(bootstrapViper(eset = eset, regulon = regulon, 
                          nes = nes, bootstraps = bootstraps, eset.filter = eset.filter, 
                          adaptive.size = adaptive.size, minsize = minsize, 
                          cores = cores, verbose = verbose))
  }
  cores1 <- 1
  if (length(regulon) > 1) {
    cores1 <- min(cores, length(regulon))
    cores <- 1
    nes <- TRUE
  }
  if (cores > 1 | cores1 > 1) 
    verbose <- FALSE
  # switch(method, scale = {
  #   tt <- t(scale(t(eset)))
  # }, rank = {
  #   tt <- t(apply(eset, 1, rank)) * punif(length(eset), 
  #                                         -0.1, 0.1)
  # }, mad = {
  #   tt <- t(apply(eset, 1, function(x) (x - median(x))/mad(x)))
  # }, ttest = {
  #   tt <- sapply(1:ncol(eset), function(i, eset) rowTtest(eset[, 
  #                                                              i] - eset[, -i])$statistic, eset = eset)
  #   colnames(tt) <- colnames(eset)
  #   rownames(tt) <- rownames(eset)
  # }, none = {
  #   tt <- eset
  # })
  # method = ttest2:
  control <- pData(ExpSet)[[pheno]] %in% control_group
  control.exp <- exprs(ExpSet)[,control]
  eset <-  as.matrix(exprs(ExpSet)[,!control])
  tt <- sapply(1:ncol(eset), function(i, eset) {
    t <- rowTtest(eset[,i] - control.exp)$statistic
    return(t)
  }, eset = eset)
  colnames(tt) <- colnames(eset)
  rownames(tt) <- rownames(eset)
  # end method
  
  if (verbose) 
    message("Computing regulons enrichment with aREA")
  res <- mclapply(regulon, function(regulon, dnull, pleiotropy, 
                                    nes, tt, eset.filter, adaptive.size, cores, pleiotropyArgs, 
                                    verbose) {
    if (eset.filter) {
      tmp <- c(names(regulon), unlist(lapply(regulon, 
                                             function(x) names(x$tfmode)), use.names = FALSE))
      tt <- tt[rownames(tt) %in% unique(tmp), ]
    }
    regulon <- lapply(regulon, function(x, genes) {
      filtro <- names(x$tfmode) %in% genes
      x$tfmode <- x$tfmode[filtro]
      if (length(x$likelihood) == length(filtro)) 
        x$likelihood <- x$likelihood[filtro]
      return(x)
    }, genes = rownames(tt))
    if (adaptive.size) 
      regulon <- regulon[sapply(regulon, function(x) {
        sum((x$likelihood/max(x$likelihood))^2)
      }) >= minsize]
    else regulon <- regulon[sapply(regulon, function(x) length(x$tfmode)) >= 
                              minsize]
    es <- aREA(tt, regulon, cores = cores, minsize = 0, 
               verbose = verbose)
    if (!nes) {
      if (pleiotropy) 
        warning("No pleiotropy correction implemented when raw es is returned.", 
                call. = FALSE)
      return(es$es)
    }
    if (is.null(dnull)) 
      nes <- es$nes
    else {
      if (verbose) 
        message("\nEstimating NES with null model")
      tmp <- aREA(dnull, regulon, cores = cores, minsize = 0, 
                  verbose = verbose)$es
      if (ncol(tmp) > 499) {
        nes <- t(sapply(1:nrow(tmp), function(i, tmp, 
                                              es) {
          aecdf(tmp[i, ], symmetric = TRUE)(es[i, ])$nes
        }, tmp = tmp, es = es$es))
        rownames(nes) <- rownames(es$nes)
      }
      else {
        nes <- es$es/sqrt(frvarna(tmp)[, 1])
      }
    }
    if (pleiotropy) {
      pb <- NULL
      if (verbose) {
        message("\nComputing pleiotropy for ", ncol(nes), 
                " samples.")
        message("\nProcess started at ", date())
      }
      if (cores > 1) {
        nes <- mclapply(1:ncol(nes), function(i, ss, 
                                              nes, regulon, args, dnull) {
          nes <- nes[, i]
          sreg <- shadowRegulon(ss[, i], nes, regulon, 
                                regulators = args[[1]], shadow = args[[2]], 
                                targets = args[[3]], penalty = args[[4]], 
                                method = args[[5]])
          if (!is.null(sreg)) {
            if (is.null(dnull)) 
              tmp <- aREA(ss[, i], sreg, minsize = 5, 
                          cores = 1)$nes[, 1]
            else {
              tmp <- aREA(cbind(ss[, i], dnull), sreg, 
                          minsize = 5, cores = 1)$es
              tmp <- apply(tmp, 1, function(x) aecdf(x[-1], 
                                                     symmetric = TRUE)(x[1])$nes)
            }
            nes[match(names(tmp), names(nes))] <- tmp
          }
          return(nes)
        }, ss = tt, nes = nes, regulon = regulon, args = pleiotropyArgs, 
        dnull = dnull, mc.cores = cores)
        nes <- sapply(nes, function(x) x)
      }
      else {
        if (verbose) 
          pb <- txtProgressBar(max = ncol(nes), style = 3)
        nes <- sapply(1:ncol(nes), function(i, ss, nes, 
                                            regulon, args, dnull, pb) {
          nes <- nes[, i]
          sreg <- shadowRegulon(ss[, i], nes, regulon, 
                                regulators = args[[1]], shadow = args[[2]], 
                                targets = args[[3]], penalty = args[[4]], 
                                method = args[[5]])
          if (!is.null(sreg)) {
            if (is.null(dnull)) 
              tmp <- aREA(ss[, i], sreg, minsize = 5)$nes[, 
                                                          1]
            else {
              tmp <- aREA(cbind(ss[, i], dnull), sreg, 
                          minsize = 5)$es
              tmp <- apply(tmp, 1, function(x) aecdf(x[-1], 
                                                     symmetric = TRUE)(x[1])$nes)
            }
            nes[match(names(tmp), names(nes))] <- tmp
          }
          if (is(pb, "txtProgressBar")) 
            setTxtProgressBar(pb, i)
          return(nes)
        }, ss = tt, nes = nes, regulon = regulon, args = pleiotropyArgs, 
        dnull = dnull, pb = pb)
      }
      if (verbose) 
        message("\nProcess ended at ", date(), "\n")
      if (is.null(nrow(nes))) 
        nes <- matrix(nes, length(nes), 1, dimnames = list(names(nes), 
                                                           NULL))
      colnames(nes) <- colnames(eset)
      return(nes)
    }
    return(nes)
  }, dnull = dnull, pleiotropy = pleiotropy, nes = nes, tt = tt, 
  eset.filter = eset.filter, adaptive.size = adaptive.size, 
  cores = cores, pleiotropyArgs = pleiotropyArgs, verbose = verbose, 
  mc.cores = cores1)
  if (length(res) == 1) 
    nes <- res[[1]]
  else {
    genes <- unique(unlist(lapply(res, rownames), use.names = FALSE))
    nes <- sapply(res, function(x, genes) as.vector(x[match(genes, 
                                                            rownames(x)), ]), genes = genes)
    nes[is.na(nes)] <- 0
    if (length(mvws) == 1) {
      ws <- abs(nes)^mvws
    }
    else {
      ws <- sigT(abs(nes), mvws[2], mvws[1])
    }
    nes <- matrix(rowSums(nes * ws)/rowSums(ws), length(genes), 
                  ncol(res[[1]]), dimnames = list(genes, colnames(res[[1]])))
  }
  if (is.null(pdata)) 
    return(nes)
  # return(ExpressionSet(assayData = nes, phenoData = pdata))
  return(nes)
}

Calculate_Activity <- function(ExpMat.logcpm, meta.info, regulon, method = "viper", 
                               minsize = 4, eset.filter = FALSE, cores = 1, verbose = TRUE,
                               pheno, control_group) {
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
                       eset.filter = eset.filter, cores = cores, verbose = verbose)
    activity <- viper.all@assayData[["exprs"]]
  } else if (method == "viper_enhanced") {
    activity <- Viper_ttest_Enhanced(ExpSet, pheno, control_group, regulon, method = "ttest", minsize = minsize,
                       eset.filter = eset.filter, cores = cores, verbose = verbose)
  } else {
    cat("WARNING: unknown method!")
  }
  return(activity)
}

Calculate_Activity_and_Rank_of_TF <- function(meta_info, ExpMat_comb_logcpm, ExpMat_ctrl_logcpm,
                                              regulon, regulon_type) {
  # use to save activity rank of focus TF
  rank_by_group <- data.frame(matrix(nrow = 0, ncol = 6))
  colnames(rank_by_group) <- c("sample.name", "combination", "TF", "rank", "percentage", "regulon_type")
  for (i in 1:nrow(meta_info)) {
    # get info of each combination
    TF_list = unlist(strsplit(meta_info[i, "combination"], split = "\\+"))
    combination = meta_info[i, "combination"]
    focus_combination = meta_info[i, "sample.name"]
    if ("GFP" %in% TF_list) {
      TF_list = TF_list[1]
      combination = TF_list[1]
    }
    
    # focus gene with top 80% exp
    focus_gene_comb <- rownames(ExpMat_comb_logcpm[ExpMat_comb_logcpm[, focus_combination] != 0,])
    focus_gene_comb <- rownames(ExpMat_comb_logcpm[order(ExpMat_comb_logcpm[focus_gene_comb, focus_combination], 
                                                         decreasing = TRUE),])[1:as.integer(0.8*length(focus_gene_comb))]
    focus_gene_ctrl <- rownames(subset(ExpMat_ctrl_logcpm, rowSums(ExpMat_ctrl_logcpm) != 0))
    focus_gene <- intersect(focus_gene_comb, focus_gene_ctrl)
    # merge with control
    ExpMat_all_logcpm <- ExpMat_ctrl_logcpm[intersect(rownames(ExpMat_ctrl_logcpm),
                                                      focus_gene), ]
    ExpMat_all_logcpm <- as.data.frame(ExpMat_all_logcpm)
    ExpMat_all_logcpm[, focus_combination] <- ExpMat_comb_logcpm[rownames(ExpMat_all_logcpm), focus_combination]
    
    # generate meta info
    meta_data <- data.frame(sample.name = colnames(ExpMat_all_logcpm),
                            type = c(rep("GFP", ncol(ExpMat_ctrl_logcpm)), 
                                     rep("combination", ncol(ExpMat_all_logcpm)-ncol(ExpMat_ctrl_logcpm))))
    row.names(meta_data) <- meta_data$sample.name
    # calculate activity by enhanced viper
    activity <- Calculate_Activity(ExpMat_all_logcpm, meta_data, regulon, method = "viper_enhanced",
                                   pheno = "type", control_group = "GFP",
                                   verbose = FALSE)
    # save TF rank
    for (TF in TF_list) {
      if (!TF %in% rownames(activity)) { 
        # regulon_E: "PAX5, "MSGN1", "HNF4A", "FERD3L", "PTF1A"
        # regulon_TD: "CDX1", "MSGN1", "FERD3L"
        # regulon_TDE: "MSGN1", "FERD3L"
        cat(TF, "not in activity\n")
        next
      }
      tmp.df <- data.frame(sample.name = focus_combination,
                           combination = combination,
                           TF = TF,
                           rank = rank(-activity[,])[TF],
                           percentage = rank(-activity[,])[TF] * 100 / length(activity),
                           regulon_type = regulon_type)
      rank_by_group <- rbind(rank_by_group, tmp.df)
    }
  }
  return(rank_by_group)
}

# ----
# 0 read data --------
ExpMat_comb <- read.csv(paste0(data_folder, "TFAtlas_combinatorial_bulk.csv"), 
                        stringsAsFactors = FALSE,
                        header = T, row.names = 1, encoding = "UTF-8")
dim(ExpMat_comb) # [1] 37302    46
ExpMat_ctrl <- read.csv(paste0(data_folder, "TFAtlas_GFP_bulk.csv"), 
                        stringsAsFactors = FALSE,
                        header = T, row.names = 1, encoding = "UTF-8")
dim(ExpMat_ctrl) # [1] 37528     3
# meta info
meta_info <- read.csv(paste0(data_folder, "meta_info_bulk.csv"), stringsAsFactors = FALSE,
                      header = T, row.names = 1, encoding = "UTF-8")
dim(meta_info) # [1] 46  2
# regulon
regulon_E_Viper <- readRDS(file = "../data/regulon/regulon_E_Viper.rds")
regulon_TD_Viper <- readRDS(file = "../data/regulon/regulon_TD_Viper.rds")
regulon_TDE_Viper <- readRDS(file = "../data/regulon/regulon_TDE_Viper.rds")

# ----
# 1 normalize data --------
# combinatorial
ExpMat_comb_cpm <- apply(ExpMat_comb,2,function(x){
  x/sum(x) * 10^6
})
ExpMat_comb_logcpm <- log2(ExpMat_comb_cpm + 1)
# control
ExpMat_ctrl_cpm <- apply(ExpMat_ctrl,2,function(x){
  x/sum(x) * 10^6
})
ExpMat_ctrl_logcpm <- log2(ExpMat_ctrl_cpm + 1)

# ----
# 2 calculate activity by group and save TF rank --------
rank_by_group <- Calculate_Activity_and_Rank_of_TF(meta_info, ExpMat_comb_logcpm, ExpMat_ctrl_logcpm,
                                                   regulon_E_Viper, regulon_type = "regulon_E")
rank_by_group <- rbind(rank_by_group,
                       Calculate_Activity_and_Rank_of_TF(meta_info, ExpMat_comb_logcpm, ExpMat_ctrl_logcpm,
                                                         regulon_TD_Viper, regulon_type = "regulon_TD"))
rank_by_group <- rbind(rank_by_group,
                       Calculate_Activity_and_Rank_of_TF(meta_info, ExpMat_comb_logcpm, ExpMat_ctrl_logcpm,
                                                         regulon_TDE_Viper, regulon_type = "regulon_TDE"))
write.csv(rank_by_group,
          file = paste0(data_folder, "activity_rank_of_focus_TF.csv"),
          row.names = TRUE, na = 'NA', fileEncoding="UTF-8")
rank_by_group <- read.csv(paste0(data_folder, "activity_rank_of_focus_TF.csv"), 
                          stringsAsFactors = FALSE,
                          header = T, row.names = 1, encoding = "UTF-8")

# ----
# 3 calculate recall by percentage step --------
recall_by_group <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(recall_by_group) <- c("percentage", "num", "type")
for (i in seq(0, max(rank_by_group$percentage), 0.5)) {
  for (type in unique(rank_by_group$regulon_type)) {
    num <- nrow(rank_by_group[rank_by_group$regulon_type == type
                                        & rank_by_group$percentage <= i,])
    num <- num / nrow(rank_by_group[rank_by_group$regulon_type == type,])
    tmp.df <- data.frame(percentage = i,
                         num = num,
                         type = type)
    recall_by_group <- rbind(recall_by_group, tmp.df)
  }
}

write.csv(recall_by_group,
          file = paste0(data_folder, "activity_recall_of_focus_TF.csv"),
          row.names = TRUE, na = 'NA', fileEncoding="UTF-8")
recall_by_group <- read.csv(paste0(data_folder, "activity_recall_of_focus_TF.csv"), 
                            stringsAsFactors = FALSE,
                            header = T, row.names = 1, encoding = "UTF-8")
# Fig. 4D
p <- ggplot(recall_by_group, 
            aes(x = percentage, y = num, 
                group = type, color = type)) +
  geom_step(size = 1.5) +
  labs(title = "recalled TF in top activity",
       x = "Activity rank percentage",
       y = "Recalled gene ratio") +
  # scale_color_aaas(name="Database",
  #                  breaks=c("regulon_E", "regulon_TD", "regulon_TDE"),
  #                  labels=c("regulon_E", "regulon_TD", "regulon_TDE")) +
  scale_color_manual(name = "Database",
                     values = colors_regulon,
                     labels = names(colors_regulon)) +
  theme_test() +
  theme(legend.justification=c(1.05,-0.1), 
        legend.position=c(1,0),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        plot.title = element_text(size = 18),
        # legend.position = c(0.15, 0.2),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        panel.border = element_rect(linewidth = 2),
        panel.background = element_blank(),)
ggsave(paste0(image_folder, "recalled_TF_ratio.png"), 
       p, width = 4, height = 4, dpi = 300)









