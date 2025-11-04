# '=== Header Start ==='
# Title:       Section4_Cluster4_GO_Analysis
# Author:      Wanqi Li
# Date:        20250309
# Purpose:     Analysis and figure for GO enrichment of cluster 4
# '=== Header End ==='

# config --------
setwd("E:/wq/zhiyuanlab/20230717-CRED/CRED_project/code/")
data_folder <- "../data/Dataset_TFAtlas/isoform/"
image_folder <- "../data/Dataset_TFAtlas/image/isoform/"

# import --------
library('clusterProfiler')
library('org.Hs.eg.db')
library('ggplot2')

# ----
# 0 read data --------
Jaccard.cluster <- read.csv(paste0(data_folder, "Jaccard_index_isoform_VS_isoform_6cluster.csv"), 
                            stringsAsFactors = FALSE,
                            header = T, row.names = 1, encoding = "UTF-8")
dim(Jaccard.cluster) # [1] 3250 3252
unique.interaction.cluster <- read.csv(paste0(data_folder, "../../Section5&6_Isoform/unique_interaction_isoform_byTF_clusterinfo.csv"), 
                            stringsAsFactors = FALSE,
                            header = T, row.names = 1, encoding = "UTF-8")
dim(unique.interaction.cluster) # [1] 1717   18

# save all TF names
TF_list <- unique.interaction.cluster$TFname
write.csv(TF_list, file = paste0(data_folder, "TF_list.csv"),
          na = 'NA', fileEncoding="UTF-8")
# get TF ID
TF_list <- read.csv(paste0(data_folder, "TF_list.csv"), 
                    stringsAsFactors = FALSE,
                    header = T, row.names = 1, encoding = "UTF-8")
rownames(TF_list) <- TF_list$TF_name

# ----
# 1 focus identical TFs in cluster 4 --------
identical_TF_list <- unique.interaction.cluster[unique.interaction.cluster$count >1 & unique.interaction.cluster$in.cluster4.ratio == 1, "TFname"]
identical_TF_ID <- TF_list[identical_TF_list,]
colnames(identical_TF_ID) <- c("identical_TF",
                               "identical_TF_ID")
rownames(identical_TF_ID) <- 1:length(identical_TF_list)

# generate n background TFlist with same length
num_bglist <- 50
background_df <- data.frame(matrix(ncol = 0, nrow = length(identical_TF_list)))
for (i in 1:num_bglist) {
  set.seed(10 * i)
  background_list <- sample(rownames(TF_list), length(identical_TF_list))
  background_df_i <- TF_list[background_list,]
  colnames(background_df_i) <- c(paste0("background_", i),
                                    paste0("background_", i, "_ID"))
  rownames(background_df_i) <- 1:length(identical_TF_list)
  background_df <- cbind(background_df, background_df_i)
}

TF_list_bg <- cbind(identical_TF_ID, background_df)
write.csv(TF_list_bg, file = paste0(data_folder, "TF_list_identical+background.csv"),
          na = 'NA', fileEncoding="UTF-8")

# ----
# 2 do GO enrichment --------
# change to entrezID
# CaseGeneSet <- bitr(identical_TF, 'SYMBOL', "ENTREZID", "org.Hs.eg.db")[, "ENTREZID"]
# often fail to mapping, use online result from metascape
TF_list_bg_ID <- read.csv(paste0(data_folder, "TF_list_identical+background.csv"), 
                          stringsAsFactors = FALSE,
                          header = T, row.names = 1, encoding = "UTF-8")
identical.go <- enrichGO(
  gene = TF_list_bg_ID$identical_TF_ID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = T
)

# dotplot(identical.go)


identical.go.res <- identical.go@result
identical.go.res$GeneRatio <- apply(identical.go.res, 1, function(row) {
  parts <- strsplit(row['GeneRatio'], "/")[[1]]
  return(as.numeric(parts[1]) / as.numeric(parts[2]))
})
write.csv(identical.go.res, file = paste0(data_folder, "GO_res_identical.csv"),
          na = 'NA', fileEncoding="UTF-8")
identical.go.res <- identical.go.res[identical.go.res$Count >= 3,]

focus_go_list <- rownames(identical.go.res)
length(focus_go_list) # [1] 23

# get enrichment result for focus go on background TFs
background.go_res <- data.frame(matrix(ncol = ncol(identical.go.res), nrow = 0))
colnames(background.go_res) <- colnames(identical.go.res)
for (i in 1:num_bglist) {
  name_bg_ID <- paste0("background_", i, "_ID")
  go <- enrichGO(
    gene = TF_list_bg_ID[,name_bg_ID],
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    readable = T
  )
  for (focus_go in focus_go_list) {
    if (focus_go %in% rownames(go@result)) {
      # print(focus_go)
      background.go_res <- rbind(background.go_res, go@result[focus_go,])
    }
  }
}

background.go_res$GeneRatio <- apply(background.go_res, 1, function(row) {
  parts <- strsplit(row['GeneRatio'], "/")[[1]]
  return(as.numeric(parts[1]) / as.numeric(parts[2]))
})

write.csv(background.go_res, file = paste0(data_folder, "GO_res_background.csv"),
          na = 'NA', fileEncoding="UTF-8")

p <- ggplot() +
  geom_boxplot(data = background.go_res, aes(x = GeneRatio, y = Description), color = "gray", alpha = 0.5) + 
  geom_point(data = identical.go.res, aes(x = GeneRatio, y = Description, size = Count, color = pvalue), alpha = 1) + 
  scale_size_continuous(range = c(3, 6)) +
  scale_color_gradient(low = "#d53e4f", high = "#3288bd") +
  labs(x = "GeneRatio", y = NULL, size = "Count", color = "p Value") +
  theme_test() +
  theme(text = element_text(size = 8),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 8),
        axis.title = element_text(size = 15),
        legend.position = "right",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        panel.grid.major.x = element_line(size = 0.5, linetype = "dashed", color = "gray"),
        panel.spacing = unit(0.1, "lines"))
p
ggsave(paste0(image_folder, "go_analysis_raw.pdf"), p, 
       width = 12, height = 6)

# ----
# 3 select significant terms --------
# save terms with gene ratio greater than Q3
# remove replicate gene terms
focus_go_list_filter <- c("GO:0051098", "GO:0031667", "GO:0006367", "GO:0048706", "GO:0098732", "GO:0030856", "GO:0048511")
identical.go.res_filter <- identical.go.res[focus_go_list_filter,]
identical.go.res_filter$Description <- factor(identical.go.res_filter$Description, levels = rev(identical.go.res_filter$Description))
background.go_res_filter <- background.go_res[background.go_res$ID %in% focus_go_list_filter,]

p <- ggplot() +
  geom_boxplot(data = background.go_res_filter, aes(x = GeneRatio, y = Description), color = "gray", alpha = 0.5) + 
  geom_point(data = identical.go.res_filter, aes(x = GeneRatio, y = Description, size = Count, color = -log10(pvalue)), alpha = 1) + 
  scale_size_continuous(range = c(6, 12)) +
  scale_color_gradient(low = "#3288bd", high = "#d53e4f") +
  labs(x = "Gene Ratio", y = NULL, size = "Count", color = "-Log10(P Value)",
       title = "GO enrichment for identical TFs") +
  theme_test() +
  theme(text = element_text(size = 8),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 18),
        axis.title = element_text(size = 15),
        legend.position = "right",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        panel.grid.major.x = element_line(size = 0.5, linetype = "dashed", color = "gray"),
        panel.border = element_rect(linewidth = 2),
        panel.spacing = unit(0.1, "lines"))
p
ggsave(paste0(image_folder, "go_analysis.png"), p, 
       width = 10, height = 5)
