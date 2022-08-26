
library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)

#microglia subcluster heatmap
mat <- read.csv("output/files/micro_hmp_mat.csv", row.names = 1) %>% as.matrix()
md <- readRDS("output/seurat_metadata_subset.rds")
rownames(md) <- gsub("-", ".", rownames(md))

homeo_clusters <- c("micro.3", "micro.4", "micro.25")
dam_clusters <- c("micro.0", "micro.1", "micro.5", "micro.6", "micro.7", "micro.9", "micro.12",
                  "micro.13", "micro.17", "micro.23", "micro.36")

md <- md[colnames(mat),] %>% 
  select(orig.ident, Treatment, Region, seurat_clusters, celltype)

md <- md %>% 
  mutate(DAM = case_when(celltype %in% homeo_clusters ~ "No",
                         celltype %in% dam_clusters ~ "Yes"))

#add enrichment scores
tmp <- read.csv("output/files/enrichment_microglial_clusters.csv", row.names = 1)
md <- md %>% left_join(tmp, by = "celltype")

# what's the value range in the matrix
mat <- (mat - rowMeans(mat))/matrixStats::rowSds(mat)
quantile(mat, c(0.1, 0.95))

col_fun = circlize::colorRamp2(c(-1, -0.7, -0.5, -0.3, -0.1, 0, 0.3, 0.7, 1.1, 1.5, 2), rev(brewer.pal(11, "PuOr")))

range(md[,c(7:9)])
col_cells = circlize::colorRamp2(c(-1, -0.7, -0.5, -0.3, -0.1, 0, 1, 3, 5, 6, 7), rev(brewer.pal(11, "PuOr")))

dark2 <- brewer.pal(8, "Dark2")

ha <- HeatmapAnnotation(Treatment = md$Treatment,
                        DAM = md$DAM,
                        Region = md$Region,
                        log2FC_hp_vs_cx = md$enrichment_hp_cx,
                        log2FC_hp_RML_vs_Mock = md$enrichment_hp_RML_Mock,
                        log2FC_cx_RML_vs_Mock = md$enrichment_cx_RML_Mock,
                        col = list(Treatment = c("Mock" = "navy", "RML" = "firebrick"),
                                   DAM = c("No" = dark2[8], "Yes" = dark2[7]),
                                   Region = c("hp" = dark2[6], "cx" = dark2[5]),
                                   log2FC_hp_vs_cx = col_cells,
                                   log2FC_hp_RML_vs_Mock = col_cells,
                                   log2FC_cx_RML_vs_Mock = col_cells))

mkrs <- read.delim("output/files/hmap_micro_mkrs.txt", header = FALSE)$V1
mkrs_plot <- rownames(mat) %in% mkrs

hmp <- Heatmap(mat,
               col = col_fun,
               name = "Expression",
               show_column_names = FALSE,
               column_split = md$seurat_clusters, 
               row_km = 8,
               row_dend_reorder = TRUE,
               column_dend_reorder = TRUE,
               top_annotation = ha,
               row_names_gp = gpar(fontsize = 8),
               use_raster = TRUE,
               raster_quality = 10) +
  rowAnnotation(foo = anno_mark(at=which(mkrs_plot), labels = rownames(mat)[mkrs_plot], labels_gp = gpar(fontsize=8)))
pdf("output/plots/heatmaps/micro_clean_hmap.pdf", width=8.5, height=11)
hmp <- draw(hmp, merge_legend = TRUE)
dev.off()

rownames(mat)[order(rownames(mat))] %>% writeClipboard()

rownames(mat[row_order(hmp)[[1]],]) %>% writeClipboard()
rownames(mat[row_order(hmp)[[2]],]) %>% writeClipboard()
rownames(mat[row_order(hmp)[[3]],]) %>% writeClipboard()
rownames(mat[row_order(hmp)[[4]],]) %>% writeClipboard()
rownames(mat[row_order(hmp)[[5]],]) %>% writeClipboard()
rownames(mat[row_order(hmp)[[6]],]) %>% writeClipboard()
rownames(mat[row_order(hmp)[[7]],]) %>% writeClipboard()
rownames(mat[row_order(hmp)[[8]],]) %>% writeClipboard()
#rownames(mat[row_order(hmp)[[9]],]) %>% writeClipboard()
#rownames(mat[row_order(hmp)[[10]],]) %>% writeClipboard()

write.table(rownames(mat[row_order(hmp)[[1]],]), "output/files/heatmap cluster genes/micro_cluster_1.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rownames(mat[row_order(hmp)[[2]],]), "output/files/heatmap cluster genes/micro_cluster_2.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rownames(mat[row_order(hmp)[[3]],]), "output/files/heatmap cluster genes/micro_cluster_3.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rownames(mat[row_order(hmp)[[4]],]), "output/files/heatmap cluster genes/micro_cluster_4.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rownames(mat[row_order(hmp)[[5]],]), "output/files/heatmap cluster genes/micro_cluster_5.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rownames(mat[row_order(hmp)[[6]],]), "output/files/heatmap cluster genes/micro_cluster_6.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rownames(mat[row_order(hmp)[[7]],]), "output/files/heatmap cluster genes/micro_cluster_7.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rownames(mat[row_order(hmp)[[8]],]), "output/files/heatmap cluster genes/micro_cluster_8.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(rownames(mat[row_order(hmp)[[9]],]), "output/files/heatmap cluster genes/micro_cluster_9.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(rownames(mat[row_order(hmp)[[10]],]), "output/files/heatmap cluster genes/micro_cluster_10.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


###Astrocyte heatmap
mat <- read.csv("output/files/astro_hmp_mat.csv", row.names = 1) %>% as.matrix()
md <- readRDS("output/seurat_metadata_astrocyte_subset.rds")
rownames(md) <- gsub("-", ".", rownames(md))

md <- md[colnames(mat),] %>% select(orig.ident, Treatment, Region, seurat_clusters, celltype)

tmp <- read.csv("output/files/enrichment_astrocyte_clusters.csv", row.names = 1)
md$celltype <- md$seurat_clusters
tmp$celltype <- as.factor(tmp$celltype)
md <- md %>% left_join(tmp, by = "celltype")

# what's the value range in the matrix
mat <- (mat - rowMeans(mat))/matrixStats::rowSds(mat)
quantile(mat, c(0.1, 0.95))

col_fun = circlize::colorRamp2(c(-1, -0.7, -0.5, -0.3, -0.1, 0, 0.3, 0.7, 1.1, 1.5, 2), rev(brewer.pal(11, "PuOr")))

range(md[,c(6:8)])
col_cells = circlize::colorRamp2(c(-5, -4, -3, -2, -2, 0, 1, 2, 3, 3.5, 4), rev(brewer.pal(11, "PuOr")))

dark2 <- brewer.pal(8, "Dark2")

ha <- HeatmapAnnotation(Treatment = md$Treatment,
                        Region = md$Region,
                        log2FC_hp_vs_cx = md$enrichment_hp_cx,
                        log2FC_hp_RML_vs_Mock = md$enrichment_hp_RML_Mock,
                        log2FC_cx_RML_vs_Mock = md$enrichment_cx_RML_Mock,
                        col = list(Treatment = c("Mock" = "navy", "RML" = "firebrick"),
                                   Region = c("hp" = dark2[6], "cx" = dark2[5]),
                                   log2FC_hp_vs_cx = col_cells,
                                   log2FC_hp_RML_vs_Mock = col_cells,
                                   log2FC_cx_RML_vs_Mock = col_cells))

mkrs <- read.delim("output/files/hmap_astro_mkrs.txt", header = FALSE)$V1
mkrs_plot <- rownames(mat) %in% mkrs

hmp <- Heatmap(mat,
               col = col_fun,
               name = "Expression",
               show_column_names = FALSE,
               column_split = md$seurat_clusters, 
               row_km = 7,
               row_dend_reorder = TRUE,
               column_dend_reorder = TRUE,
               top_annotation = ha,
               row_names_gp = gpar(fontsize = 8),
               use_raster = TRUE,
               raster_quality = 10) +
  rowAnnotation(foo = anno_mark(at=which(mkrs_plot), labels = rownames(mat)[mkrs_plot], labels_gp = gpar(fontsize=8)))

pdf("output/plots/heatmaps/astro_clean_hmap.pdf", width=8.5, height=11)
hmp <- draw(hmp, merge_legend = TRUE)
dev.off()

writeClipboard(rownames(mat)[order(rownames(mat))])

rownames(mat[row_order(hmp)[[1]],]) %>% writeClipboard()
rownames(mat[row_order(hmp)[[2]],]) %>% writeClipboard()
rownames(mat[row_order(hmp)[[3]],]) %>% writeClipboard()
rownames(mat[row_order(hmp)[[4]],]) %>% writeClipboard()
rownames(mat[row_order(hmp)[[5]],]) %>% writeClipboard()
rownames(mat[row_order(hmp)[[6]],]) %>% writeClipboard()
rownames(mat[row_order(hmp)[[7]],]) %>% writeClipboard()

write.table(rownames(mat[row_order(hmp)[[1]],]), "output/files/heatmap cluster genes/astro_cluster_1.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rownames(mat[row_order(hmp)[[2]],]), "output/files/heatmap cluster genes/astro_cluster_2.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rownames(mat[row_order(hmp)[[3]],]), "output/files/heatmap cluster genes/astro_cluster_3.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rownames(mat[row_order(hmp)[[4]],]), "output/files/heatmap cluster genes/astro_cluster_4.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rownames(mat[row_order(hmp)[[5]],]), "output/files/heatmap cluster genes/astro_cluster_5.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rownames(mat[row_order(hmp)[[6]],]), "output/files/heatmap cluster genes/astro_cluster_6.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rownames(mat[row_order(hmp)[[7]],]), "output/files/heatmap cluster genes/astro_cluster_7.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)



##Neuron subcluster heatmap
mat <- read.csv("output/files/neuron_hmp_mat.csv", row.names = 1) %>% as.matrix()
md <- readRDS("output/seurat_metadata_neuron_subsets_2.rds")
rownames(md) <- gsub("-", ".", rownames(md))

md <- md[colnames(mat),] %>% select(orig.ident, Treatment, Region, seurat_clusters, celltype)
md <- md %>% mutate(type = case_when(seurat_clusters %in% 1 ~ "npc",
                                     seurat_clusters %in% c(2,3,5,12,14) ~ "diff",
                                     seurat_clusters %in% 7 ~ "cr",
                                     seurat_clusters %in% c(0,11) ~ "mature",
                                     seurat_clusters %in% c(4,6,8,9,15) ~ "ex",
                                     seurat_clusters %in% c(10,13) ~ "inh")
                    )

#add enrichment scores
tmp <- read.csv("output/files/enrichment_neuron_clusters.csv", row.names = 1)
md <- md %>% left_join(tmp, by = "celltype")


# what's the value range in the matrix
mat <- (mat - rowMeans(mat))/matrixStats::rowSds(mat)
quantile(mat, c(0.1, 0.95))

col_fun = circlize::colorRamp2(c(-1, -0.7, -0.5, -0.3, -0.1, 0, 0.3, 0.7, 1.1, 1.5, 2), rev(brewer.pal(11, "PuOr")))

range(md[,c(7:9)])
col_cells = circlize::colorRamp2(c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5), rev(brewer.pal(11, "PuOr")))

tcls <- brewer.pal(8, "Dark2")

ha <- HeatmapAnnotation(Treatment = md$Treatment,
                        Type = md$type,
                        Region = md$Region,
                        log2FC_hp_vs_cx = md$enrichment_hp_cx,
                        log2FC_hp_RML_vs_Mock = md$enrichment_hp_RML_Mock,
                        log2FC_cx_RML_vs_Mock = md$enrichment_cx_RML_Mock,
                        col = list(Treatment = c("Mock" = "navy", "RML" = "firebrick"),
                                   Type = c("npc" = tcls[1], "diff" = tcls[2],
                                            "cr" = tcls[3], "mature" = tcls[4],
                                            "ex" = tcls[5], "inh" = tcls[6]),
                                   Region = c("hp" = tcls[7], "cx" = tcls[8]),
                                   log2FC_hp_vs_cx = col_cells,
                                   log2FC_hp_RML_vs_Mock = col_cells,
                                   log2FC_cx_RML_vs_Mock = col_cells))

mkrs <- read.delim("output/files/hmap_neuron_mkrs.txt", header = FALSE)$V1
mkrs_plot <- rownames(mat) %in% mkrs

hmp <- Heatmap(mat,
               col = col_fun,
               name = "Expression",
               show_column_names = FALSE,
               column_split = md$seurat_clusters, 
               row_km = 5,
               row_dend_reorder = TRUE,
               column_dend_reorder = TRUE,
               top_annotation = ha,
               row_names_gp = gpar(fontsize = 8),
               use_raster = TRUE,
               raster_quality = 10) +
  rowAnnotation(foo = anno_mark(at=which(mkrs_plot), labels = rownames(mat)[mkrs_plot], labels_gp = gpar(fontsize=8)))
pdf("output/plots/heatmaps/neuron_clean_hmap.pdf", width=8.5, height=11)
hmp <- draw(hmp, merge_legend = TRUE)
dev.off()

rownames(mat)[order(rownames(mat))] %>% writeClipboard()

rownames(mat[row_order(hmp)[[1]],]) %>% writeClipboard()
rownames(mat[row_order(hmp)[[2]],]) %>% writeClipboard()
rownames(mat[row_order(hmp)[[3]],]) %>% writeClipboard()
rownames(mat[row_order(hmp)[[4]],]) %>% writeClipboard()
rownames(mat[row_order(hmp)[[5]],]) %>% writeClipboard()

write.table(rownames(mat[row_order(hmp)[[1]],]), "output/files/heatmap cluster genes/neu_cluster_1.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rownames(mat[row_order(hmp)[[2]],]), "output/files/heatmap cluster genes/neu_cluster_2.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rownames(mat[row_order(hmp)[[3]],]), "output/files/heatmap cluster genes/neu_cluster_3.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rownames(mat[row_order(hmp)[[4]],]), "output/files/heatmap cluster genes/neu_cluster_4.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(rownames(mat[row_order(hmp)[[5]],]), "output/files/heatmap cluster genes/neu_cluster_5.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


rownames(mat[row_order(hmp)[[3]],]) %>% union(rownames(mat[row_order(hmp)[[4]],])) %>% writeClipboard()


# ##Mature neuron subcluster heatmap
# mat <- read.csv("output/files/neuron_9_15_vs_10_13_hmp_mat.csv", row.names = 1) %>% as.matrix()
# md <- readRDS("output/seurat_metadata_neuron_subsets_2.rds")
# rownames(md) <- gsub("-", ".", rownames(md))
# 
# md <- md[colnames(mat),] %>% select(orig.ident, Treatment, Region, seurat_clusters, celltype)
# 
# md <- md %>% 
#   mutate(RML_status = case_when(seurat_clusters %in% c("9", "15") ~ "sustained",
#                                 seurat_clusters %in% c("10", "13") ~ "depleted"))
# 
# 
# # what's the value range in the matrix
# mat <- (mat - rowMeans(mat))/matrixStats::rowSds(mat)
# quantile(mat, c(0.1, 0.95))
# 
# col_fun = circlize::colorRamp2(c(-1, -0.7, -0.5, -0.3, -0.1, 0, 0.3, 0.7, 1.1, 1.5, 2), rev(brewer.pal(11, "PuOr")))
# 
# ha <- HeatmapAnnotation(Treatment = md$Treatment,
#                         RML_status = md$RML_status,
#                         col = list(Treatment = c("Mock" = "navy", "RML" = "firebrick"),
#                                    RML_status = c("sustained" = "lightgrey", "depleted" = "black")))
# 
# mkrs <- read.delim("output/files/hmap__mature_neuron_mkrs.txt", header = FALSE)$V1
# mkrs_plot <- rownames(mat) %in% mkrs
# 
# hmp <- Heatmap(mat,
#                col = col_fun,
#                name = "Expression",
#                show_column_names = FALSE,
#                column_split = md$seurat_clusters, 
#                row_km = 3,
#                top_annotation = ha,
#                row_names_gp = gpar(fontsize = 8)) +
#   rowAnnotation(foo = anno_mark(at=which(mkrs_plot), labels = rownames(mat)[mkrs_plot], labels_gp = gpar(fontsize=8)))
# pdf("output/plots/heatmaps/mature_neuron_clean_hmap.pdf", width=8.5, height=11)
# hmp <- draw(hmp, merge_legend = TRUE)
# dev.off()
# 
# rownames(mat)[order(rownames(mat))] %>% writeClipboard()
# 
# rownames(mat[row_order(hmp)[[1]],]) %>% writeClipboard()
# rownames(mat[row_order(hmp)[[2]],]) %>% writeClipboard()
# rownames(mat[row_order(hmp)[[3]],]) %>% writeClipboard()
# 
# write.table(rownames(mat[row_order(hmp)[[1]],]), "output/files/heatmap cluster genes/mature_neu_cluster_1.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
# write.table(rownames(mat[row_order(hmp)[[2]],]), "output/files/heatmap cluster genes/mature_neu_cluster_2.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
# write.table(rownames(mat[row_order(hmp)[[3]],]), "output/files/heatmap cluster genes/mature_neu_cluster_3.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
