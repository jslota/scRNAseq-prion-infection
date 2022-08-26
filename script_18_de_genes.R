###n DE genes in each cluster
#2022-06-01
#Jessy Slota

library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)

fls <- Sys.glob("output/files/DE res/response*")
clusters <- rep(gsub(".csv", "", gsub("output/files/DE res/response_", "", fls)),2)
clusters <- clusters[order(clusters)]
out <- data.frame(cluster=clusters,
                  dir=rep(c("up", "dn"), length(fls)),
                  n=0)
de <- list()
for (i in Sys.glob("output/files/DE res/response*")) {
  print(i)
  x <- gsub(".csv", "", gsub("output/files/DE res/response_", "", i))
  
  #Reads DE results
  tmp <- read.csv(i)
  #filter out non-protein coding and technical genes
  nonsense <- c(tmp$X[grep("^RP", tmp$X)], 
                tmp$X[grep("^Gm[0-9]", tmp$X)],
                tmp$X[grep("^AC[0-9]", tmp$X)],
                tmp$X[grep("[0-9]Rik", tmp$X)],
                tmp$X[grep("AY[0-9]", tmp$X)],
                tmp$X[grep("Malat1", tmp$X)],
                tmp$X[grep("^mt-", tmp$X)],
                tmp$X[grep("^Rp[sl]", tmp$X)],
                tmp$X[grep("Hb[^(p)]", tmp$X)],
                tmp$X[grep("^Ndu", tmp$X)],
                tmp$X[grep("^Cox", tmp$X)])

  out[out$cluster==x & out$dir=="up",]$n <- tmp  %>% 
    filter(p_val_adj < 0.05, avg_log2FC > 0.5, pct.1 > 0.5, pct.2 > 0.25, (pct.1-pct.2) > -0.1, !X %in% nonsense) %>%
    nrow()
  
  tmp.1 <- tmp  %>% 
    filter(p_val_adj < 0.05, avg_log2FC > 0.5, pct.1 > 0.5, pct.2 > 0.25, (pct.1-pct.2) > -0.1, !X %in% nonsense)
  
  if(nrow(tmp.1) > 0) {
    de[[paste0(x,"1")]] <- tmp.1
    de[[paste0(x,"1")]]$cluster <- x
  }

  out[out$cluster==x & out$dir=="dn",]$n <- tmp  %>% 
    filter(p_val_adj < 0.05, avg_log2FC < -0.5, pct.1 > 0.25, pct.2 > 0.5, (pct.1-pct.2) < 0.1, !X %in% nonsense) %>%
    nrow()
  
  tmp.1 <- tmp  %>% 
    filter(p_val_adj < 0.05, avg_log2FC < -0.5, pct.1 > 0.25, pct.2 > 0.5, (pct.1-pct.2) < 0.1, !X %in% nonsense)
  
  if(nrow(tmp.1) > 0) {
    de[[paste0(x,"2")]] <- tmp.1
    de[[paste0(x,"2")]]$cluster <- x
  }

}
de <- do.call(rbind, de)

ggplot(out, aes(x=cluster, y=n, fill=dir, label=n)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(position = position_dodge2(width=0.9)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(11, "PuOr")[c(11,2)]) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90))
ggsave("output/plots/n_DE_genes.pdf", width = 8, height = 4, units = "in")


plt_dt <- de %>% reshape2::dcast(formula = X~cluster, value.var = "avg_log2FC", fill = 0)
rownames(plt_dt) <- plt_dt$X
plt_dt <- plt_dt[,-1]
plt_dt <- as.matrix(plt_dt)

col_fun = circlize::colorRamp2(c(-10, -7, -5, -3, -1, 0, 1, 3, 5, 7, 10), rev(brewer.pal(11, "PuOr")))

#rml <- readxl::read_excel("output/files/idmap.xlsx") %>% pull(symbol)
# 
#valid <- rownames(plt_dt) %>% intersect(rml)
#gene_types %>% filter(cat == "micro") %>% rownames() %>% intersect(valid) %>% writeClipboard()


mkrs <- read.delim("output/files/hmap_de_gene_mkrs.txt", header = FALSE)$V1
mkrs_plot <- rownames(plt_dt) %in% mkrs

colnames(plt_dt)
cell_types <- data.frame(row.names = colnames(plt_dt),
                         cell_type = c("astro", "astro","endo", "endo", "endo", "endo", "other", "g.neu", "g.neu",
                                       "g.neu", "g.neu", "im.neu", "im.neu", "im.neu", "other", "micro", "micro", "micro",
                                       "micro", "micro", "micro", "micro", "micro", "micro",  "micro", "micro", "micro",
                                       "opc", "endo", "pvm", "pvm", "endo", "endo"))

cols <- viridisLite::viridis(8)
ha <- HeatmapAnnotation(cell_type = cell_types$cell_type,
                        col = list(cell_type = c("astro" = cols[1], "endo" = cols[2], "other" = cols[3],
                                                  "g.neu" = cols[4], "im.neu" = cols[4], "micro" = cols[6],
                                                  "opc" = cols[7], "pvm" = cols[8])))


gene_types <- data.frame(row.names = rownames(plt_dt),
                         cat = rep("other", nrow(plt_dt)))

tmp <- de %>% filter(cluster %in% c("pvm.22", "pvm.29")) %>% pull(X) %>% unique()
gene_types[tmp,] = "pvm"

tmp <- de %>% filter(cluster %in% c("opc.19")) %>% pull(X) %>% unique()
gene_types[tmp,] = "opc"

tmp <- de %>% filter(cluster %in% c("endo.11", "endo.2", "endo.27", "endo.38", "peri.14", "smc.21", "vlmc.30")) %>% pull(X) %>% unique()
gene_types[tmp,] = "endo"

tmp <- de %>% filter(cluster %in% c("im.neu.18", "im.neu.24", "im.neu.32")) %>% pull(X) %>% unique()
gene_types[tmp,] = "im.neu"

tmp <- de %>% filter(cluster %in% c("g.neu.16", "g.neu.31", "g.neu.33", "g.neu.34")) %>% pull(X) %>% unique()
gene_types[tmp,] = "g.neu"

tmp <- de %>% filter(cluster %in% c("micro.0", "micro.1", "micro.12", "micro.13", "micro.17", "micro.25", "micro.3", "micro.4", "micro.5", "micro.6", "micro.7", "micro.7", "micro.9")) %>% pull(X) %>% unique()
gene_types[tmp,] = "micro"

tmp <- de %>% filter(cluster %in% c("astro.10", "astro.20")) %>% pull(X) %>% unique()
gene_types[tmp,] = "astro"

hmp <- Heatmap(plt_dt,
               col = col_fun,
               name = "log2(FC)",
               column_split = cell_types$cell_type,
               row_split = gene_types$cat,
               column_dend_reorder = TRUE,
               row_dend_reorder = TRUE,
               top_annotation = ha) +
  rowAnnotation(foo = anno_mark(at=which(mkrs_plot),
                                labels = rownames(plt_dt)[mkrs_plot],
                                labels_gp = gpar(fontsize=8)))
pdf("output/plots/heatmaps/DE_gene_hmap.pdf", width=8.5, height=11)
hmp <- draw(hmp, merge_legend = TRUE)
dev.off()

de %>%
  filter(cluster == "g.neu.16", avg_log2FC > 0) %>%
  arrange(desc(avg_log2FC)) %>%
  pull(X) %>%
  writeClipboard()

de %>%
  filter(cluster == "g.neu.16", avg_log2FC < 0) %>%
  arrange(avg_log2FC) %>%
  pull(X) %>%
  writeClipboard()

gene_types %>% filter(cat == "micro") %>% rownames() %>% write.table("output/files/heatmap cluster genes/DE_micro.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
gene_types %>% filter(cat == "astro") %>% rownames() %>% write.table("output/files/heatmap cluster genes/DE_astro.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
gene_types %>% filter(cat == "pvm") %>% rownames() %>% write.table("output/files/heatmap cluster genes/DE_pvm.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
gene_types %>% filter(cat == "opc") %>% rownames() %>% write.table("output/files/heatmap cluster genes/DE_opc.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
gene_types %>% filter(cat == "endo") %>% rownames() %>% write.table("output/files/heatmap cluster genes/DE_endo.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
gene_types %>% filter(cat == "g.neu") %>% rownames() %>% write.table("output/files/heatmap cluster genes/DE_g.neu.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
gene_types %>% filter(cat == "im.neu") %>% rownames() %>% write.table("output/files/heatmap cluster genes/DE_im.neu.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


################################################################################
###########################   neuron subclusters   #############################
################################################################################



fls <- Sys.glob("output/files/Neu DE res/response*")
clusters <- rep(gsub(".csv", "", gsub("output/files/Neu DE res/response_", "", fls)),2)
clusters <- clusters[order(clusters)]
out <- data.frame(cluster=clusters,
                  dir=rep(c("up", "dn"), length(fls)),
                  n=0)
de <- list()
for (i in fls) {
  print(i)
  x <- gsub(".csv", "", gsub("output/files/Neu DE res/response_", "", i))
  
  #Reads DE results
  tmp <- read.csv(i)
  #filter out non-protein coding and technical genes
  nonsense <- c(tmp$X[grep("^RP", tmp$X)], 
                tmp$X[grep("^Gm[0-9]", tmp$X)],
                tmp$X[grep("^AC[0-9]", tmp$X)],
                tmp$X[grep("[0-9]Rik", tmp$X)],
                tmp$X[grep("AY[0-9]", tmp$X)],
                tmp$X[grep("Malat1", tmp$X)],
                tmp$X[grep("^mt-", tmp$X)],
                tmp$X[grep("^Rp[sl]", tmp$X)],
                tmp$X[grep("Hb[^(p)]", tmp$X)],
                tmp$X[grep("^Ndu", tmp$X)],
                tmp$X[grep("^Cox", tmp$X)])
  
  out[out$cluster==x & out$dir=="up",]$n <- tmp  %>% 
    filter(p_val_adj < 0.05, avg_log2FC > 0.5, pct.1 > 0.5, pct.2 > 0.25, (pct.1-pct.2) > -0.1, !X %in% nonsense) %>%
    nrow()
  
  tmp.1 <- tmp  %>% 
    filter(p_val_adj < 0.05, avg_log2FC > 0.5, pct.1 > 0.5, pct.2 > 0.25, (pct.1-pct.2) > -0.1, !X %in% nonsense)
  
  if(nrow(tmp.1) > 0) {
    de[[paste0(x,"1")]] <- tmp.1
    de[[paste0(x,"1")]]$cluster <- x
  }
  
  out[out$cluster==x & out$dir=="dn",]$n <- tmp  %>% 
    filter(p_val_adj < 0.05, avg_log2FC < -0.5, pct.1 > 0.25, pct.2 > 0.5, (pct.1-pct.2) < 0.1, !X %in% nonsense) %>%
    nrow()
  
  tmp.1 <- tmp  %>% 
    filter(p_val_adj < 0.05, avg_log2FC < -0.5, pct.1 > 0.25, pct.2 > 0.5, (pct.1-pct.2) < 0.1, !X %in% nonsense)
  
  if(nrow(tmp.1) > 0) {
    de[[paste0(x,"2")]] <- tmp.1
    de[[paste0(x,"2")]]$cluster <- x
  }
  
}
de <- do.call(rbind, de)

ggplot(out, aes(x=cluster, y=n, fill=dir, label=n)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(position = position_dodge2(width=0.9)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(11, "PuOr")[c(11,2)]) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90))
ggsave("output/plots/n_DE_genes_neuron.pdf", width = 8, height = 4, units = "in")

#Heatmap
plt_dt <- de %>% reshape2::dcast(formula = X~cluster, value.var = "avg_log2FC", fill = 0)
rownames(plt_dt) <- plt_dt$X
plt_dt <- plt_dt[,-1]
plt_dt <- as.matrix(plt_dt)

col_fun = circlize::colorRamp2(c(-10, -7, -5, -3, -1, 0, 1, 3, 5, 7, 10), rev(brewer.pal(11, "PuOr")))

# rml <- readxl::read_excel("output/files/idmap.xlsx") %>% pull(symbol)
# # 
# valid <- rownames(plt_dt) %>% intersect(rml)
# gene_types %>% filter(cat == "micro") %>% rownames() %>% intersect(valid) %>% writeClipboard()


mkrs <- read.delim("output/files/hmap_de_gene_neuron_mkrs.txt", header = FALSE)$V1
mkrs_plot <- rownames(plt_dt) %in% mkrs
 
colnames(plt_dt)
cell_types <- data.frame(row.names = colnames(plt_dt),
                         cell_type = c("cr.neu", "diff.neu", "diff.neu", "diff.neu", "diff.neu", "diff.neu", "ex.neu", "ex.neu",
                                       "ex.neu", "ex.neu", "ex.neu", "inh.neu", "inh.neu", "m.neu", "m.neu", "npc"))


unique(cell_types$cell_type)
cols <- viridisLite::viridis(6)
ha <- HeatmapAnnotation(cell_type = cell_types$cell_type,
                        col = list(cell_type = c("cr.neu" = cols[1], "diff.neu" = cols[2], "ex.neu" = cols[3],
                                                 "inh.neu" = cols[4], "m.neu" = cols[4], "npc" = cols[6])))


gene_types <- data.frame(row.names = rownames(plt_dt),
                         cat = rep("other", nrow(plt_dt)))

tmp <- de %>% filter(cluster %in% c("m.neu.0", "m.neu.11")) %>% pull(X) %>% unique()
gene_types[tmp,] = "m.neu"

tmp <- de %>% filter(cluster %in% c("diff.neu.12", "diff.neu.14", "diff.neu.2", "diff.neu.3", "diff.neu.5")) %>% pull(X) %>% unique()
gene_types[tmp,] = "diff.neu"

tmp <- de %>% filter(cluster %in% c("inh.neu.10", "inh.neu.13")) %>% pull(X) %>% unique()
gene_types[tmp,] = "inh.neu"

tmp <- de %>% filter(cluster %in% c("ex.neu.15", "ex.neu.4", "ex.neu.6", "ex.neu.8", "ex.neu.9")) %>% pull(X) %>% unique()
gene_types[tmp,] = "ex.neu"

tmp <- de %>% filter(cluster %in% c("npc.1")) %>% pull(X) %>% unique()
gene_types[tmp,] = "npc"

tmp <- de %>% filter(cluster %in% c("cr.neu.7")) %>% pull(X) %>% unique()
gene_types[tmp,] = "cr.neu"

hmp <- Heatmap(plt_dt,
               col = col_fun,
               name = "log2(FC)",
               column_split = cell_types$cell_type,
               row_split = gene_types$cat,
               column_dend_reorder = TRUE,
               row_dend_reorder = TRUE,
               top_annotation = ha) +
  rowAnnotation(foo = anno_mark(at=which(mkrs_plot),
                                labels = rownames(plt_dt)[mkrs_plot],
                                labels_gp = gpar(fontsize=8)))
pdf("output/plots/heatmaps/Neuron_DE_gene_hmap.pdf", width=8.5, height=11)
hmp <- draw(hmp, merge_legend = TRUE)
dev.off()

gene_types %>% filter(cat == "m.neu") %>% rownames() %>% write.table("output/files/heatmap cluster genes/neu_DE_m.neu.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
gene_types %>% filter(cat == "diff.neu") %>% rownames() %>% write.table("output/files/heatmap cluster genes/neu_DE_diff.neu.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
gene_types %>% filter(cat == "inh.neu") %>% rownames() %>% write.table("output/files/heatmap cluster genes/neu_DE_inh.neu.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
gene_types %>% filter(cat == "ex.neu") %>% rownames() %>% write.table("output/files/heatmap cluster genes/neu_DE_ex.neu.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
gene_types %>% filter(cat == "npc") %>% rownames() %>% write.table("output/files/heatmap cluster genes/neu_DE_npc.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
gene_types %>% filter(cat == "cr.neu") %>% rownames() %>% write.table("output/files/heatmap cluster genes/neu_DE_cr.neu.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


de %>%
  filter(cluster == "m.neu.0", avg_log2FC > 0) %>%
  arrange(desc(avg_log2FC)) %>%
  pull(X) %>%
  writeClipboard()

de %>%
  filter(cluster == "npc.1", avg_log2FC < 0) %>%
  arrange(avg_log2FC) %>%
  pull(X) %>%
  writeClipboard()


de %>% filter(avg_log2FC > 0) %>% count(X) %>% arrange(desc(n)) %>% slice(1:50) %>% pull(X) %>% writeClipboard()
de %>% filter(avg_log2FC < 0) %>% count(X) %>% arrange(desc(n)) %>% slice(1:50) %>% pull(X) %>% writeClipboard()

de %>% filter(X == "Fgf14")

common <- union(de %>% filter(avg_log2FC > 0) %>% count(X) %>% arrange(desc(n)) %>% filter(n>2) %>% pull(X),
                de %>% filter(avg_log2FC < 0) %>% count(X) %>% arrange(desc(n)) %>% filter(n>2) %>% pull(X))

hmp <- Heatmap(plt_dt[common,],
               col = col_fun,
               name = "log2(FC)",
               column_split = cell_types$cell_type,
               column_dend_reorder = TRUE,
               row_dend_reorder = TRUE,
               top_annotation = ha)
pdf("output/plots/heatmaps/Neuron_common_DE_gene_hmap.pdf", width=8.5, height=11)
hmp <- draw(hmp, merge_legend = TRUE)
dev.off()

###################################################################################
#Get DE res
tmp <- list()
for (i in Sys.glob("output/files/DE res/response*")) {
  print(i)
  x <- gsub(".csv", "", gsub("output/files/DE res/response_", "", i))
  
  y <- read.csv(i) %>% filter(X %in% mkrs, p_val_adj < 0.05, abs(avg_log2FC) > 0.5, pct.1 > 0.25, pct.2 > 0.25)
  
  if(nrow(y) > 1) {
    tmp[[x]] <- y
    tmp[[x]]$cluster <- x
  }
  
}
tmp <- do.call(rbind, tmp)

plt_dt <- tmp %>% reshape2::dcast(formula = X~cluster, value.var = "avg_log2FC", fill = 0)
rownames(plt_dt) <- plt_dt$X
plt_dt <- plt_dt[,-1]
plt_dt <- as.matrix(plt_dt)


col_fun = circlize::colorRamp2(c(-10, -7, -5, -3, -1, 0, 1, 3, 5, 7, 10), rev(brewer.pal(11, "PuOr")))

rml <- readxl::read_excel("output/files/idmap.xlsx") %>% pull(symbol)

valid <- rownames(plt_dt) %>% intersect(rml)
plt_dt <- plt_dt[valid,]
mkrs <- read.delim("output/files/hmap_de_gene_mkrs.txt", header = FALSE)$V1
mkrs_plot <- rownames(plt_dt) %in% mkrs

hmp <- Heatmap(plt_dt,
               col = col_fun,
               name = "log2(FC)") +
  rowAnnotation(foo = anno_mark(at=which(mkrs_plot),
                                labels = rownames(plt_dt)[mkrs_plot],
                                labels_gp = gpar(fontsize=8)))
pdf("output/plots/heatmaps/DE_gene_hmap.pdf", width=8.5, height=11)
hmp <- draw(hmp, merge_legend = TRUE)
dev.off()



####Comparing between lists
tmp_fls <- fls
dt_up <- matrix(nrow = length(fls), ncol = length(fls), dimnames = list(fls, fls))
dt_dn <- matrix(nrow = length(fls), ncol = length(fls), dimnames = list(fls, fls))

for (i in fls) {
  for (j in tmp_fls) {
    print(paste0(i, " vs ", j))
    #Reads DE results
    tmp <- read.csv(i)
    tmpi <- tmp
    #filter out non-protein coding and technical genes
    nonsensei <- c(tmp$X[grep("^RP", tmp$X)], 
                   tmp$X[grep("^Gm[0-9]", tmp$X)],
                   tmp$X[grep("^AC[0-9]", tmp$X)],
                   tmp$X[grep("[0-9]Rik", tmp$X)],
                   tmp$X[grep("AY[0-9]", tmp$X)],
                   tmp$X[grep("Malat1", tmp$X)],
                   tmp$X[grep("^mt-", tmp$X)],
                   tmp$X[grep("^Rp[sl]", tmp$X)],
                   tmp$X[grep("Hb[^(p)]", tmp$X)],
                   tmp$X[grep("^Ndu", tmp$X)],
                   tmp$X[grep("^Cox", tmp$X)])
    
    tmp <- read.csv(j)
    tmpj <- tmp
    #filter out non-protein coding and technical genes
    nonsensej <- c(tmp$X[grep("^RP", tmp$X)], 
                   tmp$X[grep("^Gm[0-9]", tmp$X)],
                   tmp$X[grep("^AC[0-9]", tmp$X)],
                   tmp$X[grep("[0-9]Rik", tmp$X)],
                   tmp$X[grep("AY[0-9]", tmp$X)],
                   tmp$X[grep("Malat1", tmp$X)],
                   tmp$X[grep("^mt-", tmp$X)],
                   tmp$X[grep("^Rp[sl]", tmp$X)],
                   tmp$X[grep("Hb[^(p)]", tmp$X)],
                   tmp$X[grep("^Ndu", tmp$X)],
                   tmp$X[grep("^Cox", tmp$X)])
    rm(tmp)
    
    a <- tmpi %>% filter(p_val_adj < 0.05, avg_log2FC > 0.5, pct.1 > 0.5, pct.2 > 0.25, (pct.1-pct.2) > 0.1, !X %in% nonsensei) %>% pull(X)
    b <- tmpj %>% filter(p_val_adj < 0.05, avg_log2FC > 0.5, pct.1 > 0.5, pct.2 > 0.25, (pct.1-pct.2) > 0.1, !X %in% nonsensej) %>% pull(X)
    dt_up[i,j] <- length(intersect(a, b))/length(union(a, b))
    
    a <- tmpi %>% filter(p_val_adj < 0.05, avg_log2FC < -0.5, pct.1 > 0.25, pct.2 > 0.5, (pct.1-pct.2) < -0.1, !X %in% nonsensei) %>% pull(X)
    b <- tmpj %>% filter(p_val_adj < 0.05, avg_log2FC < -0.5, pct.1 > 0.25, pct.2 > 0.5, (pct.1-pct.2) < -0.1, !X %in% nonsensej) %>% pull(X)
    dt_dn[i,j] <- length(intersect(a, b))/length(union(a, b))
  }
  tmp_fls <- tmp_fls[-grep(i, tmp_fls)]
}

colnames(dt_up) <- gsub(".csv", "", gsub("output/files/DE res/response_", "", colnames(dt_up)))
rownames(dt_up) <- gsub(".csv", "", gsub("output/files/DE res/response_", "", rownames(dt_up)))
colnames(dt_dn) <- gsub(".csv", "", gsub("output/files/DE res/response_", "", colnames(dt_dn)))
rownames(dt_dn) <- gsub(".csv", "", gsub("output/files/DE res/response_", "", rownames(dt_dn)))

tmp <- dt_up %>% 
  reshape2::melt() %>%
  na.omit()

plt_cls <- colorRampPalette(RColorBrewer::brewer.pal(9, "Oranges"))(100)

ggplot(tmp, aes(x=Var1, y=Var2, fill=value, label=round(value,2))) +
  geom_tile() +
  geom_text(size=2) +
  scale_fill_gradientn(colours = plt_cls) +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0), 
        axis.title = element_blank())

tmp <- dt_dn %>% 
  reshape2::melt() %>%
  na.omit()

plt_cls <- colorRampPalette(RColorBrewer::brewer.pal(9, "Purples"))(100)

ggplot(tmp, aes(x=Var1, y=Var2, fill=value, label=round(value,2))) +
  geom_tile() +
  geom_text(size=2) +
  scale_fill_gradientn(colours = plt_cls) +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0), 
        axis.title = element_blank())





library(dplyr)
library(ggplot2)

fls <- Sys.glob("output/files/Neu DE res/response*")
clusters <- rep(gsub(".csv", "", gsub("output/files/Neu DE res/response_", "", fls)),2)
clusters <- clusters[order(clusters)]
out <- data.frame(cluster=clusters,
                  dir=rep(c("up", "dn"), length(fls)),
                  n=0)

for (i in Sys.glob("output/files/Neu DE res/response*")) {
  print(i)
  x <- gsub(".csv", "", gsub("output/files/Neu DE res/response_", "", i))
  
  #Reads DE results
  tmp <- read.csv(i)
  #filter out non-protein coding and technical genes
  nonsense <- c(tmp$X[grep("^RP", tmp$X)], 
                tmp$X[grep("^Gm[0-9]", tmp$X)],
                tmp$X[grep("^AC[0-9]", tmp$X)],
                tmp$X[grep("[0-9]Rik", tmp$X)],
                tmp$X[grep("AY[0-9]", tmp$X)],
                tmp$X[grep("Malat1", tmp$X)],
                tmp$X[grep("^mt-", tmp$X)],
                tmp$X[grep("^Rp[sl]", tmp$X)],
                tmp$X[grep("Hb[^(p)]", tmp$X)],
                tmp$X[grep("^Ndu", tmp$X)],
                tmp$X[grep("^Cox", tmp$X)])
  
  out[out$cluster==x & out$dir=="up",]$n <- tmp  %>% 
    filter(p_val_adj < 0.05, avg_log2FC > 0.5, pct.1 > 0.5, pct.2 > 0.25, (pct.1-pct.2) > 0.1, !X %in% nonsense) %>%
    nrow()
  
  out[out$cluster==x & out$dir=="dn",]$n <- tmp  %>% 
    filter(p_val_adj < 0.05, avg_log2FC < -0.5, pct.1 > 0.25, pct.2 > 0.5, (pct.1-pct.2) < -0.1, !X %in% nonsense) %>%
    nrow()
}

ggplot(out, aes(x=cluster, y=n, fill=dir, label=n)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text() +
  scale_fill_manual(values = RColorBrewer::brewer.pal(11, "PuOr")[c(11,2)]) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90))


tmp_fls <- fls
dt_up <- matrix(nrow = length(fls), ncol = length(fls), dimnames = list(fls, fls))
dt_dn <- matrix(nrow = length(fls), ncol = length(fls), dimnames = list(fls, fls))

for (i in fls) {
  for (j in tmp_fls) {
    print(paste0(i, " vs ", j))
    #Reads DE results
    tmp <- read.csv(i)
    tmpi <- tmp
    #filter out non-protein coding and technical genes
    nonsensei <- c(tmp$X[grep("^RP", tmp$X)], 
                   tmp$X[grep("^Gm[0-9]", tmp$X)],
                   tmp$X[grep("^AC[0-9]", tmp$X)],
                   tmp$X[grep("[0-9]Rik", tmp$X)],
                   tmp$X[grep("AY[0-9]", tmp$X)],
                   tmp$X[grep("Malat1", tmp$X)],
                   tmp$X[grep("^mt-", tmp$X)],
                   tmp$X[grep("^Rp[sl]", tmp$X)],
                   tmp$X[grep("Hb[^(p)]", tmp$X)],
                   tmp$X[grep("^Ndu", tmp$X)],
                   tmp$X[grep("^Cox", tmp$X)])
    
    tmp <- read.csv(j)
    tmpj <- tmp
    #filter out non-protein coding and technical genes
    nonsensej <- c(tmp$X[grep("^RP", tmp$X)], 
                   tmp$X[grep("^Gm[0-9]", tmp$X)],
                   tmp$X[grep("^AC[0-9]", tmp$X)],
                   tmp$X[grep("[0-9]Rik", tmp$X)],
                   tmp$X[grep("AY[0-9]", tmp$X)],
                   tmp$X[grep("Malat1", tmp$X)],
                   tmp$X[grep("^mt-", tmp$X)],
                   tmp$X[grep("^Rp[sl]", tmp$X)],
                   tmp$X[grep("Hb[^(p)]", tmp$X)],
                   tmp$X[grep("^Ndu", tmp$X)],
                   tmp$X[grep("^Cox", tmp$X)])
    rm(tmp)
    
    a <- tmpi %>% filter(p_val_adj < 0.05, avg_log2FC > 0.5, pct.1 > 0.5, pct.2 > 0.25, (pct.1-pct.2) > 0.1, !X %in% nonsensei) %>% pull(X)
    b <- tmpj %>% filter(p_val_adj < 0.05, avg_log2FC > 0.5, pct.1 > 0.5, pct.2 > 0.25, (pct.1-pct.2) > 0.1, !X %in% nonsensej) %>% pull(X)
    dt_up[i,j] <- length(intersect(a, b))/length(union(a, b))
    
    a <- tmpi %>% filter(p_val_adj < 0.05, avg_log2FC < -0.5, pct.1 > 0.25, pct.2 > 0.5, (pct.1-pct.2) < -0.1, !X %in% nonsensei) %>% pull(X)
    b <- tmpj %>% filter(p_val_adj < 0.05, avg_log2FC < -0.5, pct.1 > 0.25, pct.2 > 0.5, (pct.1-pct.2) < -0.1, !X %in% nonsensej) %>% pull(X)
    dt_dn[i,j] <- length(intersect(a, b))/length(union(a, b))
  }
  tmp_fls <- tmp_fls[-grep(i, tmp_fls)]
}

colnames(dt_up) <- gsub(".csv", "", gsub("output/files/Neu DE res/response_", "", colnames(dt_up)))
rownames(dt_up) <- gsub(".csv", "", gsub("output/files/Neu DE res/response_", "", rownames(dt_up)))
colnames(dt_dn) <- gsub(".csv", "", gsub("output/files/Neu DE res/response_", "", colnames(dt_dn)))
rownames(dt_dn) <- gsub(".csv", "", gsub("output/files/Neu DE res/response_", "", rownames(dt_dn)))

tmp <- dt_up %>% 
  reshape2::melt() %>%
  na.omit()

plt_cls <- colorRampPalette(RColorBrewer::brewer.pal(9, "Oranges"))(100)

ggplot(tmp, aes(x=Var1, y=Var2, fill=value, label=round(value,2))) +
  geom_tile() +
  geom_text(size=2) +
  scale_fill_gradientn(colours = plt_cls) +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0), 
        axis.title = element_blank())

tmp <- dt_dn %>% 
  reshape2::melt() %>%
  na.omit()

plt_cls <- colorRampPalette(RColorBrewer::brewer.pal(9, "Purples"))(100)

ggplot(tmp, aes(x=Var1, y=Var2, fill=value, label=round(value,2))) +
  geom_tile() +
  geom_text(size=2) +
  scale_fill_gradientn(colours = plt_cls) +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0), 
        axis.title = element_blank())