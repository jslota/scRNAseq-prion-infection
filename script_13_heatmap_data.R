###analysis of microglia cell clusters
#2022-06-01
#Jessy Slota

require(Seurat)
require(dplyr)

timestamp()

parallel::mclapply

dt.full <- readRDS("output/seurat_final_filtered.rds")
print("dataset is loaded...")
print(gc())
timestamp()

print(levels(Idents(dt.full)))
Idents(dt.full) <- "celltype"
print(levels(Idents(dt.full)))

###Microglia analysis
dt.sct <- subset(dt.full,
                 idents = c("micro.0", "micro.1", "micro.3", "micro.4", "micro.5", "micro.6",
                            "micro.7", "micro.9", "micro.12", "micro.13", "micro.17", "micro.23",
                            "micro.25", "micro.36")
)

print("dataset is subset...")
print(gc())
timestamp()

#microglial RML markers
mkrs <- read.delim("output/files/common_micro_mrks.txt", header = FALSE)$V1
for (i in Sys.glob("output/files/microglia subcluster markers/*")) {
  print(i)
  tmp <- read.csv(i)

  nonsense <- c(tmp$X[grep("^RP", tmp$X)],
                tmp$X[grep("^Gm[0-9]", tmp$X)],
                tmp$X[grep("^AC[0-9]", tmp$X)],
                tmp$X[grep("[0-9]Rik", tmp$X)])

  mkrs <- union(mkrs, tmp %>%
                  filter(p_val_adj < 1e-05, avg_log2FC > 2, pct.1 > 0.5, (pct.1-pct.2) > 0.2, !X %in% nonsense) %>%
                  arrange(desc(avg_log2FC)) %>%
                  slice(1:25) %>%
                  pull(X))
}

#Heatmap
p <- DoHeatmap(subset(dt.sct, downsample = 100), features = mkrs, size = 3, group.by = "Treatment")
pdf("output/plots/heatmaps/microglial_hmp.pdf", width = 11, height = 8.5)
p
dev.off()

#save gene expression matrix for ComplexHeatmap
mat <- subset(dt.sct, downsample = 100)[["integrated"]]@data[mkrs, ] %>% as.matrix()
## scale the rows
mat <- t(scale(t(mat)))
write.csv(mat, "output/files/micro_hmp_mat.csv")

print("microglia heatmap saved...")
print(gc())
timestamp()

###Astrocyte makers
dt.sct <- readRDS("output/seurat_astro_subset.rds")

print("dataset is subset...")
print(gc())
timestamp()

#astrocyte markers
mkrs <- character()
for (i in Sys.glob("output/files/astrocyte subset markers/*")) {
  print(i)
  tmp <- read.csv(i)
  
  nonsense <- c(tmp$X[grep("^RP", tmp$X)],
                tmp$X[grep("^Gm[0-9]", tmp$X)],
                tmp$X[grep("^AC[0-9]", tmp$X)],
                tmp$X[grep("[0-9]Rik", tmp$X)])
  
  mkrs <- union(mkrs, tmp %>%
                  filter(p_val_adj < 1e-05, avg_log2FC > 2, pct.1 > 0.7, (pct.1-pct.2) > -0.1, !X %in% nonsense) %>%
                  arrange(desc(avg_log2FC)) %>%
                  #slice(1:25) %>%
                  pull(X))
}

#Heatmap
p <- DoHeatmap(subset(dt.sct, downsample = 100), features = mkrs, size = 3, group.by = "Treatment")
pdf("output/plots/heatmaps/astro_hmp.pdf", width = 11, height = 8.5)
p
dev.off()

#save gene expression matrix for ComplexHeatmap
mat <- subset(dt.sct, downsample = 100)[["integrated"]]@data[mkrs, ] %>% as.matrix()
## scale the rows
mat <- t(scale(t(mat)))
write.csv(mat, "output/files/astro_hmp_mat.csv")

print("astrocyte heatmap saved...")
print(gc())
timestamp()

##Neurons analysis
dt.sct <- readRDS("output/seurat_neuron_subset_2.rds")

print("dataset is subset...")
print(gc())
timestamp()

#neuron markers
mkrs <- character()
for (i in Sys.glob("output/files/neuron subcluster markers/*")) {
  print(i)
  tmp <- read.csv(i)
  
  nonsense <- c(tmp$X[grep("^RP", tmp$X)],
                tmp$X[grep("^Gm[0-9]", tmp$X)],
                tmp$X[grep("^AC[0-9]", tmp$X)],
                tmp$X[grep("[0-9]Rik", tmp$X)])
  
  mkrs <- union(mkrs, tmp %>% 
                  filter(p_val_adj < 0.05, avg_log2FC > 0.5, pct.1 > 0.5, pct.2 > 0.2, (pct.1-pct.2) > 0.2, !X %in% nonsense) %>% 
                  arrange(desc(pct.1), desc(pct.1-pct.2), p_val_adj) %>%
                  slice(1:50) %>%
                  pull(X))
}

#Heatmap
p <- DoHeatmap(subset(dt.sct, downsample = 100), features = mkrs, size = 3, group.by = "Treatment")
pdf("output/plots/heatmaps/neuron_hmp.pdf", width = 11, height = 8.5)
p
dev.off()

#save gene expression matrix for ComplexHeatmap
mat <- subset(dt.sct, downsample = 100)[["integrated"]]@data[mkrs, ] %>% as.matrix()
## scale the rows
mat <- t(scale(t(mat)))
write.csv(mat, "output/files/neuron_hmp_mat.csv")

##Neurons analysis
dt.sct <- readRDS("output/seurat_neuron_subset_2.rds")

dt.sct <- subset(dt.sct,
                 idents = c("ex.neu.9", "inh.neu.10", "inh.neu.13", "ex.neu.15")
)

print("dataset is subset...")
print(gc())
timestamp()

#neuron markers
mkrs <- character()
tmp <- read.csv("output/files/Neu DE res/response_9_15_vs_10_13.csv")
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
              tmp$X[grep("^Cox", tmp$X)],
              tmp$X[grep("^Uqcr[1-10]", tmp$X)])

mkrs <- mkrs %>% 
  union(tmp  %>% 
          filter(p_val_adj < 0.05, avg_log2FC > 0.5, pct.1 > 0.5, pct.2 > 0.25, (pct.1-pct.2) > -0.1, !X %in% nonsense) %>% 
          arrange(desc(avg_log2FC)) %>%
          pull(X)) %>%
  union(tmp  %>% 
          filter(p_val_adj < 0.05, avg_log2FC < -0.5, pct.1 > 0.25, pct.2 > 0.5, (pct.1-pct.2) < 0.1, !X %in% nonsense) %>% 
          arrange(avg_log2FC) %>%
          pull(X))

#save gene expression matrix for ComplexHeatmap
mat <- dt.sct[["integrated"]]@data[mkrs, ] %>% as.matrix()
## scale the rows
mat <- t(scale(t(mat)))
write.csv(mat, "output/files/neuron_9_15_vs_10_13_hmp_mat.csv")

print("neuron heatmaps saved... job complete")
print(gc())
timestamp()