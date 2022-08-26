###Identify clusters on integrated dataset
#2022-06-01
#Jessy Slota

require(Seurat)
require(sctransform)
require(patchwork)
require(dplyr)

timestamp()

parallel::mclapply

dt.combined.sct <- readRDS("output/seurat_integrated.rds")

print("dataset is loaded")
timestamp()
print(gc())

dt.combined.sct <- RunPCA(dt.combined.sct)
dt.combined.sct <- RunUMAP(dt.combined.sct, reduction = "pca", dims = 1:50)
dt.combined.sct <- FindNeighbors(dt.combined.sct, reduction = "pca", dims = 1:50)
dt.combined.sct <- FindClusters(dt.combined.sct)

print("data is clustered")
timestamp()
print(gc())

p1 <- DimPlot(dt.combined.sct, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(dt.combined.sct, reduction = "umap", label = TRUE, repel = TRUE)
pdf("output/plots/integrated_datasets.pdf", width=11, height = 8.5)
p1 + p2
dev.off()

all.markers <- FindAllMarkers(dt.combined.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers %>% write.csv("output/files/all_markers.csv")
all.markers %>%
  group_by(cluster) %>%
  slice_max(n = 25, order_by = avg_log2FC) %>%
  write.csv("output/files/top_25_markers.csv")

print("markers found")
timestamp()
print(gc())

saveRDS(dt.combined.sct, "output/seurat_integrated_clustered.rds")

print("dataset is saved... job complete")
print(gc())
timestamp()