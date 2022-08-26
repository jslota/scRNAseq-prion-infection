###analysis of microglia cell clusters
#2022-06-01
#Jessy Slota

require(Seurat)
require(sctransform)
require(patchwork)
require(dplyr)

timestamp()

parallel::mclapply

dt.sct <- readRDS("output/seurat_final_filtered.rds")
print("dataset is loaded...")
print(gc())
timestamp()

print(levels(Idents(dt.sct)))
Idents(dt.sct) <- "celltype"
print(levels(Idents(dt.sct)))
dt.sct <- subset(dt.sct,
                 idents = c("micro.0", "micro.1", "micro.3", "micro.4", "micro.5", "micro.6",
                            "micro.7", "micro.9", "micro.12", "micro.13", "micro.17", "micro.23",
                            "micro.25", "micro.36")
                 )

print("dataset is subset...")
print(gc())
timestamp()

dt.sct <- RunPCA(dt.sct)
dt.sct <- RunUMAP(dt.sct, reduction = "pca", dims = 1:50)

print("UMAP complete...")
print(gc())
timestamp()

pdf("output/plots/microglia_clusters.pdf", width=11, height = 8.5)
DimPlot(dt.sct, label = TRUE)
dev.off()

pdf("output/plots/microglia_clusters_split.pdf", width=11, height = 8.5)
DimPlot(dt.sct, label = TRUE, split.by = "Treatment")
dev.off()

print("finding all microglia markers...")
print(gc())
timestamp()


dam_clusters <- c("micro.0", "micro.1", "micro.3", "micro.4", "micro.5", "micro.6",
                  "micro.7", "micro.9", "micro.12", "micro.13", "micro.17", "micro.23",
                  "micro.25", "micro.36")
for (i in dam_clusters) {
  print(paste0("Finding markers for ", i))
  response <- FindMarkers(dt.sct, ident.1 = i)
  write.csv(response, paste0("output/files/microglia subcluster markers/", i, "_markers.csv"))
}

print("markers found... job complete")
print(gc())
timestamp()
