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
                 idents = c("g.neu.16", "im.neu.18", "im.neu.24", "g.neu.31",
                            "im.neu.32", "g.neu.33", "g.neu.34")
)

print("dataset is subset...")
print(gc())
timestamp()

dt.sct <- RunPCA(dt.sct)
dt.sct <- RunUMAP(dt.sct, reduction = "pca", dims = 1:50)
dt.sct <- FindNeighbors(dt.sct, reduction = "pca", dims = 1:50)
dt.sct <- FindClusters(dt.sct, resolution = 0.4)

print("UMAP complete...")
print(gc())
timestamp()

pdf("output/plots/neuron_clusters.pdf", width=11, height = 8.5)
DimPlot(dt.sct, label = TRUE)
dev.off()

pdf("output/plots/neuron_clusters_split.pdf", width=11, height = 8.5)
DimPlot(dt.sct, label = TRUE, split.by = "Treatment")
dev.off()

## extract meta data
md <- dt.sct@meta.data
saveRDS(md, "output/seurat_metadata_neuron_subsets.rds")
print("metadata is saved...")
timestamp()
print(gc())

print("finding all neuron markers...")
print(gc())
timestamp()

clusters <-  levels(Idents(dt.sct))

for (i in clusters) {
  print(paste0("Finding markers for ", i))

  response <- FindMarkers(dt.sct, ident.1 = i)
  write.csv(response, paste0("output/files/neuron subcluster markers/", i, "_markers.csv"))
}

print("saving dataset...")
saveRDS(dt.sct, "output/seurat_neuron_subset.rds")
print("dataset is saved... job complete")

print("markers found... job complete")
print(gc())
timestamp()
