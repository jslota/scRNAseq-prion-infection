###analysis of microglia cell clusters
#2022-06-01
#Jessy Slota

require(Seurat)
require(sctransform)
require(patchwork)
require(dplyr)

timestamp()

parallel::mclapply

dt.sct <- readRDS("output/seurat_neuron_subset.rds")
print("dataset is loaded...")
print(gc())
timestamp()

print(levels(Idents(dt.sct)))

dt.sct <- subset(dt.sct,
                 idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8",
                            "10", "11", "12", "13", "14")
)

print("dataset is subset...")
print(gc())
timestamp()

dt.sct <- RunPCA(dt.sct)
dt.sct <- RunUMAP(dt.sct, reduction = "pca", dims = 1:50)
dt.sct <- FindNeighbors(dt.sct, reduction = "pca", dims = 1:50)
dt.sct <- FindClusters(dt.sct, resolution = 0.6)


dt.sct <- RenameIdents(dt.sct, `0` = "m.neu.0", `1` = "npc.1", `2` = "diff.neu.2",
                       `3` = "diff.neu.3", `4` = "ex.neu.4", `5` = "diff.neu.5", `6` = "ex.neu.6",
                       `7` = "cr.neu.7", `8` = "ex.neu.8", `9` = "ex.neu.9", `10` = "inh.neu.10",
                       `11` = "m.neu.11", `12` = "diff.neu.12", `13` = "inh.neu.13", `14` = "diff.neu.14",
                       `15` = "ex.neu.15")
dt.sct$celltype <- Idents(dt.sct)

print("UMAP complete...")
print(gc())
timestamp()

pdf("output/plots/neuron_subset_clusters.pdf", width=11, height = 8.5)
DimPlot(dt.sct, label = TRUE)
dev.off()

pdf("output/plots/neuron_subset_clusters_split.pdf", width=11, height = 6)
DimPlot(dt.sct, label = TRUE, split.by = "Treatment")
dev.off()

## extract meta data
md <- dt.sct@meta.data
saveRDS(md, "output/seurat_metadata_neuron_subsets_2.rds")
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
saveRDS(dt.sct, "output/seurat_neuron_subset_2.rds")
print("dataset is saved... job complete")

print("markers found... job complete")
print(gc())
timestamp()