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
                 idents = c("astro.10", "astro.20")
)

print("dataset is subset...")
print(gc())
timestamp()

dt.sct <- RunPCA(dt.sct)
dt.sct <- RunUMAP(dt.sct, reduction = "pca", dims = 1:50)
dt.sct <- FindNeighbors(dt.sct, reduction = "pca", dims = 1:50)
dt.sct <- FindClusters(dt.sct, resolution = 0.8)

print("UMAP complete...")
print(gc())
timestamp()

# pdf("output/plots/astrocyte_subset_clusters.pdf", width=6, height = 5)
# DimPlot(dt.sct, label = TRUE)
# dev.off()
# 
# pdf("output/plots/astrocyte_subset_clusters_split.pdf", width=11, height = 6)
# DimPlot(dt.sct, label = TRUE, split.by = "Treatment")
# dev.off()

## extract meta data
md <- dt.sct@meta.data
saveRDS(md, "output/seurat_metadata_astrocyte_subset.rds")
print("metadata is saved...")
timestamp()
print(gc())

print("finding all astrocyte markers...")
print(gc())
timestamp()

clusters <-  levels(Idents(dt.sct))

for (i in clusters) {
  print(paste0("Finding markers for ", i))
  
  response <- FindMarkers(dt.sct, ident.1 = i)
  write.csv(response, paste0("output/files/astrocyte subset markers/", i, "_markers.csv"))
}

print("saving dataset...")
saveRDS(dt.sct, "output/seurat_astro_subset.rds")
print("dataset is saved...")

print("... job complete")
print(gc())
timestamp()