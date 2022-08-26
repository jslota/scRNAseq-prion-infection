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

pdf("output/plots/prnp/prnp_featureplots.pdf", width=6, height = 5)
FeaturePlot(
  dt.sct,
  c("Prnp"),
  dims = c(1,2),
  raster = TRUE,
  raster.dpi = c(1200, 1200)
)
dev.off()


pdf("output/plots/prnp/prnp_violinplot.pdf", width = 12, height = 6)
VlnPlot(dt.sct,
        features = c("Prnp"),
        group.by = "celltype",  pt.size = 0)
dev.off()

pdf("output/plots/all_mkr_violinplot.pdf", width = 11, height = 16)
VlnPlot(dt.sct,
        features = c("P2ry12", "Tmem119", "Aif1", "Gfap",
                     "Slc1a2", "Rbfox3", "Grin2b", "Dcx",
                     "Mki67", "Reln", "Cd163", "Pdgfra",
                     "Spag17", "Cldn5", "Il34", "Pln"),
        group.by = "celltype",  pt.size = 0)
dev.off()

pdf("output/plots/filtered_clusters_featureplots.pdf", width=10, height = 8)
FeaturePlot(
  dt.sct,
  c("P2ry12", "Gfap", "Rbfox3", "Cd163", "Pdgfra", "Spag17", "Dcx", "Mki67", "Cldn5"),
  dims = c(1,2)
)
dev.off()

pdf("output/plots/filtered_clusters_nonraster.pdf", width = 10, height = 8)
DimPlot(dt.sct, group.by = "celltype", label = TRUE, raster = FALSE, shuffle = TRUE) + NoLegend()
dev.off()

print("job complete...")
print(gc())
timestamp()