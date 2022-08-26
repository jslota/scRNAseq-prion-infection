require(Seurat)
require(sctransform)
require(patchwork)
require(dplyr)
require(viridisLite)

timestamp()

parallel::mclapply

dt.sct <- readRDS("output/seurat_micro_subset.rds")
print("dataset is loaded...")
print(gc())
timestamp()

pdf("output/plots/microglia/microglia_clusters.pdf", width=11, height = 8.5)
DimPlot(dt.sct, label = TRUE)
dev.off()

pdf("output/plots/microglia/microglia_clusters_split.pdf", width=11, height = 6)
DimPlot(dt.sct, label = TRUE,
        split.by = "Treatment",
        raster = TRUE,
        raster.dpi = c(1200, 1200))
dev.off()

pdf("output/plots/microglia/microglia_featureplots.pdf", width=14, height = 8.5)
FeaturePlot(
  dt.sct,
  c("P2ry12", "Cx3cr1", "Tmem119", "Aif1", 
    "C1qa", "Trem2", "Tyrobp", "Foxp1", 
    "Rhob",  "Il1b", "Serpine1", "Il12b",
    "Cxcl10", "Trim30a", "Cxcl2", "Cd74"),
  dims = c(1,2),
  raster = TRUE,
  raster.dpi = c(1200, 1200)
)
dev.off()

pdf("output/plots/microglia/microglia_featureplots_validate_mkrs.pdf", width=8, height = 6)
FeaturePlot(
  dt.sct,
  c("Tmem119", "Aif1", "Nav2", "Fos", "Cxcl10", "Ftl1", "Cd74"),
  dims = c(1,2),
  raster = TRUE,
  raster.dpi = c(1200, 1200)
)
dev.off()

pdf("output/plots/microglia/microglia_marker_violinplot.pdf", width = 14, height = 18)
VlnPlot(dt.sct,
        features = c("Nav2", "P2ry12", "Cx3cr1", "Tmem119",
                     "Jun", "Fos", "Il1a", "Il12b",
                     "Trim30a", "Oasl2", "Cxcl10", "Bst2",
                     "Aif1", "Ftl1", "Fau", "Trem2",
                     "Cd74", "H2-Aa", "Cd52", "Ccl6",
                     "Serpine1", "Cst7", "Cd14", "Nfkbia"),
        group.by = "celltype",  pt.size = 0)
dev.off()

pdf("output/plots/microglia/microglia_violinplot_validate_mkrs.pdf", width = 9, height = 9)
VlnPlot(dt.sct,
        features = c("Tmem119", "Aif1", "Nav2", "Fos", "Cxcl10", "Ftl1", "Cd74"),
        group.by = "celltype",  pt.size = 0)
dev.off()

pdf("output/plots/microglia/microglia_featureplots_astro_mkrs.pdf", width=10, height = 8)
FeaturePlot(
  dt.sct,
  c("Gfap", "Aqp4", "Vim", "Apod", "Serpina3n", "Serping1", "Timp1", "S100b", "Aldoc", "Nes", "Nrxn1"),
  dims = c(1,2),
  raster = TRUE,
  raster.dpi = c(1200, 1200)
)
dev.off()

dt.sct <- RenameIdents(dt.sct, `micro.0` = "prolif", `micro.1` = "inter", `micro.3` = "homeo", `micro.4` = "homeo",
                       `micro.5` = "inter", `micro.6` = "inter", `micro.7` = "inter", `micro.9` = "prolif",
                       `micro.12` = "prolif", `micro.13` = "phago", `micro.17` = "MHC", `micro.23` = "IFN",
                       `micro.25` = "homeo", `micro.36` = "MHC")

pdf("output/plots/trajectory plots/microglia_subtypes.pdf", width=6, height = 5)
DimPlot(dt.sct, cols = viridis(6),
        label = TRUE,
        raster = TRUE,
        raster.dpi = c(1200, 1200))
dev.off()

print("plots saved... job complete")
print(gc())
timestamp()