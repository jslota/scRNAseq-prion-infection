require(Seurat)
require(sctransform)
require(patchwork)
require(dplyr)

timestamp()

parallel::mclapply

dt.sct <- readRDS("output/seurat_neuron_subset_2.rds")
print("dataset is loaded...")
print(gc())
timestamp()

pdf("output/plots/neurons/neuron_featureplots.pdf", width=8, height = 6)
FeaturePlot(
  dt.sct,
  c("Prnp", "Mki67", "Reln", "Dcx", "Cck",
    "Slc17a7", "Gad1", "Rbfox3", "Calb1"),
  dims = c(1,2)
)
dev.off()

pdf("output/plots/neurons/neuron_featureplots_validate_mkrs.pdf", width=9, height = 3)
FeaturePlot(
  dt.sct,
  c("Rbfox3", "Slc17a7", "Gad1"),
  ncol=3,
  dims = c(1,2)
)
dev.off()

pdf("output/plots/neurons/neuron_featureplots_split.pdf", width=8, height = 16)
FeaturePlot(
  dt.sct,
  c("Nrgn", "Prnp", "Tafa1", "Foxp2", "Il1rapl2", "Dlx6os1", "Lhx1"),
  dims = c(1,2),
  split.by = "Treatment"
)
dev.off()

dt.sct <- subset(dt.sct,
                 idents = c("m.neu.0", "ex.neu.4", "ex.neu.6", "ex.neu.8",
                            "ex.neu.9", "inh.neu.10", "m.neu.11", "inh.neu.13", "ex.neu.15")
)

pdf("output/plots/neurons/neuron_marker_violinplot.pdf", width = 11, height = 8.5)
VlnPlot(dt.sct,
        features = c("Rbfox3", "Slc17a7", "Slc17a6", "Gad1", "Gad2", "Prnp"),
        group.by = "celltype",  pt.size = 0)
dev.off()

pdf("output/plots/neurons/neuron_violinplot_validate_mrks.pdf", width = 8, height = 4)
VlnPlot(dt.sct,
        features = c("Rbfox3", "Slc17a7", "Gad1"),
        group.by = "celltype",  pt.size = 0)
dev.off()

dt.sct <- readRDS("output/seurat_astro_subset.rds")
print("dataset is loaded...")
print(gc())
timestamp()

dt.sct$celltype <- Idents(dt.sct)

pdf("output/plots/astrocytes/astrocyte_subset_clusters.pdf", width=6, height = 5)
DimPlot(dt.sct, label = TRUE)
dev.off()

pdf("output/plots/astrocytes/astrocyte_subset_clusters_split.pdf", width=11, height = 6)
DimPlot(dt.sct, label = TRUE, split.by = "Treatment")
dev.off()

pdf("output/plots/astrocytes/astrocyte_marker_violinplot.pdf", width = 11, height = 12)
VlnPlot(dt.sct,
        features = c("Gfap", "Vim", "Aqp4", "S100b",
                    "Nrxn1", "Nrxn3", "Slc1a3", "Slc1a2", 
                    "S100a6", "S100a1", "Prdx1", "Hopx",
                    "Glis3", "Cadm1", "Zbtb20", "Maml2",
                    "Prnp",  "Clic6", "Flt1", "Cst7"),
        group.by = "celltype",  pt.size = 0)
dev.off()

print("job compelete...")
print(gc())
timestamp()