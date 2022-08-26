###analysis of microglia cell clusters
#2022-06-01
#Jessy Slota

require(Seurat)
require(sctransform)
require(patchwork)
require(dplyr)

timestamp()

parallel::mclapply

dt.sct <- readRDS("output/seurat_integrated_clustered_sctype.rds")
print("dataset is loaded...")
print(gc())
timestamp()

dt.sct <- RenameIdents(dt.sct, `0` = "micro.0", `1` = "micro.1", `2` = "endo.2",
                       `3` = "micro.3", `4` = "micro.4", `5` = "micro.5", `6` = "micro.6",
                       `7` = "micro.7", `8` = "d.neu.8", `9` = "micro.9", `10` = "astro.10",
                       `11` = "endo.11", `12` = "micro.12", `13` = "micro.13", `14` = "peri.14",
                       `15` = "pvm.15", `16` = "g.neu.16", `17` = "micro.17", `18` = "im.neu.18",
                       `19` = "opc.19", `20` = "astro.20", `21` = "smc.21", `22` = "pvm.22",
                       `23` = "micro.23", `24` = "im.neu.24", `25` = "micro.25", `26` = "epen.26",
                       `27` = "endo.27", `28` = "lymph.28", `29` = "pvm.29", `30` = "vlmc.30",
                       `31` = "g.neu.31", `32` = "im.neu.32", `33` = "g.neu.33", `34` = "g.neu.34",
                       `35` = "d.neu.35", `36` = "micro.36", `37` = "opc.37", `38` = "endo.38")

dt.sct <- subset(dt.sct,
                 idents = c("micro.0", "micro.1", "endo.2",
                            "micro.3", "micro.4", "micro.5", "micro.6",
                            "micro.7", "micro.9", "astro.10",
                            "endo.11", "micro.12", "micro.13", "peri.14",
                            "g.neu.16", "micro.17", "im.neu.18",
                            "opc.19", "astro.20", "smc.21", "pvm.22",
                            "micro.23", "im.neu.24", "micro.25", "epen.26",
                            "endo.27", "lymph.28", "pvm.29", "vlmc.30",
                            "g.neu.31", "im.neu.32", "g.neu.33", "g.neu.34",
                            "micro.36", "opc.37", "endo.38")
                 )

print("dataset is subset...")
print(gc())
timestamp()

#Filter expression of technical genes
dt.sct <- dt.sct[!grepl("Malat1", rownames(dt.sct)), ]
dt.sct <- dt.sct[!grepl("^mt-", rownames(dt.sct)), ]
dt.sct <- dt.sct[!grepl('^Rp[sl]', rownames(dt.sct)), ]
dt.sct <- dt.sct[!grepl("^Hb[^(p)]", rownames(dt.sct)), ]
dt.sct <- dt.sct[!grepl("^Ndu", rownames(dt.sct)), ]
dt.sct <- dt.sct[!grepl("^Cox", rownames(dt.sct)), ]
dt.sct <- dt.sct[!grepl("^Uqcr[1-10]", rownames(dt.sct)), ]

dt.sct <- RunPCA(dt.sct)
dt.sct <- RunUMAP(dt.sct, reduction = "pca", dims = 1:50)

print("UMAP complete...")
print(gc())
timestamp()

pdf("output/plots/filtered_clusters.pdf", width=11, height = 8.5)
DimPlot(dt.sct, label = TRUE)
dev.off()

pdf("output/plots/filtered_clusters_split.pdf", width=11, height = 8.5)
DimPlot(dt.sct, label = TRUE, split.by = "Treatment")
dev.off()

pdf("output/plots/filtered_clusters_featureplots.pdf", width=14, height = 8.5)
FeaturePlot(
  dt.sct,
  c("P2ry12", "Gfap", "Rbfox3", "Cd163", "Pdgfra", "Spag17", "Dcx", "Mki67", "Cldn5"),
  dims = c(1,2)
)
dev.off()

clusters <- levels(Idents(dt.sct))

dt.sct$celltype.Treatment <- paste(Idents(dt.sct), dt.sct$Treatment, sep = "_")
dt.sct$celltype <- Idents(dt.sct)
Idents(dt.sct) <- "celltype.Treatment"

## extract meta data
md <- dt.sct@meta.data
saveRDS(md, "output/seurat_metadata_subset.rds")
print("metadata is saved...")
timestamp()
print(gc())

print("calculating DE genes for each cluster... ")
print(clusters)
timestamp()
print(gc())

#DE analysis for each cluster
for (i in clusters) {
  if (nrow(md[md$celltype.Treatment==paste0(i, "_RML"),])>3 && nrow(md[md$celltype.Treatment==paste0(i, "_Mock"),])>3) {
    print(paste0("calculating DE genes for cluster ", i))
    print(paste0(nrow(md[md$celltype.Treatment==paste0(i, "_RML"),]), " cells in cluster ", paste0(i, "_RML")))
    print(paste0(nrow(md[md$celltype.Treatment==paste0(i, "_Mock"),]), " cells in cluster ", paste0(i, "_Mock")))
    timestamp()
    response <- FindMarkers(dt.sct, ident.1 = paste0(i, "_RML"), ident.2 = paste0(i, "_Mock"))
    write.csv(response, paste0("output/files/DE res/response_", i, ".csv"))
  } else  if (nrow(md[md$celltype.Treatment==paste0(i, "_RML"),])<3) {
    print(paste0("Skipping DE analysis in cluster ", i, "... No cells in cluster ", paste0(i, "_RML")))
  } else  if (nrow(md[md$celltype.Treatment==paste0(i, "_Mock"),])<3) {
    print(paste0("Skipping DE analysis in cluster ", i, "... No cells in cluster ", paste0(i, "_Mock")))
  }
}

print("finished DE analysis")
print(gc())
timestamp()

print("saving dataset...")
saveRDS(dt.sct, "output/seurat_final_filtered.rds")
print("dataset is saved... job complete")
timestamp()
