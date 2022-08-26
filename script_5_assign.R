###assign clusters on integrated dataset
#2022-06-01
#Jessy Slota

require(Seurat)
require(sctransform)
require(patchwork)
require(dplyr)

timestamp()

parallel::mclapply

dt.combined.sct <- readRDS("output/seurat_integrated_clustered_sctype.rds")
print("dataset is loaded...")
print(gc())
timestamp()

dt.combined.sct <- RenameIdents(dt.combined.sct, `0` = "micro.0", `1` = "micro.1", `2` = "endo.2",
                                `3` = "r.micro.3", `4` = "r.micro.4", `5` = "r.micro.5", `6` = "micro.6",
                                `7` = "r.micro.7", `8` = "d.neu.8", `9` = "dc.9", `10` = "astro.10",
                                `11` = "endo.11", `12` = "micro.12", `13` = "r.micro.13", `14` = "peri.14",
                                `15` = "pvm.15", `16` = "g.neu.16", `17` = "dc.17", `18` = "im.neu.18",
                                `19` = "opc.19", `20` = "astro.20", `21` = "smc.21", `22` = "pvm.22",
                                `23` = "r.micro.23", `24` = "im.neu.24", `25` = "micro.25", `26` = "astro.26",
                                `27` = "endo.27", `28` = "lymph.28", `29` = "dc.29", `30` = "vlmc.30",
                                `31` = "g.neu.31", `32` = "im.neu.32", `33` = "g.neu.33", `34` = "g.neu.34",
                                `35` = "d.neu.35", `36` = "dc.36", `37` = "oligo.37", `38` = "endo.38")

pdf("output/plots/assigned_clusters.pdf", width=11, height = 8.5)
DimPlot(dt.combined.sct, label = TRUE)
dev.off()

pdf("output/plots/assigned_clusters_split.pdf", width=11, height = 8.5)
DimPlot(dt.combined.sct, label = TRUE, split.by = "Treatment")
dev.off()

clusters <- levels(Idents(dt.combined.sct))
print(clusters)

dt.combined.sct$celltype.Treatment <- paste(Idents(dt.combined.sct), dt.combined.sct$Treatment, sep = "_")
dt.combined.sct$celltype <- Idents(dt.combined.sct)
Idents(dt.combined.sct) <- "celltype.Treatment"

## extract meta data
md <- dt.combined.sct@meta.data
saveRDS(md, "output/seurat_metadata_subset.rds")
print("metadata is saved...")
timestamp()
print(gc())

write.csv(table(Idents(dt.combined.sct)), "output/files/number_cells_cluster.csv")

print("finished analysis... job complete")
print(gc())
timestamp()

