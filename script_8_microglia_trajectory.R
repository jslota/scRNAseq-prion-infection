library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
set.seed(1234)

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

#Convert seurate to monocle3
dt.cds <- as.cell_data_set(dt.sct)

## Calculate size factors using built-in function in monocle3
dt.cds <- estimate_size_factors(dt.cds)

## Add gene names into CDS
dt.cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(dt.sct[["integrated"]])

dt.cds <- cluster_cells(cds = dt.cds, reduction_method = "UMAP")
dt.cds <- learn_graph(dt.cds)

print("dataset is converted to CellDataSet...")
print(gc())
timestamp()

hsc <- readLines("output/files/micro_4_root_cells.txt")
dt.cds <- order_cells(dt.cds, reduction_method = "UMAP", root_cells = hsc)

pdf("output/plots/trajectory plots/microglia_trajectories.pdf", width=11, height = 8.5)
plot_cells(
  cds = dt.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)
dev.off()

pdf("output/plots/trajectory plots/microglia_trajectories_small.pdf", width=6, height = 5)
plot_cells(
  cds = dt.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)
dev.off()

#add pseudotime to seurat metadata
dt.sct <- AddMetaData(
  object = dt.sct,
  metadata = dt.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "Pseudotime"
)

## extract meta data
md <- dt.sct@meta.data
saveRDS(md, "output/seurat_metadata_micro_subset.rds")
print("metadata is saved...")
timestamp()
print(gc())

print("saving dataset...")
saveRDS(dt.sct, "output/seurat_micro_subset.rds")
print("dataset is saved...")

print("... job complete")
print(gc())
timestamp()

# res <- graph_test(dt.cds, neighbor_graph="principal_graph", cores=16)
# write.csv(res, "output/files/trajectory_mkrs.csv")
# 
# print("trajectory markers saved...")
# print(gc())
# timestamp()
# 
# save_monocle_objects(cds=dt.cds, directory_path='output/microglia_monocle_cds', comment='microglia cds. Stored 2022-07-26.')
# 
# print("monocle dataset saved... trajectory analysis complete...")
# print(gc())
# timestamp()

