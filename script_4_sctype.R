###Cell type classification with sctype
#2022-06-06
#Jessy Slota

# load libraries and functions
lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# get cell-type-specific gene sets from our in-built database (DB)
gs_list = gene_sets_prepare("ScTypeDB_custom.xlsx", "Brain") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

parallel::mclapply
timestamp()

# assign cell types
dt.combined.sct = readRDS("output/seurat_integrated_clustered.rds") #load example scRNA-seq matrix
es.max = sctype_score(scRNAseqData = dt.combined.sct[["integrated"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

print("celltypes assigned")
timestamp()
print(gc())

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(dt.combined.sct@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(dt.combined.sct@meta.data[dt.combined.sct@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(dt.combined.sct@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])
write.csv(sctype_scores[order(sctype_scores$cluster),], "output/files/SCtype_scores.csv", row.names = FALSE)

dt.combined.sct@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  dt.combined.sct@meta.data$customclassif[dt.combined.sct@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

pdf("output/plots/integrated_datasets_sctype.pdf", width=11, height = 8.5)
DimPlot(dt.combined.sct, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif') 
dev.off()
timestamp()
print(gc())

saveRDS(dt.combined.sct, "output/seurat_integrated_clustered_sctype.rds")
print("dataset is saved... job complete")
timestamp()