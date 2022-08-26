###analysis of microglia cell clusters
#2022-06-01
#Jessy Slota

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

clusters <- levels(Idents(dt.sct))

dt.sct$celltype.Treatment <- paste(Idents(dt.sct), dt.sct$Treatment, sep = "_")
dt.sct$celltype <- Idents(dt.sct)
Idents(dt.sct) <- "celltype.Treatment"

## extract meta data
md <- dt.sct@meta.data


#DE analysis for each cluster
for (i in clusters) {
  if (nrow(md[md$celltype.Treatment==paste0(i, "_RML"),])>3 && nrow(md[md$celltype.Treatment==paste0(i, "_Mock"),])>3) {
    print(paste0("calculating DE genes for cluster ", i))
    print(paste0(nrow(md[md$celltype.Treatment==paste0(i, "_RML"),]), " cells in cluster ", paste0(i, "_RML")))
    print(paste0(nrow(md[md$celltype.Treatment==paste0(i, "_Mock"),]), " cells in cluster ", paste0(i, "_Mock")))
    timestamp()
    response <- FindMarkers(dt.sct, ident.1 = paste0(i, "_RML"), ident.2 = paste0(i, "_Mock"))
    write.csv(response, paste0("output/files/Neu DE res/response_", i, ".csv"))
  } else  if (nrow(md[md$celltype.Treatment==paste0(i, "_RML"),])<3) {
    print(paste0("Skipping DE analysis in cluster ", i, "... No cells in cluster ", paste0(i, "_RML")))
  } else  if (nrow(md[md$celltype.Treatment==paste0(i, "_Mock"),])<3) {
    print(paste0("Skipping DE analysis in cluster ", i, "... No cells in cluster ", paste0(i, "_Mock")))
  }
}

print("finished DE analysis")
print(gc())
timestamp()

#DE full cells in main cluster
Idents(dt.sct) <- dt.sct$celltype

dt.sct <- subset(dt.sct,
                 idents = c("m.neu.0", "ex.neu.4", "ex.neu.8", "ex.neu.9",
                            "inh.neu.10", "m.neu.11", "inh.neu.13", "ex.neu.15")
)

response <- FindMarkers(dt.sct, ident.1 = "RML", ident.2 = "Mock", group.by = "Treatment")
write.csv(response, paste0("output/files/Neu DE res/extra_response_main.csv"))


#DE compare clusters 9+15 (sustained) with 10+13 (depleted)
Idents(dt.sct) <- dt.sct$celltype

dt.sct <- subset(dt.sct,
                 idents = c("ex.neu.9", "inh.neu.10", "inh.neu.13", "ex.neu.15")
)

response <- FindMarkers(dt.sct, ident.1 = c("ex.neu.9", "ex.neu.15"), ident.2 = c("inh.neu.10", "inh.neu.13"))
write.csv(response, paste0("output/files/Neu DE res/extra_response_9_15_vs_10_13.csv"))


print("markers found... job complete")
print(gc())
timestamp()
