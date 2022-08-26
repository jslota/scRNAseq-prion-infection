library(Seurat)
library(DropletUtils)
library(Matrix)
library(dplyr)

timestamp()
parallel::mclapply

dt.sct <- readRDS("output/seurat_final_filtered.rds")
print("dataset is loaded...")
print(gc())
timestamp()

write10xCounts(x = dt.sct@assays$RNA@counts, path = "output/files_upload/raw_matrix/")
write10xCounts(x = GetAssayData(dt.sct, assay = "integrated", slot = "data" ), path = "output/files_upload/normalized_matrix/")

md <- readRDS("output/seurat_metadata_subset.rds")
tmp <- (gsub("\\..*", "", md$celltype))
unique(tmp)

out <- data.frame(NAME = rownames(md),
                  biosample_id = gsub("RML", "mouse", gsub("PBS", "mouse", md$orig.ident)),
                  donor_id = gsub("HP", "", gsub("CX", "", gsub("RML", "mouse", gsub("PBS", "mouse", md$orig.ident)))),
                  species = "NCBI:txid10090",
                  species__ontology_label = "Mus musculus",
                  disease = gsub("Mock", "PATO_0000461", gsub("RML", "MONDO:0005429", md$Treatment)),
                  disease__ontology_label = gsub("Mock", "normal", gsub("RML", "prion disease", md$Treatment)),
                  organ = gsub("cx", "UBERON:0001851", gsub("hp", "UBERON:0002421", md$Region)),
                  organ__ontology_label = gsub("cx", "cortex", gsub("hp", "hippocampal formation", md$Region)),
                  library_preparation_protocol = "EFO:0030003",
                  library_preparation_protocol__ontology_label = "10x 3' transcription profiling",
                  sex = "female",
                  cell_type =gsub("lymph", "CL:0000542",
                                  gsub("opc", "CL:0002453",
                                       gsub("vlmc", "CL:4023051",
                                            gsub("pvm", "CL:0000881",
                                                 gsub("peri", "CL:0000669",
                                                      gsub("g", "CL:0000540",
                                                           gsub("im", "PATO:0001501",
                                                                gsub("epen", "CL:0000065",
                                                                     gsub("endo", "CL:2000044",
                                                                          gsub("astro", "CL:0000127",
                                                                               gsub("smc", "CL:0002590",
                                                                                    gsub("micro", "CL:0000129", tmp)))))))))))),
                  cell_type__ontology_label = gsub("endo", "brain microvascular endothelial cell",
                                                   gsub("micro", "microglial cell",
                                                        gsub("opc", "oligodendrocyte precursor cell",
                                                             gsub("vlmc", "vascular leptomeningeal cell",
                                                                  gsub("pvm", "perivascular macrophage",
                                                                       gsub("peri", "pericyte",
                                                                            gsub("g", "neuron",
                                                                                 gsub("im", "immature",
                                                                                      gsub("epen", "ependymal cell",
                                                                                           gsub("astro", "astrocyte",
                                                                                                gsub("smc", "smooth muscle cell of the brain vasculature",
                                                                                                     gsub("lymph", "lymphocyte", tmp)))))))))))),
                  Treatment = md$Treatment,
                  Region = md$Region,
                  seurat_clusters = md$seurat_clusters,
                  celltype_cluster = md$celltype,
                  celltype_cluster_treatment = md$celltype.Treatment)

write.csv(out, "output/files_upload/metadata.csv", row.names = FALSE)

dt.sct <- RunUMAP(dt.sct, dims = 1:50, n.components = 3)

out <- dt.sct[["umap"]]@cell.embeddings %>% as.data.frame()
write.csv(out, "output/files_upload/umap_embeddings.csv")

print("files saved... job complete")
print(gc())
timestamp()