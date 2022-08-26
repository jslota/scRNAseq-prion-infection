
#trying to integrate two samples with soupX

library(SoupX)
library(Seurat)
library(sctransform)
library(dplyr)
library(HGNChelper)
library(DoubletFinder)

parallel::mclapply
timestamp()

fls <- c("/Drives/L/HGPD/Jessy_Slota/10x_genomics/cell_ranger_runs/PBS48CX_count_res/outs/",
         "/Drives/L/HGPD/Jessy_Slota/10x_genomics/cell_ranger_runs/PBS60CX_count_res_merged/outs/",
         "/Drives/L/HGPD/Jessy_Slota/10x_genomics/cell_ranger_runs/PBS61CX_count_res/outs/",
         "/Drives/L/HGPD/Jessy_Slota/10x_genomics/cell_ranger_runs/PBS73CX_count_res/outs/",
         "/Drives/L/HGPD/Jessy_Slota/10x_genomics/cell_ranger_runs/RML122CX_count_res_merged/outs/",
         "/Drives/L/HGPD/Jessy_Slota/10x_genomics/cell_ranger_runs/RML132CX_count_res_merged/outs/",
         "/Drives/L/HGPD/Jessy_Slota/10x_genomics/cell_ranger_runs/RML133CX_count_res_merged/outs/",
         "/Drives/L/HGPD/Jessy_Slota/10x_genomics/cell_ranger_runs/RML134CX_count_res_merged/outs/",
         "/Drives/L/HGPD/Jessy_Slota/10x_genomics/cell_ranger_runs/RML138CX_count_res_merged/outs/",
         "/Drives/L/HGPD/Jessy_Slota/10x_genomics/cell_ranger_runs/RML140CX_count_res_merged/outs/",
         "/Drives/L/HGPD/Jessy_Slota/10x_genomics/cell_ranger_runs/RML142CX_count_res_merged/outs/",
         "/Drives/L/HGPD/Jessy_Slota/10x_genomics/cell_ranger_runs/RML145CX_count_res_merged/outs/",
         "/Drives/L/HGPD/Jessy_Slota/10x_genomics/cell_ranger_runs/PBS25HP_count_res/outs/",
         "/Drives/L/HGPD/Jessy_Slota/10x_genomics/cell_ranger_runs/PBS48HP_count_res/outs/",
         "/Drives/L/HGPD/Jessy_Slota/10x_genomics/cell_ranger_runs/RML122HP_count_res_merged/outs/",
         "/Drives/L/HGPD/Jessy_Slota/10x_genomics/cell_ranger_runs/RML132HP_count_res/outs/",
         "/Drives/L/HGPD/Jessy_Slota/10x_genomics/cell_ranger_runs/RML133HP_count_res_merged/outs/",
         "/Drives/L/HGPD/Jessy_Slota/10x_genomics/cell_ranger_runs/RML138HP_count_res_merged/outs/",
         "/Drives/L/HGPD/Jessy_Slota/10x_genomics/cell_ranger_runs/RML140HP_count_res/outs/",
         "/Drives/L/HGPD/Jessy_Slota/10x_genomics/cell_ranger_runs/RML142HP_count_res_merged/outs/",
         "/Drives/L/HGPD/Jessy_Slota/10x_genomics/cell_ranger_runs/RML145HP_count_res/outs/")

for (i in fls) {
  print(i)
  x <- gsub("_.*", "", gsub("/Drives/L/HGPD/Jessy_Slota/10x_genomics/cell_ranger_runs/", "", i))
  print(paste0("beginning QC pre-processing... ", x))
  
  sc <- load10X(i)
  sc <- autoEstCont(sc)
  out <- adjustCounts(sc)
  
  print(paste0("soup removal complete... ", x))
  timestamp()
  print(gc())
  
  srat = CreateSeuratObject(out, project = x)
  rm(out, sc)
  print(gc())
  
  print(paste0("initial...", x))
  print(dim(srat))
  
  #percent mitochondrial reads to regress out
  srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^mt-")
  #Add percent ribo-protein reads
  srat <- PercentageFeatureSet(srat, "^Rp[sl]", col.name = "percent.ribo")
  # Percentage hemoglobin genes - includes all genes starting with HB except HBP.
  srat <- PercentageFeatureSet(srat, "^Hb[^(p)]", col.name = "percent.hb")
  
  #filter the dataset
  srat <- subset(srat, subset = nFeature_RNA > 1000 & percent.mt < 20 & percent.ribo > 1 & percent.hb < 20)
  
  #Apply sctransform normalization 
  #BiocManager::install("glmGamPoi")
  srat <- SCTransform(srat, method = "glmGamPoi", vars.to.regress = "percent.mt")
  
  print(paste0("dataset is normalized... ", x))
  timestamp()
  print(gc())
  
  # These are now standard steps in the Seurat workflow for visualization and clustering
  srat <- RunPCA(srat)
  srat <- RunUMAP(srat, dims = 1:30)
  
  # define the expected number of doublet cells.
  nExp <- round(ncol(srat) * 0.076)  # expect 7.6% doublets
  srat <- doubletFinder_v3(srat, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10, sct = TRUE)
  
  # name of the DF prediction can change, so extract the correct column name.
  DF.name = colnames(srat@meta.data)[grepl("DF.classification", colnames(srat@meta.data))]
  
  #remove all predicted doublets
  srat = srat[, srat@meta.data[, DF.name] == "Singlet"]
  print(paste0("after doublets removed... ", x))
  print(dim(srat))
  
  saveRDS(srat, paste0("datasets/", x, ".rds"))
  print(paste0("dataset saved... ", x))
  print(gc())
}

print("job complete... all datasets are qc'd")