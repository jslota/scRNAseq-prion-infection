
require(Seurat)
require(sctransform)
require(patchwork)
require(dplyr)

timestamp()

parallel::mclapply

#load data
PBS48CX <- readRDS("datasets/PBS48CX.rds")
PBS60CX <- readRDS("datasets/PBS60CX.rds")
PBS61CX <- readRDS("datasets/PBS61CX.rds")
PBS73CX <- readRDS("datasets/PBS73CX.rds")
RML122CX <- readRDS("datasets/RML122CX.rds")
RML132CX <- readRDS("datasets/RML132CX.rds")
RML133CX <- readRDS("datasets/RML133CX.rds")
RML134CX <- readRDS("datasets/RML134CX.rds")
RML138CX <- readRDS("datasets/RML138CX.rds")
RML140CX <- readRDS("datasets/RML140CX.rds")
RML142CX <- readRDS("datasets/RML142CX.rds")
RML145CX <- readRDS("datasets/RML145CX.rds")

#load data
PBS25HP <- readRDS("datasets/PBS25HP.rds")
PBS48HP <- readRDS("datasets/PBS48HP.rds")
RML122HP <- readRDS("datasets/RML122HP.rds")
RML132HP <- readRDS("datasets/RML132HP.rds")
RML133HP <- readRDS("datasets/RML133HP.rds")
RML138HP <- readRDS("datasets/RML138HP.rds")
RML140HP <- readRDS("datasets/RML140HP.rds")
RML142HP <- readRDS("datasets/RML142HP.rds")
RML145HP <- readRDS("datasets/RML145HP.rds")


#add meta data
PBS48CX$Treatment <- "Mock"
PBS60CX$Treatment <- "Mock"
PBS61CX$Treatment <- "Mock"
PBS73CX$Treatment <- "Mock"
RML122CX$Treatment <- "RML"
RML132CX$Treatment <- "RML"
RML133CX$Treatment <- "RML"
RML134CX$Treatment <- "RML"
RML138CX$Treatment <- "RML"
RML140CX$Treatment <- "RML"
RML142CX$Treatment <- "RML"
RML145CX$Treatment <- "RML"

PBS25HP$Treatment <- "Mock"
PBS48HP$Treatment <- "Mock"
RML122HP$Treatment <- "RML"
RML132HP$Treatment <- "RML"
RML133HP$Treatment <- "RML"
RML138HP$Treatment <- "RML"
RML140HP$Treatment <- "RML"
RML142HP$Treatment <- "RML"
RML145HP$Treatment <- "RML"

PBS48CX$Region <- "cx"
PBS60CX$Region <- "cx"
PBS61CX$Region <- "cx"
PBS73CX$Region <- "cx"
RML122CX$Region <- "cx"
RML132CX$Region <- "cx"
RML133CX$Region <- "cx"
RML134CX$Region <- "cx"
RML138CX$Region <- "cx"
RML140CX$Region <- "cx"
RML142CX$Region <- "cx"
RML145CX$Region <- "cx"

PBS25HP$Region <- "hp"
PBS48HP$Region <- "hp"
RML122HP$Region <- "hp"
RML132HP$Region <- "hp"
RML133HP$Region <- "hp"
RML138HP$Region <- "hp"
RML140HP$Region <- "hp"
RML142HP$Region <- "hp"
RML145HP$Region <- "hp"

print("datasets are loaded")
print(gc())
timestamp()

dt.list <- merge(PBS48CX, c(PBS60CX, PBS61CX, PBS73CX,
                            RML122CX, RML132CX, RML133CX, RML134CX, RML138CX, RML140CX, RML142CX, RML145CX,
                            PBS25HP, PBS48HP,
                            RML122HP, RML132HP, RML133HP, RML138HP, RML140HP, RML142HP, RML145HP),
                 add.cell.ids = c("PBS48CX", "PBS60CX", "PBS61CX", "PBS73CX",
                                  "RML122CX", "RML132CX", "RML133CX", "RML134CX", "RML138CX", "RML140CX", "RML142CX", "RML145CX",
                                  "PBS25HP", "PBS48HP",
                                  "RML122HP", "RML132HP", "RML133HP", "RML138HP", "RML140HP", "RML142HP", "RML145HP"))

rm(PBS48CX,PBS60CX, PBS61CX, PBS73CX,
   RML122CX, RML132CX, RML133CX, RML134CX, RML138CX, RML140CX, RML142CX, RML145CX,
   PBS25HP, PBS48HP,
   RML122HP, RML132HP, RML133HP, RML138HP, RML140HP, RML142HP, RML145HP)
print(gc())
print("datasets merged")

# split the dataset into a list of two seurat objects (stim and CTRL)
dt.list <- SplitObject(dt.list, split.by = "orig.ident")

# normalize and identify variable features for each dataset independently
# use sctransform for normalization
dt.list <- lapply(X = dt.list, FUN = function(X) {
  X <- SCTransform(object = X, method = "glmGamPoi", vars.to.regress = c("percent.mt"))
}
)
features <- SelectIntegrationFeatures(object.list = dt.list, nfeatures = 3000)
dt.list <- PrepSCTIntegration(object.list = dt.list, anchor.features = features)

dt.list <- lapply(X = dt.list, FUN = function(x) {
  x <- RunPCA(x, features = features, verbose = FALSE)
}
)

print("dataset is normalized")
timestamp()
print(gc())

#find anchors
dt.anchors <- FindIntegrationAnchors(object.list = dt.list,
                                     normalization.method = "SCT", 
                                     anchor.features = features,
                                     reference = c(1,12),
                                     reduction = "rpca",
                                     dims = 1:50)
dt.combined.sct <- IntegrateData(anchorset = dt.anchors,
                                 normalization.method = "SCT", 
                                 dims = 1:50)

# dt.combined.sct <- ScaleData(dt.combined.sct, verbose = FALSE)
# dt.combined.sct <- RunPCA(dt.combined.sct, verbose = FALSE)
# dt.combined.sct <- RunUMAP(dt.combined.sct, dims = 1:50)
# 
# pdf("output/plots/rpca_integrated_datasets.pdf", width=11, height = 8.5)
# DimPlot(dt.combined.sct, group.by = "orig.ident")
# dev.off()

print("dataset is integrated")
timestamp()
print(gc())

saveRDS(dt.combined.sct, "output/seurat_integrated.rds")

print("dataset is saved")
print(gc())
timestamp()

print("job complete... datasets are integrated")