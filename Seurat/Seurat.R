# This script is a generic walkthrough of how to use and run Seurat.
# This accounts for multiple samples being integrated.
# It also shows how to add automatic annotations with SingleR 
# and output results to the CellxGene Browser.

# Input requires count matrix/tables(usually output from Cellranger Count from 10x Genomics)
# This script runs the SCTransform Vignette from Seurat

# $ module load python/3.6.4-anaconda
# $ module load R/4.0.2-gccmkl
# $ module add rstudio-sktop/1.1.456
# $ rstudio

##############################

#Set MiB limit needed for SCT Integration Step
Mb = 245000
options(future.globals.maxSize= Mb*1024^2)
options(java.parameters = "-Xmx64g")

library(Seurat)
library(Matrix)
library(ggplot2)
library(gridExtra)
library(viridis)
library(SingleR)
library(genefu)
library(tidyverse)
library(scater)
library(clustree)
library(celldex)
library(BiocParallel)
library(xlsx)
library(autothresholdr)
library(pheatmap)
library(RColorBrewer)

#Set workign directory if desired
setwd('')

##############################
##### Read in data

#Read 10x Cellranger Count output
sample1.data <- Read10X(data.dir="/path/to/directory/sample1/outs/filtered_feature_bc_matrix")
sample1 <- CreateSeuratObject(counts = sample1.data, project = "10x")

#Read Count Table
sample2.data <- read.csv("/path/to/directory/sample2/sample2.counts_data.txt", sep = ",", header = TRUE, row.names = 1)
sample2 <- CreateSeuratObject(counts = sample2.data, project = "CountTable")

#Read Pre-Generated Seurat Object
sample3 <- readRDS("/path/to/directory/sample3/sample3.rds")

#Merge Seurat Objects
Project <- merge(sample1, y = c(sample2, sample3),
                 add.cell.ids = c("sample1", "sample2", "sample3"),
                 project = "sc-Seurat")

#Save raw data object
saveRDS(Project, file = "/path/to/directory/project-raw.rds")



##############################
##############################
##############################



##############################
##### Filter and Integrate samples together using SCTransform. If you have a large dataset, you may reach a memory error. Please use the Seurat_XL.R script that implements the rPCA method of integration.

scObj <- readRDS("/path/to/directory/project-raw.rds")

#Add percent Mito and Run Filtering if needed. If using mouse data, percent mito may need to change pattern to "^mt-"
Idents(object = scObj) <- "orig.ident"
scObj[["percent.mt"]] <- PercentageFeatureSet(scObj, pattern = "^MT-")
VlnPlot(scObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

##### Normalize data
# SCT Integration
scObj.list <- SplitObject(scObj, split.by = "SampleID")
for (i in 1:length(scObj.list)) {
  scObj.list[[i]] <- SCTransform(scObj.list[[i]], method = "glmGamPoi", verbose = FALSE)
}
features <- SelectIntegrationFeatures(object.list = scObj.list, nfeatures = 3000)
scObj.list <- PrepSCTIntegration(object.list = scObj.list, anchor.features = features, verbose = FALSE)
scObj.anchors <- FindIntegrationAnchors(object.list = scObj.list, normalization.method = "SCT", anchor.features = features, verbose = FALSE)
scObj <- IntegrateData(anchorset = scObj.anchors, normalization.method = "SCT", verbose = FALSE)

##### Dimension Reduction
# TSNE and UMAP, add PCA
scObj <- RunPCA(scObj, verbose = FALSE, assay = "integrated")
scObj <- RunTSNE(scObj, dims = 1:30, verbose = FALSE, assay = "integrated")
scObj <- RunUMAP(scObj, dims = 1:30, verbose = FALSE, assay = "integrated")

#Run clustering with resolution 0.1
res <- c(0.1, 1.0)
scObj <- FindNeighbors(object = scObj, dims = 1:30, assay = "integrated")
scObj <- FindClusters(object = scObj, resolution = res, assay = "integrated")

### Add cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
scObj <- CellCycleScoring(scObj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE, assay = "SCT")

#Save filtered data object
saveRDS(scObj, file = "/path/to/directory/project-filtered.rds")



##############################
##############################
##############################



##############################
##### Annotate Clusters using automated naming with Human Primary Cell Atlas (HPCA)
#SingleR Annotations
SCE.scObj <- as.SingleCellExperiment(scObj)

### HumanPrimaryCellAtlasData
hpca.se <- HumanPrimaryCellAtlasData()
common_hpca.scObj <- intersect(rownames(SCE.scObj), rownames(hpca.se))

hpca.sescObj <- hpca.se[common_hpca.scObj,]
hpca.scObj <- SCE.scObj[common_hpca.scObj,]
hpca.scObj <- logNormCounts(hpca.scObj)

pred.hpca.scObj <- SingleR(test = hpca.scObj, ref = hpca.sescObj, labels = hpca.sescObj$label.main,
                               method = "cluster", clusters = hpca.scObj$integrated_snn_res.0.1, BPPARAM = MulticoreParam(workers=20))
# Add labels to seurat object
scObj[["HPCA_Labels"]] <- pred.hpca.scObj$labels[scObj$integrated_snn_res.0.1]

saveRDS(scObj, file="scObj-annotated.rds")
scObj <- readRDS("/project/shared/baker_cantarel/prune/analysis/single_cell/scObj/scObj-annotated.rds")
