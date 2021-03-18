# This script is part of the Seurat.R script.
# This will show you how to add automated labels and custom labels to your Seurat object.
# Labels generated using HPCA, MICD, and manual, with all output results to the CellxGene Browser.

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
##### Annotate Clusters using automated naming with Human Primary Cell Atlas (HPCA)
#SingleR Annotations
SCE.scObj <- as.SingleCellExperiment(scObj)

### HumanPrimaryCellAtlasData
hpca.se <- HumanPrimaryCellAtlasData()
common_hpca.scObj <- intersect(rownames(SCE.scObj), rownames(hpca.se))

hpca.sescObj <- hpca.se[common_hpca.scObj,]
hpca.scObj <- SCE.scObj[common_hpca.scObj,]
hpca.scObj <- logNormCounts(hpca.scObj)

# Main Label Groups based on 0.1 Resolution Clusterng
pred.hpca.scObj <- SingleR(test = hpca.scObj,
                           ref = hpca.sescObj,
                           labels = hpca.sescObj$label.main,
                           method = "cluster",
                           clusters = hpca.scObj$integrated_snn_res.0.1,
                           BPPARAM = MulticoreParam(workers=20))
scObj[["HPCA_Main_Labels"]] <- pred.hpca.scObj$labels[scObj$integrated_snn_res.0.1]

# Fine Label Groups based on 0.1 Resolution Clustering
pred.hpca.scObj <- SingleR(test = hpca.scObj,
                           ref = hpca.sescObj,
                           labels = hpca.sescObj$label.fine,
                           method = "cluster",
                           clusters = hpca.scObj$integrated_snn_res.0.1,
                           BPPARAM = MulticoreParam(workers=20))
scObj[["HPCA_Fine_Labels"]] <- pred.hpca.scObj$labels[scObj$integrated_snn_res.0.1]


##############################
##### Annotate Clusters using automated naming with Monaco Immune Data Set (MIDS)
#MonacoImmuneDataSet
micd.se <- MonacoImmuneData()
common_micd.scObj <- intersect(rownames(SCE.scObj), rownames(micd.se))

micd.sescObj <- micd.se[common_micd.scObj,]
micd.scObj <- SCE.scObj[common_micd.scObj,]
micd.scObj <- logNormCounts(micd.scObj)

# Main Label Groups based on 1.0 Resolution Clusterng
pred.micd.scObj <- SingleR(test = micd.scObj,
                           ref = micd.sescObj,
                           labels = micd.sescObj$label.main,
                           method = "cluster",
                           clusters = micd.scObj$integrated_snn_res.1,
                           BPPARAM = MulticoreParam(workers=20))
scObj[["MICD_Main_Labels"]] <- pred.micd.scObj$labels[scObj$integrated_snn_res.1]


# Fine Label Groups based on 1.0 Resolution Clusterng
pred.micd.scObj <- SingleR(test = micd.scObj,
                           ref = micd.sescObj,
                           labels = micd.sescObj$label.fine,
                           method = "cluster",
                           clusters = micd.scObj$integrated_snn_res.1,
                           BPPARAM = MulticoreParam(workers=20))
scObj[["MICD_Fine_Labels"]] <- pred.micd.scObj$labels[scObj$integrated_snn_res.1]

saveRDS(scObj, file="scObj-annotated.rds")
scObj <- readRDS("/project/shared/baker_cantarel/prune/analysis/single_cell/scObj/scObj-annotated.rds")
