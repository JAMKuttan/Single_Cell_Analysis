# This is a tutorial on how to convert Seurat objects to CellPhoneDB ready files on BioHPC.
# This asusmes you have already run Seurat processing and have a usable Seurat object.
# CellPhoneDB is only able to run on human data currently.
# This script is for non-human data samples. Here mouse is used as an example. CellPhoneDB suggests usng the human homologs for analysis. There are some drawbacks but the analysis has been shown to be usable.
# This script generates homologs using biomaRt but other methods are available.

library(Seurat)
library(tidyverse)
library(biomaRt)

#Set working drectory if desired
setwd("")

#Convert Gene Symbol to Ensembl Gene ID
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

convertMouseGeneList <- function(x){
  require(biomaRt)

  genesV2 = getLDS(attributes = c("mgi_symbol"),
               filters = "mgi_symbol", values = x, mart = mouse,
               attributesL = c("ensembl_gene_id"), martL = human,
               uniqueRows = TRUE)
  genesV2$MGI.symbol <- toupper(genesV2$MGI.symbol)
  genesV2 <- genesV2 %>% rename(Gene = Gene.stable.ID)
  return(genesV2)
}

#Import Seurat object
sample <- readRDS("/path/to/directory/sample.rds")

#Generate counts data
sample_counts <- as.data.frame(sample@assays$RNA@data)
sample_counts$Mouse_Genes <- toupper(rownames(sample_counts))
mm2hu <- convertMouseGeneList(rownames(sample_counts))
sample_counts <- merge(x = sample_counts,
                       y = mm2hu,
                       by.x = "Mouse_Genes",
                       by.y = "MGI.symbol",
                       all.x = TRUE)
sample_counts <- sample_counts[which(!is.na(sample_counts$Gene)),]
sample_counts <- sample_counts[,!(names(sample_counts) %in% c("Mouse_Genes"))] %>%
  dplyr::select(Gene, everything())

#Write counts.txt
write.table(sample_counts, file = "/path/to/directory/sample_counts.txt",
            quote = FALSE,
            sep = '\t',
            row.names = FALSE,
            col.names = TRUE)

#Generate meta data (this is to assign whatever identidies you require, here the cluster labels from resolution 1.0 are used)
sample_meta <- data.frame(Cell = rownames(sample@meta.data),
                          cell_type = sample@meta.data$RNA_snn_res.1)

#Write meta.txt
write.table(sample_meta, file = "/path/to/directory/sample_meta.txt",
            quote = FALSE,
            sep = '\t',
            row.names = FALSE,
            col.names = TRUE)
