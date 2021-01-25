# This is a tutorial on how to convert Seurat objects to CellPhoneDB ready files on BioHPC.
# This asusmes you have already run Seurat processing and have a usable Seurat object.
# CellPhoneDB is only able to run on human data currently. If your data is any organism other than human then please see the cpdb_mouse.R script on how to prepare your data.

library(Seurat)
library(tidyverse)
library(biomaRt)

#Set working drectory if desired
setwd("")

#Import Seurat object
sample <- readRDS("/path/to/directory/sample.rds")

#Generate counts data
sample_counts <- as.data.frame(sample@assays$RNA@data)
sample_counts$Symbol <- rownames(sample_counts)

#Convert Gene Symbol to Ensembl Gene ID
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
symbol2ID <- getBM(mart = human,
                   attributes = c("hgnc_symbol", "ensembl_gene_id"),
                   filters = "hgnc_symbol",
                   values = rownames(sample_counts),
                   uniqueRows = TRUE)
symbol2ID <- symbol2ID %>% rename(Gene = ensembl_gene_id)
sample_counts <- merge(x = sample_counts,
                       y = symbol2ID,
                       by.x = "Symbol",
                       by.y = "hgnc_symbol",
                       all.x = TRUE)
sample_counts <- sample_counts[,!(names(sample_counts) %in% c("Symbol"))] %>%
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
