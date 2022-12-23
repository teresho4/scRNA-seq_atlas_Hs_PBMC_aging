library(stringr)
library(Seurat)
library(dplyr)
library(tidyr)
dir.create('output', showWarnings = FALSE)

PROJECT=unlist(str_split(readLines('../project_names'), ','))
PROJECT <- PROJECT[which(PROJECT !='')]

POOL=unlist(str_split(readLines('./sample_names'), ','))
POOL <- POOL[which(POOL !='')]

seurat_obj1 <- Read10X(paste0("/storage1/fs1/martyomov/Active/NextGenImmunology/BSL3/", PROJECT, "/align/GEX+FB/", POOL, "/outs/filtered_feature_bc_matrix/", sep=""))
pbmc.hashtag <- CreateSeuratObject(counts = seurat_obj1$`Antibody Capture`[21:26,],  assay = "HTO")
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")
pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.99)
meta <- pbmc.hashtag@meta.data
patients <- unique(meta$HTO_classification)[!grepl("_",unique(meta$HTO_classification))]
patients <- patients[!grepl("Negative", patients)]
for(pat in patients){
    names <- meta %>% filter(HTO_classification == pat) %>% row.names()
    dir.create(paste('output/', pat, sep=''), showWarnings = FALSE)
    write.table(names, paste('/storage1/fs1/martyomov/Active/NextGenImmunology/BSL3/', PROJECT,'/sorc/snakemake/', POOL,'/output/', pat, '/', 'subset_barcodes.tsv',sep=''), 
                sep='\t', quote = F, row.names = F, col.names = F)
  }


