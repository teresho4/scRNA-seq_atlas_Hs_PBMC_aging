library(Seurat)
library(harmony)

### CD4+ T cells ###

norm_counts <- readRDS('datasets_submission/cd4_t_cells/cd4_rna.rds')
meta <- read.csv('datasets_submission/cd4_t_cells/cd4_metadata.csv', row.names = 'X')
norm_counts_hto <-  readRDS('datasets_submission/cd4_t_cells/cd4_hto.rds')
data <- CreateSeuratObject(counts = norm_counts, meta.data = meta)
data@assays$RNA@data <- norm_counts
data[['HTO']] <-  CreateAssayObject(norm_counts_hto)
data@assays$HTO@data <- norm_counts_hto
data <- FindVariableFeatures(object = data, selection.method = 'mean.var.plot', mean.cutoff = c(0.01, Inf), dispersion.cutoff = c(0.55, Inf))
data@assays$RNA@var.features <- data@assays$RNA@var.features[!grepl("^TRA|^TRB|^IGH|^IGK|^IGL", data@assays$RNA@var.features)]
data <- ScaleData(object = data, features = VariableFeatures(object = data), vars.to.regress = c("nCount_RNA", "percent.mt"))
data <- RunPCA(object = data)
data <- RunHarmony(object = data, group.by.vars = c("Batch"), assay.use = "RNA", max.iter.harmony = 20)
data <- RunUMAP(data, dims = 1:30, reduction = "harmony")
data <- FindNeighbors(data, dims = 1:30, reduction = "harmony", k.param = 15)
data <- FindClusters(data, resolution = 0.5)

### CD4+ helper memory T cells ###

norm_counts <- readRDS('gex_hto_data/cd4_t_helper_memory_cells/cd4_helper_memory_rna.rds')
meta <- read.csv('gex_hto_data/cd4_t_helper_memory_cells/cd4_helper_memory_metadata.csv', row.names = 'X')
norm_counts_hto <-  readRDS('gex_hto_data/cd4_t_helper_memory_cells/cd4_helper_memory_hto.rds')
data <- CreateSeuratObject(counts = norm_counts, meta.data = meta)
data@assays$RNA@data <- norm_counts
data[['HTO']] <-  CreateAssayObject(norm_counts_hto)
data@assays$HTO@data <- norm_counts_hto
data <- FindVariableFeatures(object = data, selection.method = 'mean.var.plot', mean.cutoff = c(0.01, Inf), dispersion.cutoff = c(0.62, Inf))
data@assays$RNA@var.features <- data@assays$RNA@var.features[!grepl("^TRA|^TRB|^IGH|^IGK|^IGL", data@assays$RNA@var.features)]
data <- ScaleData(object = data, features = VariableFeatures(object = data), vars.to.regress = c("nCount_RNA", "percent.mt"))
data <- RunPCA(object = data)
data <- RunHarmony(object = data, group.by.vars = c("Batch"), assay.use = "RNA", max.iter.harmony = 20)
data <- RunUMAP(data, dims = 1:20, reduction = "harmony")
data <- FindNeighbors(data, dims = 1:20, reduction = "harmony", k.param = 15)
data <- FindClusters(data, resolution = 0.3)

### Conventional CD8+ T cells ###

norm_counts <- readRDS('gex_hto_data/conventional_cd8_t_cells/conventional_cd8_rna.rds')
meta <- read.csv('gex_hto_data/conventional_cd8_t_cells/conventional_cd8_metadata.csv', row.names = 'X')
norm_counts_hto <-  readRDS('gex_hto_data/conventional_cd8_t_cells/conventional_cd8_hto.rds')
data <- CreateSeuratObject(counts = norm_counts, meta.data = meta)
data@assays$RNA@data <- norm_counts
data[['HTO']] <-  CreateAssayObject(norm_counts_hto)
data@assays$HTO@data <- norm_counts_hto
data <- FindVariableFeatures(object = data, selection.method = 'mean.var.plot', mean.cutoff = c(0.01, Inf), dispersion.cutoff = c(0.65, Inf))
data@assays$RNA@var.features <- data@assays$RNA@var.features[!grepl("^TRA|^TRB|^IGH|^IGK|^IGL", data@assays$RNA@var.features)]
data <- ScaleData(object = data, features = VariableFeatures(object = data), vars.to.regress = c("nCount_RNA", "percent.mt"))
data <- RunPCA(object = data)
data <- RunHarmony(object = data, group.by.vars = c("Batch"), assay.use = "RNA", max.iter.harmony = 20)
data <- RunUMAP(data, dims = 1:30, reduction = "harmony")
data <- FindNeighbors(data, dims = 1:30, reduction = "harmony", k.param = 15)
data <- FindClusters(data, resolution = 0.8)

###  Myeloid cells ###

norm_counts <- readRDS('gex_hto_data/myeloid_cells/myeloid_cells_rna.rds')
meta <- read.csv('gex_hto_data/myeloid_cells/myeloid_cells_metadata.csv', row.names = 'X')
norm_counts_hto <-  readRDS('gex_hto_data/myeloid_cells/myeloid_cells_hto.rds')
data <- CreateSeuratObject(counts = norm_counts, meta.data = meta)
data@assays$RNA@data <- norm_counts
data[['HTO']] <-  CreateAssayObject(norm_counts_hto)
data@assays$HTO@data <- norm_counts_hto
data <- FindVariableFeatures(object = data, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, Inf), dispersion.cutoff = c(0.65, Inf))
data@assays$RNA@var.features <- data@assays$RNA@var.features[!grepl("^TRA|^TRB|^IGH|^IGK|^IGL", data@assays$RNA@var.features)]
data <- ScaleData(object = data, features = VariableFeatures(object = data), vars.to.regress = c("nCount_RNA", "percent.mt"))
data <- RunPCA(object = data)
data <- RunHarmony(object = data, group.by.vars = c("Batch"), assay.use = "RNA", max.iter.harmony = 20)
data <- RunUMAP(data, dims = 1:10, reduction = "harmony")
data <- FindNeighbors(data, dims = 1:10, reduction = "harmony", k.param = 15)
data <- FindClusters(data, resolution = 0.1)

### NK cells ###

norm_counts <- readRDS('gex_hto_data/nk_cells/nk_cells_rna.rds')
meta <- read.csv('gex_hto_data/nk_cells/nk_cells_metadata.csv', row.names = 'X')
norm_counts_hto <-  readRDS('gex_hto_data/nk_cells/nk_cells_hto.rds')
data <- CreateSeuratObject(counts = norm_counts, meta.data = meta)
data@assays$RNA@data <- norm_counts
data[['HTO']] <-  CreateAssayObject(norm_counts_hto)
data@assays$HTO@data <- norm_counts_hto
data <- FindVariableFeatures(object = data, selection.method = 'mean.var.plot', mean.cutoff = c(0.01, Inf), dispersion.cutoff = c(0.69, Inf))
data@assays$RNA@var.features <- data@assays$RNA@var.features[!grepl("^TRA|^TRB|^IGH|^IGK|^IGL", data@assays$RNA@var.features)]
data <- ScaleData(object = data, features = VariableFeatures(object = data), vars.to.regress = c("nCount_RNA", "percent.mt"))
data <- RunPCA(object = data)
data <- RunHarmony(object = data, group.by.vars = c("Batch"), assay.use = "RNA", max.iter.harmony = 20)
data <- RunUMAP(data, dims = 1:15, reduction = "harmony")
data <- FindNeighbors(data, dims = 1:15, reduction = "harmony", k.param = 15)
data <- FindClusters(data, resolution = 0.4)

### gd T cells ###

norm_counts <- readRDS('gex_hto_data/gd_t_cells/gd_t_cells_rna.rds')
meta <- read.csv('gex_hto_data/gd_t_cells/gd_t_cells_metadata.csv', row.names = 'X')
norm_counts_hto <-  readRDS('gex_hto_data/gd_t_cells/gd_t_cells_hto.rds')
data <- CreateSeuratObject(counts = norm_counts, meta.data = meta)
data@assays$RNA@data <- norm_counts
data[['HTO']] <-  CreateAssayObject(norm_counts_hto)
data@assays$HTO@data <- norm_counts_hto
data <- FindVariableFeatures(object = data, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, Inf), dispersion.cutoff = c(0.79, Inf))
data@assays$RNA@var.features <- data@assays$RNA@var.features[!grepl("^TRA|^TRB|^IGH|^IGK|^IGL", data@assays$RNA@var.features)]
data <- ScaleData(object = data, features = VariableFeatures(object = data), vars.to.regress = c("nCount_RNA", "percent.mt"))
data <- RunPCA(object = data)
data <- RunHarmony(object = data, group.by.vars = c("Batch"), assay.use = "RNA", max.iter.harmony = 20)
data <- RunUMAP(data, dims = 1:20, reduction = "harmony")
data <- FindNeighbors(data, dims = 1:20, reduction = "harmony", k.param = 15)
data <- FindClusters(data, resolution = 0.3)

### B cells ###

norm_counts <- readRDS('gex_hto_data/b_cells/b_cells_rna.rds')
meta <- read.csv('gex_hto_data/b_cells/b_cells_metadata.csv', row.names = 'X')
norm_counts_hto <-  readRDS('gex_hto_data/b_cells/b_cells_hto.rds')
data <- CreateSeuratObject(counts = norm_counts, meta.data = meta)
data@assays$RNA@data <- norm_counts
data[['HTO']] <-  CreateAssayObject(norm_counts_hto)
data@assays$HTO@data <- norm_counts_hto
data <- FindVariableFeatures(object = data, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
data@assays$RNA@var.features <- data@assays$RNA@var.features[!grepl("^TRA|^TRB|^IGH|^IGK|^IGL", data@assays$RNA@var.features)]
data <- ScaleData(object = data, features = VariableFeatures(object = data), vars.to.regress = c("nCount_RNA", "percent.mt"))
data <- RunPCA(object = data)
data <- RunHarmony(object = data, group.by.vars = c("Batch"), assay.use = "RNA", max.iter.harmony = 20)
data <- RunUMAP(data, dims = 1:30, reduction = "harmony")
data <- FindNeighbors(data, dims = 1:30, reduction = "harmony", k.param = 15)
data <- FindClusters(data, resolution = 0.3)

### Progenitor cells ###

norm_counts <- readRDS('gex_hto_data/progenitor_cells/progenitor_cells_rna.rds')
meta <- read.csv('gex_hto_data/progenitor_cells/progenitor_cells_metadata.csv', row.names = 'X')
norm_counts_hto <-  readRDS('gex_hto_data/progenitor_cells/progenitor_cells_hto.rds')
data <- CreateSeuratObject(counts = norm_counts, meta.data = meta)
data@assays$RNA@data <- norm_counts
data[['HTO']] <-  CreateAssayObject(norm_counts_hto)
data@assays$HTO@data <- norm_counts_hto
data <- FindVariableFeatures(object = data, selection.method = 'mean.var.plot', mean.cutoff = c(0.01, Inf), dispersion.cutoff = c(1.6, Inf))
data@assays$RNA@var.features <- data@assays$RNA@var.features[!grepl("^TRA|^TRB|^IGH|^IGK|^IGL", data@assays$RNA@var.features)]
data <- ScaleData(object = data, features = VariableFeatures(object = data), vars.to.regress = c("nCount_RNA", "percent.mt"))
data <- RunPCA(object = data)
data <- RunHarmony(object = data, group.by.vars = c("Batch"), assay.use = "RNA", max.iter.harmony = 20)
data <- RunUMAP(data, dims = 1:10, reduction = "harmony")
data <- FindNeighbors(data, dims = 1:10, reduction = "harmony")
data <- FindClusters(data, resolution = 0.3)

### MAIT cells ###

norm_counts <- readRDS('gex_hto_data/mait_cells/mait_cells_rna.rds')
meta <- read.csv('gex_hto_data/mait_cells/mait_cells_metadata.csv', row.names = 'X')
norm_counts_hto <-  readRDS('gex_hto_data/mait_cells/mait_cells_hto.rds')
data <- CreateSeuratObject(counts = norm_counts, meta.data = meta)
data@assays$RNA@data <- norm_counts
data[['HTO']] <-  CreateAssayObject(norm_counts_hto)
data@assays$HTO@data <- norm_counts_hto
data <- FindVariableFeatures(object = data, selection.method = 'mean.var.plot', mean.cutoff = c(0.01, Inf), dispersion.cutoff = c(0.91, Inf))
data@assays$RNA@var.features <- data@assays$RNA@var.features[!grepl("^TRA|^TRB|^IGH|^IGK|^IGL", data@assays$RNA@var.features)]
data <- ScaleData(object = data, features = VariableFeatures(object = data), vars.to.regress = c("nCount_RNA", "percent.mt"))
data <- RunPCA(object = data)
data <- RunHarmony(object = data, group.by.vars = c("Batch"), assay.use = "RNA", max.iter.harmony = 20)
data <- RunUMAP(data, dims = 1:15, reduction = "harmony")
data <- FindNeighbors(data, dims = 1:15, reduction = "harmony")
data <- FindClusters(data, resolution = 0.3)
