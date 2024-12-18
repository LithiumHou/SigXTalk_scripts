rm(list = ls())
gc()
library(Seurat)
library(dplyr)
# library(CellChat)
library(MASS)
library(tibble)
library(ggplot2)
library(patchwork)
setwd("/home/jiawen/myMLnet")
source("./codes/MLNconstruct.R")
source("./codes/Numericals.R")
source("./codes/Utilities.R")
source("./codes/Crosstalk_analysis.R")



##### Data preprocessing for HNSCC datasets
hnscc_expression = readRDS("./datasets/hnscc/hnscc_expression.rds")
expression = hnscc_expression$expression
sample_info = hnscc_expression$sample_info
tumors_remove = c("HN10","HN","HN12", "HN13", "HN24", "HN7", "HN8","HN23")
# tumors_remove = c("HN")
sample_info <- sample_info %>% filter(`Lymph node` == 0 & !(tumor %in% tumors_remove))
expression <- expression[sample_info$cell,] %>% t()
cell_anno <- data.frame(cell = sample_info$cell, cluster = sample_info$`non-cancer cell type`)
cell_anno$cluster[cell_anno$cluster == '0'] <- "malignant"
cell_anno$cluster[cell_anno$cluster == '-Fibroblast'] <- "Fibroblast"
cell_anno$cluster[cell_anno$cluster == 'B cell'] <- "B_cell"
cell_anno$cluster[cell_anno$cluster == 'T cell'] <- "T_cell"
SeuratObj <- CreateSeuratObject(expression)

SeuratObj[["percent.mt"]] <- PercentageFeatureSet(SeuratObj, pattern = "^MT")
rm(expression,hnscc_expression,sample_info)
Idents(SeuratObj) <- cell_anno$cluster
gc()
VlnPlot(SeuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

SeuratObj <- subset(SeuratObj, nCount_RNA < 25000 & nFeature_RNA<10000 & percent.mt<2.0)
SeuratObj <- SeuratObj %>% NormalizeData() %>% FindVariableFeatures() 
SeuratObj <- ScaleData(SeuratObj, features = rownames(SeuratObj), vars.to.regress = "percent.mt")
SeuratObj <- SCTransform(SeuratObj,vst.flavor = "v1")
SeuratObj <- RunPCA(SeuratObj) %>% RunUMAP(dims = 1:10)

DimPlot(SeuratObj, reduction = "umap")
saveRDS(SeuratObj, "/home/jiawen/myMLnet/datasets/nichenet/seurat1006.rds")

##### COVID datasets

library(dplyr)
library(Seurat)
library(ggplot2)

setwd("/home/jiawen/myMLnet")
counts <- Matrix::readMM("./datasets/COVID/expression_data.mtx")
genes <- read.csv("./datasets/COVID/lung_geneNames_upload.csv",header = F)
cells <- read.csv("./datasets/COVID/lung_cellNames.csv",header = F)
rownames(counts) <- genes$V1; colnames(counts) <- cells$V1;
cells$V1 %>% unique() %>% length()
cell_anno_all <- read.csv("./datasets/COVID/lung_metaData.txt",sep = "\t")
cell_anno_all <- cell_anno_all[-1,]
counts <- counts[,cell_anno_all$NAME]

SeuratObj_all <- CreateSeuratObject(counts,meta.data = cell_anno_all)
SeuratObj_all[["percent.mt"]] <- PercentageFeatureSet(SeuratObj_all, pattern = "^MT")

SeuratObj_covid <- subset(SeuratObj_all, donor_id == "L17cov")
VlnPlot(SeuratObj_covid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
SeuratObj_covid <- subset(SeuratObj_covid, nCount_RNA < 4000 & nFeature_RNA< 4000 & nFeature_RNA > 200 & percent.mt<4.0)

SeuratObj_covid <- SeuratObj_covid %>% NormalizeData() %>% FindVariableFeatures() 
SeuratObj_covid <- ScaleData(SeuratObj_covid,features = rownames(SeuratObj_covid), vars.to.regress = "percent.mt") %>% SCTransform()
SeuratObj_covid <- RunPCA(SeuratObj_covid) %>% RunUMAP(dims = 1:10)
Idents(SeuratObj_covid) <- SeuratObj_covid$cell_type_main

DimPlot(SeuratObj_covid, reduction = "umap")
saveRDS(SeuratObj_covid, "./datasets/COVID/seurat_covid.rds")

# The Healthy control
SeuratObj_hc <- subset(SeuratObj_all, donor_id == "C55ctr")
VlnPlot(SeuratObj_hc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
SeuratObj_hc <- subset(SeuratObj_hc, nCount_RNA < 4000 & nFeature_RNA< 3000 & nFeature_RNA > 200 & percent.mt<4.0)

SeuratObj_hc <- SeuratObj_hc %>% NormalizeData() %>% FindVariableFeatures() 
SeuratObj_hc <- ScaleData(SeuratObj_hc,features = rownames(SeuratObj_hc), vars.to.regress = "percent.mt") %>% SCTransform()
SeuratObj_hc <- RunPCA(SeuratObj_hc) %>% RunUMAP(dims = 1:10)
Idents(SeuratObj_hc) <- SeuratObj_hc$cell_type_main

DimPlot(SeuratObj_hc, reduction = "umap")
saveRDS(SeuratObj_hc, "./datasets/COVID/seurat_hc.rds")

##### The mouse lung dataset
rna_counts <- Read10X("./datasets/GSE141259/HighResolution", gene.column = 1)
cell_info <- read.csv("./datasets/GSE141259/HighResolution/cellinfo2.csv")
SeuratObj_all <- CreateSeuratObject(rna_counts,meta.data = cell_info)
table(cell_info$meta_celltype)
table(cell_info$sample_id)
Idents(SeuratObj_all) <- cell_info$meta_celltype
SeuratObj_all[["percent.mt"]] <- PercentageFeatureSet(SeuratObj_all, pattern = "^mt-")
gene_info <- rownames(rna_counts) %>% as.data.frame()
VlnPlot(SeuratObj_all, features = c("nCount_RNA","nFeature_RNA","percent.mt"))
SeuratObj_all <- subset(
  x = SeuratObj_all,
  subset = nCount_RNA < 10000 &
    nCount_RNA < 3000 &
    percent.mt < 7.5
)
all_timepoints <- list(
  c("day 2", "day 3"),
  c("day 10", "day 11"),
  c("day 15")
)
cell_info$time_point %>% table()

rm(SeuratObj)

for (i in 1:length(all_timepoints)) {

  SeuratObj <- subset(SeuratObj_all, time_point %in% all_timepoints[[i]])
  Idents(SeuratObj) <- SeuratObj$meta_celltype
  SeuratObj <- SeuratObj %>% NormalizeData() %>% FindVariableFeatures() 
  SeuratObj <- ScaleData(SeuratObj, features = rownames(SeuratObj), vars.to.regress = "percent.mt")
  SeuratObj <- SCTransform(SeuratObj) %>%  RunPCA() %>% RunUMAP(dims = 1:10)
  cell_anno <- data.frame(cell = SeuratObj$cell_barcode, cluster = SeuratObj$meta_celltype)
  all_types <- aggregate(cell_anno$cell, list(cell_anno$cluster), length)
  types_used <- (all_types %>% filter(x > 80))$Group.1
  SeuratObj <- subset(SeuratObj, meta_celltype %in% types_used)

  filen <- paste0("./datasets/GSE141259/Seurat_t",i,".rds")
  saveRDS(SeuratObj, file = filen)
}

##### The spatial dataset
SeuratObj <- readRDS("./datasets/visium/cortex.rds")
cell_anno <- Idents(SeuratObj) %>%
    as.data.frame() %>%
    rownames_to_column()
DefaultAssay(SeuratObj) <- "Spatial"
SeuratObj[["percent.mt"]] <- PercentageFeatureSet(SeuratObj, pattern = "^mt-")
VlnPlot(SeuratObj, features = c("nCount_Spatial","nFeature_Spatial","percent.mt"))
SeuratObj <- SeuratObj %>% NormalizeData() %>% FindVariableFeatures() 
SeuratObj <- ScaleData(SeuratObj, features = rownames(SeuratObj), vars.to.regress = "percent.mt")
SeuratObj <- SCTransform(SeuratObj,assay = "Spatial") %>%  RunPCA() %>% RunUMAP(dims = 1:10)
SpatialDimPlot(SeuratObj, label = TRUE, pt.size.factor = 50, label.size = 6,image.alpha = 0.5)
ncol(SeuratObj@assays$Spatial$counts)
saveRDS(SeuratObj, file = "./datasets/visium/cortex_new.rds")
gc()

