rm(list = ls())
gc()
library(Seurat)
library(dplyr)
library(CellChat)
library(MASS)
library(tibble)
library(ggplot2)
library(stringr)
options(stringsAsFactors = FALSE)

setwd("/home/jiawen/myMLnet")
source("./codes/MLNconstruct.R")
source("./codes/Numericals.R")
source("./codes/Utilities.R")
source("./codes/Crosstalk_analysis.R")

################ Load the dataset ################
### CID4465 dataset
SeuratObj <- readRDS("/home/jiawen/myMLnet/benchmark/ESICCC/input/CID4465_ser.rds")
table(SeuratObj$celltype)
DefaultAssay(SeuratObj) <- 'SCT'
allgenes <- rownames(SeuratObj)

# SeuratObj <- subset(SeuratObj, celltype %in% c("CAFs","CancerEpithelial","Myeloid","PVL","T"))
# SeuratObj[["percent.mt"]] <- PercentageFeatureSet(SeuratObj, pattern = "^MT")
# VlnPlot(SeuratObj, features = c("nCount_Spatial", "nFeature_Spatial","percent.mt"))
# SeuratObj <- subset(SeuratObj, subset = nCount_Spatial < 40000 & nFeature_Spatial < 8000 & percent.mt < 7.5)
# SeuratObj <- SeuratObj %>%
#   NormalizeData() %>%
#  FindVariableFeatures() %>%
#  ScaleData(features = rownames(SeuratObj),vars.to.regress = "percent.mt")
# cell_anno <- data.frame(cell = colnames(SeuratObj), cluster = SeuratObj$celltype)

### HNSCC dataset
SeuratObj <- readRDS("/home/jiawen/myMLnet/datasets/nichenet/seurat.rds")
DefaultAssay(SeuratObj) <- 'RNA'
allgenes <- rownames(SeuratObj)
cell_anno <- data.frame(cell = colnames(SeuratObj), cluster = SeuratObj$celltype)

################ Load the prior knowledge of gene-gene interactions ################
LRDB <- readRDS("/home/jiawen/myMLnet/pathways/Nichenet/LR_human.rds") %>% distinct(from, to, .keep_all = T) %>% Filter_DB(allgenes)
RecTFDB <- readRDS("/home/jiawen/myMLnet/pathways/Nichenet/RTF_human.rds") %>% distinct(from, to, .keep_all = T) %>% Filter_DB(allgenes)
# RecTFDB <- readRDS("/home/jiawen/myMLnet/pathways/KEGG/RTF_human.rds")
TFTGDB <- readRDS("/home/jiawen/myMLnet/pathways/stMLnet/TFT_human.rds") %>% distinct(from, to, .keep_all = T) %>% Filter_DB(allgenes)
table(RecTFDB$database)

### Specify the target type
# target_type <- "malignant"
target_type <- "CancerEpithelial"

# Analyze the CCC
LR_original <- Infer_CCI(SeuratObj, cell_anno, LRDB = LRDB, cellchat_output = T, use_spatial = F, db_use = "human")

# Find DEG
TG_used <- FindMarkers(SeuratObj, target_type, min.pct = 0.25, only.pos = T, logfc.threshold = 0.25)
TG_used <- filter(TG_used, p_val_adj<1e-3) %>% rownames()

# Find the TF that regulates DEGs
input_dir <- "/home/jiawen/myMLnet/pythoncodes/inputs"

Prepare_Input_New(SeuratObj, target_type, TGs = TG_used, CCC_results = LR_original, RecTFDB, TFTGDB, data_dir = input_dir,
    assay = "RNA", datatype = "scale.data", imputation = F,exp_threshold = 0.1,
    CCC_threshold = 0.05, Fisher_threshold = 1)

TFTG_filtered$TF %>% table()
RTF_filtered$TF %>% table()
nrow(Exp_clu)
gc()
# Run HGNN
args.project <- "benchmark"
args.target <- target_type
conda_python <- "/home/jiawen/anaconda3/envs/SigXTalk/bin/python"

system2(conda_python, args = c("/home/jiawen/myMLnet/pythoncodes/main_new.py", paste("--project",shQuote(args.project)), paste("--target_type",args.target)))

output_dir <- '/home/jiawen/myMLnet/results/'
filen <- paste0(output_dir, 'benchmark/pathways_',target_type,'.csv')
RTFTG_results <- read.csv(filen, header = T)
RTFTG_results<- RTFTG_results[RTFTG_results$pred_label > 0.75, ]
RTFTG_results <- RTFTG_results[,1:3]
RTFTG_results$TF %>% table()
ress <- Rec_To_TFTG(Exp_clu, RTFTG_results, method = "rf", cutoff = 0.05, use_tidy = T)
ress$To %>% unique() %>% length()

aggregate(RTFTG_all$coexp, list(RTFTG_all$TF), max) 
rowSums(Exp_clu[RTFTG_results$TF %>% unique(),]>0)