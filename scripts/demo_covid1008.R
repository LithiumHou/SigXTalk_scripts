rm(list = ls())
gc()
library(Seurat)
library(dplyr)
library(CellChat)
library(MASS)
library(tibble)
library(ggplot2)
# library(hdf5r)
library(future)
options(stringsAsFactors = FALSE)

setwd("/home/jiawen/myMLnet")
source("./codes/MLNconstruct.R")
source("./codes/Numericals.R")
source("./codes/Utilities.R")
source("./codes/Crosstalk_analysis.R")

# Healthy control
SeuratObj <- readRDS("./datasets/COVID/seurat_hc.rds")
SeuratObj$cluster <- SeuratObj$cell_type_main
cell_anno <- SeuratObj$cell_type_main
allgenes <-  rownames(SeuratObj@assays$RNA$data)
# load("/home/jiawen/myMLnet/pathways/LR_layer1_human.rda")
LRDB <- readRDS("/home/jiawen/myMLnet/pathways/Nichenet/LR_human.rds") %>% distinct(from, to, .keep_all = T) %>% Filter_DB(allgenes)
RecTFDB <- readRDS("/home/jiawen/myMLnet/pathways/Nichenet/RTF_human.rds") %>% distinct(from, to, .keep_all = T) %>% Filter_DB(allgenes)
# RecTFDB <- readRDS("/home/jiawen/myMLnet/pathways/KEGG/RTF_human.rds")
TFTGDB <- readRDS("/home/jiawen/myMLnet/pathways/Nichenet/TFT_human.rds") %>% distinct(from, to, .keep_all = T) %>% Filter_DB(allgenes)

LR_original <- Infer_CCI(SeuratObj, cell_anno, LRDB = LRDB, cellchat_output = T, use_spatial = F, db_use = "human")

table(SeuratObj$cluster)
args.target <- 'hc'
args.project <- 'COVID'
target_type <- "Fibroblasts"
TG_used <- FindMarkers(SeuratObj, target_type, min.pct = 0.25, only.pos = T, logfc.threshold = 0.25)
TG_used <- filter(TG_used, p_val_adj<1e-3) %>% rownames()

input_dir <- "/home/jiawen/myMLnet/pythoncodes/inputs"
Prepare_Input_New(SeuratObj, target_type, TGs = TG_used, CCC_results = LR_original, RecTFDB, TFTGDB, data_dir = input_dir,
    assay = "RNA", datatype = "scale.data", imputation = F,exp_threshold = 0.1,
    CCC_threshold = 0.1, Fisher_threshold = 1)
gc()
# args.target <- target_type
conda_python <- "/home/jiawen/anaconda3/envs/SigXTalk/bin/python"
TFTG_filtered$TG %>%unique()
system2(conda_python, args = c("/home/jiawen/myMLnet/pythoncodes/main_new.py", paste("--project",shQuote(args.project)), paste("--target_type",args.target)))

output_dir <- '/home/jiawen/myMLnet/results/'
filen <- paste0(output_dir, args.project,'/pathways_',args.target,'.csv')
RTFTG_results <- read.csv(filen, header = T)
RTFTG_results <- RTFTG_results[RTFTG_results$pred_label > 0.75, ]
RTFTG_results <- RTFTG_results[,1:3]
RTFTG_results$TG %>% unique()
Exp_clu <- Get_Exp_Clu(SeuratObj, clusterID = target_type, assay = "RNA", datatype = "data", cutoff = 0.05)
ress <- Rec_To_TFTG(Exp_clu, RTFTG_results, method = "rf", cutoff = 0.1, use_tidy = T)
filen <- paste0(output_dir, args.project,'/results_',args.target,'.csv')

write.table(ress,file = filen, quote = F, sep = ",")
gc()


###############
# COVID patient
SeuratObj <- readRDS("./datasets/COVID/seurat_covid.rds")
SeuratObj$cluster <- SeuratObj$cell_type_main
cell_anno <- SeuratObj$cell_type_main
allgenes <-  rownames(SeuratObj@assays$RNA$data)
# load("/home/jiawen/myMLnet/pathways/LR_layer1_human.rda")
LRDB <- readRDS("/home/jiawen/myMLnet/pathways/Nichenet/LR_human.rds") %>% distinct(from, to, .keep_all = T) %>% Filter_DB(allgenes)
RecTFDB <- readRDS("/home/jiawen/myMLnet/pathways/Nichenet/RTF_human.rds") %>% distinct(from, to, .keep_all = T) %>% Filter_DB(allgenes)
# RecTFDB <- readRDS("/home/jiawen/myMLnet/pathways/KEGG/RTF_human.rds")
TFTGDB <- readRDS("/home/jiawen/myMLnet/pathways/Nichenet/TFT_human.rds") %>% distinct(from, to, .keep_all = T) %>% Filter_DB(allgenes)

LR_original <- Infer_CCI(SeuratObj, cell_anno, LRDB = LRDB, cellchat_output = T, use_spatial = F, db_use = "human")

table(SeuratObj$cluster)
args.target <- 'patient'
args.project <- 'COVID'
target_type <- "Fibroblasts"
TG_used <- FindMarkers(SeuratObj, target_type, min.pct = 0.25, only.pos = T, logfc.threshold = 0.25)
TG_used <- filter(TG_used, p_val<1e-2) %>% rownames()

input_dir <- "/home/jiawen/myMLnet/pythoncodes/inputs"
Prepare_Input_New(SeuratObj, target_type, TGs = TG_used, CCC_results = LR_original, RecTFDB, TFTGDB, data_dir = input_dir,
    assay = "RNA", datatype = "scale.data", imputation = F,exp_threshold = 0.1,
    CCC_threshold = 0.1, Fisher_threshold = 1)
gc()

c(TFTG_filtered$TF,TFTG_filtered$TG,RTF_filtered$Receptor) %>% unique()

conda_python <- "/home/jiawen/anaconda3/envs/SigXTalk/bin/python"
TFTG_filtered$TG %>% unique()
system2(conda_python, args = c("/home/jiawen/myMLnet/pythoncodes/main_new.py", paste("--project",shQuote(args.project)), paste("--target_type",args.target)))
output_dir <- '/home/jiawen/myMLnet/results/'
filen <- paste0(output_dir, args.project,'/pathways_',args.target,'.csv')
RTFTG_results <- read.csv(filen, header = T)
RTFTG_results <- RTFTG_results[RTFTG_results$pred_label > 0.75, ]
RTFTG_results <- RTFTG_results[,1:3]
RTFTG_results$TG %>% unique()
Exp_clu <- Get_Exp_Clu(SeuratObj, clusterID = target_type, assay = "RNA", datatype = "data", cutoff = 0.05)
ress <- Rec_To_TFTG(Exp_clu, RTFTG_results, method = "rf", cutoff = 0.1, use_tidy = T)
filen <- paste0(output_dir, args.project,'/results_',args.target,'.csv')
write.table(ress,file = filen, quote = F, sep = ",")


