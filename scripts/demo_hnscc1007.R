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

SeuratObj <- readRDS("/home/jiawen/myMLnet/datasets/nichenet/seurat1006.rds")
cell_anno <- data.frame(cell = names(Idents(SeuratObj)), cluster = Idents(SeuratObj) %>% as.character())
allgenes <-  rownames(SeuratObj@assays$RNA$data)
# load("/home/jiawen/myMLnet/pathways/LR_layer1_human.rda")
LRDB <- readRDS("/home/jiawen/myMLnet/pathways/Nichenet/LR_human.rds") %>% distinct(from, to, .keep_all = T) %>% Filter_DB(allgenes)
RecTFDB <- readRDS("/home/jiawen/myMLnet/pathways/Nichenet/RTF_human.rds") %>% distinct(from, to, .keep_all = T) %>% Filter_DB(allgenes)
# RecTFDB <- readRDS("/home/jiawen/myMLnet/pathways/KEGG/RTF_human.rds")
TFTGDB <- readRDS("/home/jiawen/myMLnet/pathways/Nichenet/TFT_human.rds") %>% distinct(from, to, .keep_all = T) %>% Filter_DB(allgenes)

LR_original <- Infer_CCI(SeuratObj, cell_anno, LRDB = LRDB, cellchat_output = T, use_spatial = F, db_use = "human")

table(SeuratObj$cluster)
args.project <- "hnscc"
types_used <- c("CAF","Endothelial","myofibroblast","T_cell","malignant")
# types_used <- c("malignant")

output_dir <- '/home/jiawen/myMLnet/results/'
allresult <- list()
for(i in 1:length(types_used)){

    target_type <- types_used[i]
    message(paste0("Analysing the cell type ",target_type,"...\n"))

    TG_used <- FindMarkers(SeuratObj, target_type, min.diff.pct = 0.25, only.pos = T, logfc.threshold = 0.25)
    TG_used <- filter(TG_used, p_val_adj<1e-3) %>% rownames()

    input_dir <- "/home/jiawen/myMLnet/pythoncodes/inputs"
    Prepare_Input_New(SeuratObj, target_type, TGs = TG_used, CCC_results = LR_original, RecTFDB, TFTGDB, data_dir = input_dir,
        assay = "RNA", datatype = "scale.data", imputation = F,exp_threshold = 0.1,
        CCC_threshold = 0.1, Fisher_threshold = .99)

    gc()
    args.target <- target_type
    conda_python <- "/home/jiawen/anaconda3/envs/SigXTalk/bin/python"

    system2(conda_python, args = c("/home/jiawen/myMLnet/pythoncodes/main_new.py", paste("--project",shQuote(args.project)), paste("--target_type",args.target)))

    filen <- paste0(output_dir, args.project,'/pathways_',target_type,'.csv')
    RTFTG_results <- read.csv(filen, header = T)
    RTFTG_results <- RTFTG_results[RTFTG_results$pred_label > 0.75, ]
    RTFTG_results <- RTFTG_results[,1:3]
    RTFTG_results$Receptor %>% unique() %>% length()
    gc()
    Exp_clu <- Get_Exp_Clu(SeuratObj, clusterID = target_type, assay = "RNA", datatype = "data", cutoff = 0.1)
    ress <- Rec_To_TFTG(Exp_clu, RTFTG_results, method = "rf", cutoff = 0.1, use_tidy = T)
    filen <- paste0(output_dir, args.project,'/results_',target_type,'.csv')
    write.table(ress,file = filen, quote = F, sep = ",")
    allresult[[target_type]] <- ress
    gc()
}


filen <- paste0(output_dir, args.project,'/results.RDS')
saveRDS(allresult, file)

