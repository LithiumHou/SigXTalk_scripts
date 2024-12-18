rm(list = ls())
gc()
library(Seurat)
library(SeuratData)
library(dplyr)
library(tibble)
library(ggplot2)
library(foreach)
library(patchwork)
options(stringsAsFactors = FALSE)

setwd("/home/jiawen/myMLnet")
source("./codes/MLNconstruct.R")
source("./codes/Numericals.R")
source("./codes/Utilities.R")
source("./codes/Crosstalk_analysis.R")

all_timepoints <- list(
  c("day 3", "day 2"),
  c("day 10", "day 11"),
  c("day 15")
)

for (i in 1:length(all_timepoints)) {

    filen <- paste0("./datasets/GSE141259/Seurat_t",i,".rds")
    SeuratObj <- readRDS(filen)
    cell_anno <- Idents(SeuratObj) %>%
        as.data.frame() %>%
        rownames_to_column()
    SeuratObj$cluster <- cell_anno$.
    allgenes <-  rownames(SeuratObj@assays$RNA$data)

    RecTFDB <- readRDS("/home/jiawen/myMLnet/pathways/Nichenet/RTF_mouse.rds") %>% distinct(from, to, .keep_all = T) %>% Filter_DB(allgenes)
    # RecTFDB <- readRDS("/home/jiawen/myMLnet/pathways/KEGG/RTF_human.rds")
    TFTGDB <- readRDS("/home/jiawen/myMLnet/pathways/Nichenet/TFT_mouse.rds") %>% distinct(from, to, .keep_all = T) %>% Filter_DB(allgenes)

    LigRec_original <- Infer_CCI(SeuratObj, cellchat_output = T, db_use = "mouse")
    args.project<- 'GSE141259'
    args.target <- paste0("t",i %>% as.character())

    target_type <- "Krt8 ADI"
    TG_used <- FindMarkers(SeuratObj, target_type, min.pct = 0.25, only.pos = T, logfc.threshold = 0.25)
    TG_used <- filter(TG_used, p_val<1e-3) %>% rownames()
    input_dir <- "/home/jiawen/myMLnet/pythoncodes/inputs"
    Prepare_Input_New(SeuratObj, target_type, TGs = TG_used, CCC_results = LigRec_original, RTFDB = RecTFDB, TFTGDB = TFTGDB, data_dir = input_dir,
        assay = "RNA", datatype = "scale.data", imputation = F,exp_threshold = 0.1,
        CCC_threshold = 0.1, Fisher_threshold = 1,species = "mouse")
    conda_python <- "/home/jiawen/anaconda3/envs/SigXTalk/bin/python"

    system2(conda_python, args = c("/home/jiawen/myMLnet/pythoncodes/main_new.py", paste("--project",shQuote(args.project)), paste("--target_type",args.target)))
    gc()
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

}
