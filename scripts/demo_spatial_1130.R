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

SeuratObj <- readRDS("./datasets/visium/cortex_new.rds")
cell_anno <- Idents(SeuratObj) %>%
    as.data.frame() %>%
    rownames_to_column()
SeuratObj$cluster <- cell_anno$.
allgenes <-  rownames(SeuratObj@assays$Spatial$data)
RecTFDB <- readRDS("/home/jiawen/myMLnet/pathways/Nichenet/RTF_mouse.rds") %>% distinct(from, to, .keep_all = T) %>% Filter_DB(allgenes)
# RecTFDB <- readRDS("/home/jiawen/myMLnet/pathways/KEGG/RTF_human.rds")
TFTGDB <- readRDS("/home/jiawen/myMLnet/pathways/Nichenet/TFT_mouse.rds") %>% distinct(from, to, .keep_all = T) %>% Filter_DB(allgenes)

scale_factors <- list(
    spot_diameter_fullres = 89.47789776697327,
    tissue_hires_scalef = 0.17211704,
    fiducial_diameter_fullres = 144.54121946972606,
    tissue_lowres_scalef = 0.051635113
)

LigRec_original <- Infer_CCI(SeuratObj, cell_anno,
    LRDB = LRDB,
    cellchat_output = T, use_spatial = T,
    scale_factors = scale_factors, db_use = "mouse"
)
groupSize <- as.numeric(table(LigRec_original@idents))
netVisual_heatmap(LigRec_original, measure = "weight", color.heatmap = "Blues")
pathways <- LigRec_original@netP$pathways
netVisual_aggregate(LigRec_original,
    signaling = pathways, layout = "spatial",
    edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5
)
args.project<- 'spatial'
args.target <- 'L6_IT'
target_type <- "L6 IT"

TG_used <- FindMarkers(SeuratObj, target_type, min.pct = 0.25, only.pos = T, logfc.threshold = 0.25)
TG_used <- filter(TG_used, p_val<1e-3) %>% rownames()
input_dir <- "/home/jiawen/myMLnet/pythoncodes/inputs"
Prepare_Input_New(SeuratObj, target_type, TGs = TG_used, CCC_results = LigRec_original, RecTFDB, TFTGDB, data_dir = input_dir,
    assay = "Spatial", datatype = "scale.data", imputation = F,exp_threshold = 0.1,
    CCC_threshold = 0.1, Fisher_threshold = 1,species = "mouse")
gc()
conda_python <- "/home/jiawen/anaconda3/envs/SigXTalk/bin/python"
system2(conda_python, args = c("/home/jiawen/myMLnet/pythoncodes/main_new.py", paste("--project",shQuote(args.project)), paste("--target_type",args.target)))

output_dir <- '/home/jiawen/myMLnet/results/'
filen <- paste0(output_dir, args.project,'/pathways_',args.target,'.csv')
RTFTG_results <- read.csv(filen, header = T)
RTFTG_results <- RTFTG_results[RTFTG_results$pred_label > 0.75, ]
RTFTG_results <- RTFTG_results[,1:3]
RTFTG_results$TG %>% unique()
Exp_clu <- Get_Exp_Clu(SeuratObj, clusterID = target_type, assay = "Spatial", datatype = "data", cutoff = 0.05)
ress <- Rec_To_TFTG(Exp_clu, RTFTG_results, method = "rf", cutoff = 0.1, use_tidy = T)
filen <- paste0(output_dir, args.project,'/results_',args.target,'.csv')
write.table(ress,file = filen, quote = F, sep = ",")

