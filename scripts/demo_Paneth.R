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

db <- CellChatDB.mouse
recs_ccc <- db$interaction$receptor

TF_PPR_KEGG_mouse <- read.csv("/home/jiawen/myMLnet/pathways/KEGG/KEGG_PPR_mouse.csv")
paths <- filter(TF_PPR_KEGG_mouse, pathway == "Wnt signaling pathway")
all_recs <- paths$receptor %>% unique()

PCP_recs <- c("Ror1","Ror2","Ryk","Gpc4","Vangl1") %>% sort()
common_recs <- lapply(1:10, function(x){
  paste0("Fzd",x)
}) %>% unlist() %>% sort()
Canonical_recs <- c("Lrp5","Lrp6")

# The TF-TG database
TFTGDB <- readRDS("/home/jiawen/myMLnet/pathways/Nichenet/TFT_mouse.rds")

## FVR cells

SeuratObj <- readRDS("/home/jiawen/myMLnet/datasets/Paneth/All.rds")
table(SeuratObj$sample)

sample_name <- "ct7"
if(sample_name == "ct7"){
  subobj <- subset(SeuratObj, sample == "Control_7_FVR_only")
}
if(sample_name == "ct3"){
  subobj <- subset(SeuratObj, sample == "Control_3_FVR")
}
if(sample_name == "ct5"){
  subobj <- subset(SeuratObj, sample == "Control_5_FVR")
}
if(sample_name == "mutant"){
  subobj <- subset(SeuratObj, sample == "Mutant_1")
}

meta <- subobj@meta.data %>% as.data.frame()
table(meta$proma_cell_type)
Idents(subobj) <- subobj$cluster <- subobj$proma_cell_type
subobj <- NormalizeData(subobj)
target_type <- "Paneth progenitor"


# cell_anno <- data.frame(cell = names(Idents(SeuratObj)), cluster = Idents(SeuratObj) %>% as.character())
# LR_original <- Infer_CCI(SeuratObj, cell_anno, LRDB = LRDB, cellchat_output = T, use_spatial = F, db_use = "mouse")
# LR_Pairprob <- Extract_LR_Prob(LR_original, target_type = target_type, cellchat_use = T)
# LR_Pairprob$From <- stringr::str_to_title(LR_Pairprob$From)
# LR_Pairprob$To <- stringr::str_to_title(LR_Pairprob$To)
# LR_Pairprob$To %>% unique() %>% sort()


markers <- FindMarkers(subobj, ident.1 = "ISC",only.pos = T, logfc.threshold = 0.25, min.pct = 0.25)
markers <- filter(markers, p_val_adj < 0.01) %>% rownames_to_column()
TGs <- c("Cfap126","Lgr5","Fzd2","Foxa2")
TGs %in% markers$rowname
VlnPlot(subobj, features = c(TGs, "Lrp6","Lrp5"))

Exp_clu <- Get_Exp_Clu(subobj, clusterID = target_type, assay = "RNA", datatype = "counts", cutoff = 0.05)
allgenes <- rownames(Exp_clu)
TFTGDB_filtered <- TFTGDB %>% distinct(from, to, .keep_all = T) %>% Filter_DB(allgenes)
TFTG_filtered <- FindRegulator(TFTGDB_filtered, exp = Exp_clu, selected_genes = TGs,
  pv_threshold = 1)

PCP_recs  <- intersect(PCP_recs, allgenes)
common_recs <- intersect(common_recs, allgenes)
Canonical_recs <- intersect(Canonical_recs, allgenes)
all_recs <- c(PCP_recs,common_recs,Canonical_recs) %>% unique()

PCP_tfs <- lapply(PCP_recs, function(x){
  temp <- paths %>% filter(receptor == x) %>% dplyr::select(tf)
  temp
}) %>% unlist() %>% unique() %>% sort()

common_tfs <- lapply(common_recs, function(x){
  temp <- paths %>% filter(receptor == x) %>% dplyr::select(tf)
  temp
}) %>% unlist() %>% unique() %>% sort()

Canonical_tfs <- lapply(Canonical_recs, function(x){
  temp <- paths %>% filter(receptor == x) %>% dplyr::select(tf)
}) %>% unlist() %>% unique() %>% sort()

TFTG_filtered <- filter(TFTG_filtered, TF %in% common_tfs)
paths_filtered <- filter(paths, tf %in% TFTG_filtered$TF) %>% filter(receptor %in% allgenes)
paths_filtered <- data.frame(Receptor = paths_filtered$receptor, TF = paths_filtered$tf)
RecTFTG_all <- inner_join(paths_filtered, TFTG_filtered, by = "TF")
colnames(RecTFTG_all) <- c("Receptor","SSC","Target")
RecTFTG_all <- RecTFTG_all %>% distinct(Receptor, SSC, Target, .keep_all = T)
RecTFTG_all <- filter(RecTFTG_all, Receptor %in% all_recs)

Exp_clu <- Get_Exp_Clu(subobj, clusterID = target_type, assay = "RNA", datatype = "counts", cutoff = 0.05)
Exp_clu[1:5,1:5]
ress <- Rec_To_TFTG(Exp_clu, RecTFTG_all, method = "rf", cutoff = 0.1, use_tidy = F)
ress_filtered <- filter(ress, Weight > 0.05*max(ress$Weight))

tt <- gsub(' ',"_",target_type)
filen <- paste0("/home/jiawen/myMLnet/results/Paneth/result_",sample_name,'_',tt,".csv")
write.table(ress, file = filen,quote = F, sep = ",")

aggregate(ress$Weight, list(ress$Receptor), sum)
