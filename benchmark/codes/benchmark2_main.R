##### Script for benchmark PART II

# 0. load methods and datasets
library(dplyr)
library(Seurat)
setwd("/home/jiawen/myMLnet/benchmark")
source("./codes/benchmark2_methods_func.R")
source("./codes/benchmark2_evaluation_func.R")
source("/home/jiawen/myMLnet/codes/MLNconstruct.R")
source("/home/jiawen/myMLnet/codes/Numericals.R")
source("/home/jiawen/myMLnet/codes/Utilities.R")
source("/home/jiawen/myMLnet/codes/Crosstalk_analysis.R")

input_path <- "/home/jiawen/myMLnet/datasets/nichenet/seurat.rds"
sampleID <- "HNSCC"
SeuratObj <- readRDS(input_path)
SeuratObj <- NormalizeData(SeuratObj)
table(SeuratObj$celltype)

args <- commandArgs(trailingOnly = TRUE)
target_type <- as.character(args[1])
# target_type <- "T_cell"

target_genes <- FindMarkers(SeuratObj, ident.1 = target_type, min.pct = 0.25, only.pos = T)
target_genes <- target_genes[target_genes$p_val_adj < 0.01, ]
target_genes <- target_genes %>%
  rownames() %>%
  unlist() %>%
  unique()

out_path <- paste0('./result_benchmark2/',sampleID)
if (!dir.exists(out_path)) {
  dir.create(out_path)
}

# 1. Run CytoTalk
result_cytotalk <- CytoTalk_function(SeuratObj, receiver = target_type, species = "human")
result_cytotalk[,2] %>% unique() %>% length()
saveRDS(result_cytotalk, file = paste0(out_path,"/Cytotalk_",target_type,".rds"))

# 2. Run Nichenet
result_Nichenet1 <- NicheNet_database_function(SeuratObj, geneset = target_genes, receiver = target_type)
result_Nichenet1$TG %>% unique() %>% length()
saveRDS(result_Nichenet1, file = paste0(out_path,"/Nichenet_db_",target_type,".rds"))

result_Nichenet2 <- NicheNet_computational_function(SeuratObj, geneset = target_genes, receiver = target_type)
result_Nichenet2 <- filter(result_Nichenet2,result_Nichenet2$value > 0) %>% distinct(Receptor, TG ,.keep_all = T)
result_Nichenet2$TG %>% unique() %>% length()

saveRDS(result_Nichenet2, file = paste0(out_path,"/Nichenet_comp_",target_type,".rds"))

# 3. Run SigXTalk
result_SigXTalk <- SigXTalk_computational_function(SeuratObj, geneset = target_genes, receiver = target_type)
result_SigXTalk$TF %>% unique() %>% length()
saveRDS(result_SigXTalk, file = paste0(out_path,"/SigXTalk_",target_type,".rds"))


