##### Analysis of the results
# 0. Load the network
library(tidyverse)
library(Seurat)
setwd("/home/jiawen/myMLnet/benchmark")
source("./codes/benchmark2_methods_func.R")
source("./codes/benchmark2_evaluation_func.R")
source("./codes/benchmark2_generateGT.R")
source("/home/jiawen/myMLnet/codes/Utilities.R")

input_path <- "/home/jiawen/myMLnet/datasets/hnscc/seurat.rds"
sampleID <- "HNSCC"
SeuratObj <- readRDS(input_path)
SeuratObj <- NormalizeData(SeuratObj)
table(SeuratObj$celltype)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
target_type <- as.character(args[1])
db_type <- as.character(args[2])
# target_type <- "Endothelial"
# db_type <- "KEGG"
out_path <- paste0('./result_benchmark2/',sampleID)
if (!dir.exists(out_path)) {
  dir.create(out_path)
}

Ncores <- parallel::detectCores() -1
Exp_clu <- Get_Exp_Clu(SeuratObj,clusterID = target_type, assay = "RNA", datatype = "counts", cutoff = 0.1)
allgenes <- rownames(Exp_clu)
receptors <- readRDS("/home/jiawen/myMLnet/pathways/benchmark/Receptors.rds")
tfactors <- readRDS("/home/jiawen/myMLnet/pathways/benchmark/TFs.rds")

if(db_type == "KEGG"){
  GT_graph <- readRDS("/home/jiawen/myMLnet/pathways/benchmark/KEGG_graph.rds")
  filen_result <- paste0("/home/jiawen/myMLnet/benchmark/result_benchmark2/figures/KEGG_",target_type,".csv")
  filen_figure <- paste0("/home/jiawen/myMLnet/benchmark/result_benchmark2/figures/figure_KEGG_",target_type,".pdf")
  GT_path <- Multi_PPR(GT_graph, receptors, tfactors)
  GT_path <- GT_path %>% filter(tf_PPR > 0)
  GT_path <- GT_path[,c(1,3)]
}else{
  # GT_graph <- readRDS("/home/jiawen/myMLnet/pathways/benchmark/REACTOME_graph.rds")
  # GT_graph <- GT_graph[-c(46, 47)] # filter out non-receptor links
  GT_path <- readRDS("/home/jiawen/myMLnet/pathways/benchmark/REACTOME_truth.rds")
  GT_path <- GT_path[,c(1,2)]
  filen_result <- paste0("/home/jiawen/myMLnet/benchmark/result_benchmark2/figures/REACTOME_",target_type,".csv")
  filen_figure <- paste0("/home/jiawen/myMLnet/benchmark/result_benchmark2/figures/figure_REACTOME_",target_type,".pdf")
}

colnames(GT_path) <- c("receptor","tf")

GT_new <- filter(GT_path,receptor %in% allgenes)
GT_new <- filter(GT_new,tf %in% allgenes)

GT_new <- Calculate_Coexp(GT_new, exp = Exp_clu, method = "overlap", Ncores = Ncores)
coexp_thres <- 0.2
GT_filtered <- filter(GT_new, coexp > coexp_thres)
GT_filtered <- GT_filtered[,1:2]
allrec <- GT_filtered$from %>% unique() %>% sort()
alltf <- GT_filtered$to %>% unique() %>% sort()

# 1. Analysis of CytoTalk
result <- readRDS(paste0(out_path,"/Cytotalk_",target_type,".rds"))
colnames(result) <- c("from","to","value")
df_Cytotalk <- Evaluate_method(result, GT_filtered, allrec, allgenes)
df_Cytotalk$Method <- "CytoTalk"

# 2. Analysis of Nichenet
result <- readRDS(paste0(out_path,"/Nichenet_db_",target_type,".rds"))
colnames(result) <- c("from","to","value")
df_Nichenet <- Evaluate_method(result, GT_filtered, receptors, allgenes)
df_Nichenet$Method <- "Nichenet-db"

result <- readRDS(paste0(out_path,"/Nichenet_comp_",target_type,".rds"))
colnames(result) <- c("from","to","value")
df_Nichenet2 <- Evaluate_method(result, GT_filtered, receptors, allgenes)
df_Nichenet2$Method <- "Nichenet-comp"

# 3. Analysis of SigXTalk
result <- readRDS(paste0(out_path,"/SigXTalk_",target_type,".rds"))
colnames(result) <- c("from","to","value")
# result <- filter(result, value > 0.9)
df_SigXTalk <- Evaluate_method(result, GT_filtered, receptors, allgenes)
df_SigXTalk$Method <- "SigXTalk"

df_bind <- rbind(df_Cytotalk, df_Nichenet,df_Nichenet2,df_SigXTalk)
colnames(df_bind) <- c("Coverage", "-lgFisher","-lgChisq","Method")

write.table(df_bind, file = filen_result, quote = F)

t.test(df_SigXTalk$Coverage, df_Nichenet$Coverage)
t.test(df_SigXTalk$Coverage, df_Nichenet2$Coverage)

library(ggplot2)
library(patchwork) # To display 2 charts together
library(hrbrthemes)

df1 <- data.frame(Method = df_bind$Method, Metric = "Coverage", value = df_bind$Coverage)
df2 <- data.frame(Method = df_bind$Method, Metric = "-lgFisher", value = df_bind$"-lgFisher")
df3 <- data.frame(Method = df_bind$Method, Metric = "-lgChisq", value = df_bind$"-lgChisq")

df <- rbind(df1,df2,df3)
allm <- df$Method %>% unique() %>% sort(decreasing = T)
df %>% 
    mutate(name = factor(Method, levels=allm)) %>%
    mutate(Metric = factor(Metric, levels = c("Coverage","-lgFisher","-lgChisq"))) %>% 
    ggplot( aes(x=name, y=value, fill=Metric)) + 
    geom_boxplot() +
    xlab("Method") +
    ylab("Value") +
    scale_y_continuous()+
    theme_bw() +
    theme(axis.title = element_text(size = 24)) +
    theme(axis.text = element_text(size = 24)) +
    theme(legend.text = element_text(size = 24))+
    theme(legend.title = element_text(size = 24)) +
    theme(legend.position="top") 

ggsave(filen_figure,device = cairo_pdf,width = 3600,height = 2400,unit = "px")
