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

target_type <- "CAF"
db_type <- "KEGG"

out_path <- paste0('./result_benchmark2/',sampleID)
if (!dir.exists(out_path)) {
  dir.create(out_path)
}

Ncores <- parallel::detectCores() -1
Exp_clu <- Get_Exp_Clu(SeuratObj,clusterID = target_type, assay = "RNA", datatype = "counts", cutoff = 0.05)
allgenes <- rownames(Exp_clu)
receptors <- readRDS("/home/jiawen/myMLnet/pathways/benchmark/Receptors.rds")
tfactors <- readRDS("/home/jiawen/myMLnet/pathways/benchmark/TFs.rds")

GT_graph <- readRDS("/home/jiawen/myMLnet/pathways/benchmark/KEGG_graph.rds")
GT_path <- Multi_PPR(GT_graph, receptors, tfactors)
GT_path <- GT_path[order(GT_path$tf_PPR, decreasing = T),c(1,3,4)]
GT_path <- filter(GT_path, tf_PPR > 0)
GT_new <- filter(GT_path,receptor %in% allgenes)
GT_new <- filter(GT_new,tf %in% allgenes)
GT_new <- Calculate_Coexp(GT_new, exp = Exp_clu, method = "overlap", Ncores = Ncores)
GT_filtered <- filter(GT_new, coexp > 0.2)
n_path <- nrow(GT_filtered)

ppr_level <- c(0.2, 0.4, 0.6, 0.8, 1)
for (i in seq_along(ppr_level)){

    ni <- ceiling(n_path*ppr_level[i])
    GT_filtered <- GT_filtered[1:ni, 1:2]

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

    filen_result <- paste0("/home/jiawen/myMLnet/benchmark/result_benchmark2/robustness/ppr_level_",i,".csv")
    write.table(df_bind, file = filen_result, quote = F)

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

    filen_figure <- paste0("/home/jiawen/myMLnet/benchmark/result_benchmark2/figures/ppr_level_",i,".pdf")
    ggsave(filen_figure,device = cairo_pdf,width = 3600,height = 2400,unit = "px")
}