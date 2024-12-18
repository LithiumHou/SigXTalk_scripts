gc()
setwd("/home/jiawen/myMLnet/benchmark")
source("/home/jiawen/myMLnet/benchmark/codes/Simulation_SERGIO.r")

source("/home/jiawen/myMLnet/codes/MLNconstruct.R")
source("/home/jiawen/myMLnet/codes/Numericals.R")
source("/home/jiawen/myMLnet/codes/Utilities.R")
source("/home/jiawen/myMLnet/codes/Crosstalk_analysis.R")
require(dplyr)
require(rlang)
require(igraph)
# args <- commandArgs(trailingOnly = TRUE)
no.rec <- 10; no.tf <- 40; no.tg <- 150; ncell_bin <- 500;

p1 = 0.25; p2 = 0.2;
n_bins = 3;

GRN_original <- Generate_Random_MLnet(no.rec,no.tf,no.tg,p1,p2)
GRN_sergio <- GRN_Output_For_SERGIO(GRN_original,n_bins,workdir = "./datasets/")
GRN_sergio <- GRN_sergio[,1:3]

n_genes <- no.rec+no.tf+no.tg

conda_python <- "/home/jiawen/anaconda3/envs/SigXTalk/bin/python"
system2(conda_python, args = c("./SERGIO/run_sergio.py", n_genes, n_bins, ncell_bin))

# construct the TFTGDB
Recs <- 0:(no.rec-1) %>% as.character()
TFs <- no.rec:(no.rec+no.tf-1) %>% as.character()
RecTFDB <- filter(GRN_original, From %in% Recs)
TFTGDB <- filter(GRN_original, From %in% TFs)

genes <- c(Recs,TFs)
gg <- graph_from_data_frame(RecTFDB, directed = T, vertices = genes)
plot(gg)

nod <- sapply(igraph::V(gg)$name, function(z) strsplit(strsplit(z,":")[[1]][1],","))
r <- nodes_find(Recs, gg, nod)
tf <- nodes_find(TFs, gg, nod)
RTF_filtered <- myPageRank(gg, Recs, TFs)
colnames(RTF_filtered) <- c("Receptor","TF","PPR")
RTF_filtered <- na.omit(RTF_filtered)
RTF_filtered <- RTF_filtered[RTF_filtered$PPR>0,]

RTF_filtered$Receptor <- sapply(RTF_filtered$Receptor, function(x) paste0('Gene', x)) %>% unlist()
RTF_filtered$TF <- sapply(RTF_filtered$TF, function(x) paste0('Gene', x)) %>% unlist()
RTF_filtered$PPR <- as.numeric(RTF_filtered$PPR)
TFTGDB <- as.data.frame(apply(TFTGDB, 2, function(x) paste0('Gene', x)))

table(RTF_filtered$TF)
write.csv2(GRN_sergio %>% as.data.frame(), file = "/home/jiawen/myMLnet/benchmark/datasets/sergio_groundtruth.csv",quote = F,row.names = F)
write.table(RTF_filtered, file = "/home/jiawen/myMLnet/benchmark/datasets/sergio_RecTF.csv",quote = F,row.names = F,sep = ";")
write.table(TFTGDB, file = "/home/jiawen/myMLnet/benchmark/datasets/sergio_TFTG.csv",quote = F,row.names = F,sep = ";")
