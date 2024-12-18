### Generate simulation data
gc()
setwd("/home/jiawen/myMLnet/benchmark")
source("/home/jiawen/myMLnet/benchmark/codes/Simulation_SERGIO.r")

require(dplyr)
require(rlang)

args <- commandArgs(trailingOnly = TRUE)
no.rec <- 10; 
no.tf <- 40; 
no.tg <- 150;
ncell_bin <- 500

# no.rec <- ; no.tf <- args[2]; no.tg <- args[3];


p1 = as.numeric(args[1]); p2 = as.numeric(args[1]);
n_bins = 3;

GRN_original <- Generate_Random_MLnet(no.rec,no.tf,no.tg,p1,p2)
GRN_sergio <- GRN_Output_For_SERGIO(GRN_original,n_bins,workdir = "./datasets/")
GRN_sergio <- GRN_sergio[,1:3]

n_genes <- no.rec+no.tf+no.tg

conda_python <- "/home/jiawen/anaconda3/envs/SigXTalk/bin/python"
system2(conda_python, args = c("./SERGIO/run_sergio.py", n_genes, n_bins, ncell_bin))

# Export GRN
RecTF <- GRN_sergio[which(GRN_sergio$From %in% 0:(no.rec-1)),1:2]
TFTG <- GRN_sergio[which(GRN_sergio$To %in% (no.rec+no.tf):(no.rec+no.tf+no.tg-1)), 1:2]

GRN_sergio$From <- lapply(GRN_sergio$From, function(x) paste0('Gene', x)) %>% unlist()
GRN_sergio$To <- lapply(GRN_sergio$To, function(x) paste0('Gene', x)) %>% unlist()

GRN_sergio$To <- unlist(GRN_sergio$To)

RecTF <- as.data.frame(apply(RecTF, 2, function(x) paste0('Gene', x)))
TFTG <- as.data.frame(apply(TFTG, 2, function(x) paste0('Gene', x)))

write.csv2(GRN_sergio %>% as.data.frame(), file = "/home/jiawen/myMLnet/benchmark/datasets/sergio_groundtruth.csv",quote = F,row.names = F)
write.csv2(RecTF, file = "/home/jiawen/myMLnet/benchmark/datasets/sergio_RecTF.csv",quote = F,row.names = F)
write.csv2(TFTG, file = "/home/jiawen/myMLnet/benchmark/datasets/sergio_TFTG.csv",quote = F,row.names = F)

print("Simulation succeed!\n")

require(Seurat)
require(dplyr)
rna_exp <- read.csv("/home/jiawen/myMLnet/benchmark/datasets/sergio_counts.csv",header = T) %>% as.data.frame()
rownames(rna_exp) <- rna_exp$X
rna_exp <- rna_exp[,-1] %>% t()
SeuratObj <- CreateSeuratObject(rna_exp)
c1 <- data.frame(cell = paste0("Cell",0:(ncell_bin-1)),cluster = "Cluster1")
c2 <- data.frame(cell = paste0("Cell",(ncell_bin):(ncell_bin*2-1)),cluster = "Cluster2")
c3 <- data.frame(cell = paste0("Cell",(ncell_bin*2):(ncell_bin*3-1)),cluster = "Cluster3")
cell_anno <- rbind(c1,c2,c3) %>% as.data.frame()
SeuratObj$cell_label <- Idents(SeuratObj) <- cell_anno$cluster
SeuratObj <- NormalizeData(SeuratObj) %>% SCTransform()

SeuratObj <- RunPCA(SeuratObj, features = VariableFeatures(object = SeuratObj))
SeuratObj <- RunUMAP(SeuratObj,dims = 1:10)
DimPlot(SeuratObj, reduction = "umap")

coords <- SeuratObj@reductions$umap@cell.embeddings %>% as.data.frame()

rna_data <- SeuratObj@assays$RNA$data
write.table(rna_data, file = "/home/jiawen/myMLnet/benchmark/datasets/sergio_data.csv",quote = F,row.names = T, col.names = T,sep = "\t")
write.table(coords, file = "/home/jiawen/myMLnet/benchmark/datasets/sergio_coord.csv",quote = F,row.names = T,sep = "\t")
write.csv2(cell_anno, file = "/home/jiawen/myMLnet/benchmark/datasets/sergio_anno.csv",quote = F,row.names = F)

# Generate the pseudotime

system2(conda_python, args = c("./codes/pseudo_time.py", n_genes, n_bins, ncell_bin))

# Split the datasets
rna_rectf <- rna_data[1:(no.rec+no.tf),]
rna_tftg <- rna_data[(no.rec+1):(no.rec+no.tf+no.tg),]
write.csv(rna_rectf, file = "/home/jiawen/myMLnet/benchmark/Beeline/inputs/example/bench_rectf/ExpressionData.csv",quote = F,row.names = T)
write.csv(rna_tftg, file = "/home/jiawen/myMLnet/benchmark/Beeline/inputs/example/bench_tftg/ExpressionData.csv",quote = F,row.names = T)

gc()
