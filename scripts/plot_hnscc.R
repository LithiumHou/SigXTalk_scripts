rm(list = ls())
gc()
###### Analyze the results
library(Matrix)
library(dplyr)
library(tidyfst)
library(ggsci)
library(ggplot2)
library(igraph)
library(plotrix)
library(ggraph)
#library(org.Hs.eg.db)
#library(clusterProfiler)
library(ggalluvial)
library(tidyverse)
library(circlize)
library(RColorBrewer)
library(viridis)
library(hrbrthemes)
library(Seurat)
library(ggridges)
library(ComplexHeatmap)
library(patchwork)
setwd("/home/jiawen/myMLnet")
source("./codes/MLNconstruct.R")
source("./codes/Numericals.R")
source("./codes/Utilities.R")
source("./codes/Crosstalk_analysis.R")
source("./codes/Visualization.R")

# 0. Load the data and results
SeuratObj <- readRDS("./datasets/hnscc/seurat1006.rds")
target_type <- "malignant"
filen <- paste0('./results/hnscc/results_',target_type,".csv")
results <- read.table(filen,sep = ',',header = T)
results <- results[!grepl("^MT", results$Target), ]
results <- results[!grepl("^RPL", results$Target), ]
results <- results[!grepl("^RPS", results$Target), ]
max(results$Weight)
sum(results$Weight < 0.05)
results_filtered <- filter(results, Weight > 0.05*max(results$Weight))
results_filtered$Target %>% unique() %>% length()
Exp_clu <- Get_Exp_Clu(SeuratObj, clusterID = target_type, datatype = "counts")

SeuratObj[["percent.mt"]] <- PercentageFeatureSet(SeuratObj, pattern = "^MT")
VlnPlot(SeuratObj, features = "percent.mt")

CC_pair_results <- Aggregate_Causality(results_filtered, sum_type = "sum",data_type = "Target")
cell_anno <- data.frame(cell = names(Idents(SeuratObj)), cluster = Idents(SeuratObj) %>% as.character())
LR_original <- Infer_CCI(SeuratObj, LRDB = NULL, cellchat_output = T, use_spatial = F, db_use = "human")

# 1. CCC analysis and plot - overall
groupSize <- as.numeric(table(LR_original@idents))
pdf(file = "./results/hnscc/figures_new/fig3a.pdf",width = 10,height = 10)
CellChat::netVisual_circle(LR_original@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 vertex.label.cex = 2)
dev.off()

# 2. CCC analysis - targeting malignant
CCC_threshold <- 0.1
LR_Pairprob <- Extract_LR_Prob(LR_original, target_type = target_type, cellchat_use = T)
LR_Pairprob <- LR_Pairprob[which(LR_Pairprob$Weight >= CCC_threshold * max(LR_Pairprob$Weight)), ]
pdf(file = "./results/hnscc/figures_new/fig3b.pdf",width = 10,height = 10)
PlotCCI_CirclePlot(LR_Pairprob, topk = 30)
dev.off()

# 3. CrossTalk - statistics
ps <- PlotXT_Counts(results_filtered, datatype = "Target", top_percent = 1)
ps$outer
ggsave(filename = "./results/hnscc/figures_new/fig3c_outer.pdf",device = cairo_pdf,width = 3000,height = 2400,unit = "px")
ps$inner
ggsave(filename = "./results/hnscc/figures_new/fig3c_inner.pdf",device = cairo_pdf,width = 2000,height = 1800,unit = "px")

# 4. CrossTalk - a heatmap of most-expressed target gene
Counts_pathway <- Count_Crosstalk(results_filtered, KeyGenes = NULL, verbose = F, datatype = "Target")
TG_used <- sort(Counts_pathway, decreasing = T) %>% names()
TG_used <- TG_used[1:15]
PlotXT_RecTGHeatmap(CC_pair_results, Exp_clu, TG_used = TG_used,topk = 100)
ggsave(filename = "./results/hnscc/figures_new/fig3d.pdf",device = cairo_pdf,width = 2000,height = 3600,unit = "px")

# 5. CrossTalk - Sankey
TG_used <- c("KRT5","KRT6A","KRT14","KRT17","EPCAM","SFN")
PlotXT_Alluvial(results_filtered, TG_used, min_weight = 0.98)
ggsave(filename = "./results/hnscc/figures_new/fig3e.pdf",device = cairo_pdf,width = 3000,height = 4000,unit = "px")

# 6. CrossTalk - single heatmap
TG_used <- "KRT5"
pdf(file = "./results/hnscc/figures_new/fig3f.pdf",width = 10,height = 10)
PlotXT_HeatMap(results_filtered, TG_used, "TG")
dev.off()

# 7. CrossTalk - fidelity distribution
TG_used <- "KRT17"
PlotXT_Ridgeline(results_filtered, TG_used)
ggsave(filename = "./results/hnscc/figures_new/fig3g.pdf",device = cairo_pdf,width = 3200,height = 3600,unit = "px")

# 8. CrossTalk - multiple marker
# Read the datasets
types_used <- c("CAF","Endothelial","malignant","myofibroblast","T_cell")

results_list <- lapply(types_used, function(x){
    filen <- paste0('./results/hnscc/results_',x,".csv")
    temp <- read.table(filen,sep = ',',header = T)
    temp <- temp[!grepl("^MT", temp$Target), ]
    temp <- temp[!grepl("^RPL", temp$Target), ]
    temp <- temp[!grepl("^RPS", temp$Target), ]
    filter(temp, Weight > 0.05*max(temp$Weight))
})
names(results_list) <- types_used
Receptor_list <- lapply(results_list,function(x){
    x$Receptor %>% unlist() %>% unique()
})
names(Receptor_list) <- types_used
Receptor_used <- reduce(Receptor_list[1:4],intersect)
Spe_all <- Calculate_Pairwise(CC_pair_results, KeyGene = Receptor_used, type = "Spe")
PlotXT_MultiCircularBar(df = Spe_all[c('Target','Receptor','Specificity')], KeyFactors = Receptor_used)
ggsave(filename = "./results/hnscc/figures_new/fig3h.pdf",device = cairo_pdf,width = 3200,height = 3600,unit = "px")

# 9. CrossTalk - Multiple type
KeyRec <- Receptor_used
PlotXT_MultiViolin(results_list[1:4], KeyRec)
ggsave(filename = "./results/hnscc/figures_new/fig3i.pdf",device = cairo_pdf,width = 3200,height = 3600,unit = "px")

# 10. cancer marker
SeuratObj <- readRDS("./datasets/hnscc/seurat1006.rds")
markers <- TG_used
markers %in% results_filtered$Target
VlnPlot(SeuratObj, features = markers)
ggsave(filename = "./results/hnscc/figures_new/fig3j.pdf",device = cairo_pdf,width = 4200,height = 3600,unit = "px")


# 12. Clustering
KeyRec <- "CDH3"
PlotXT_UMAP(results_filtered, KeyRec, Nk = 4)
ggsave(filename = "./results/hnscc/figures_new/fig3l.pdf",device = cairo_pdf,width = 4000,height = 3600,unit = "px")

# 13. Marker-related counts
counts_filtered <- data.frame(gene = markers, pathways = Counts_pathway[markers] %>% as.numeric())
counts_filtered$colororder <- match(counts_filtered$pathways, counts_filtered$pathways %>% unique())
counts_filtered$colors <- brewer.pal(counts_filtered$pathways %>% unique() %>% length(), "Paired")[counts_filtered $colororder]
counts_filtered %>%
    mutate(gene = fct_reorder(gene, pathways)) %>%
    ggplot(aes(x = gene, y = pathways, fill = gene)) +
    geom_bar(stat = "identity", alpha = .6, width = .4) +
    coord_flip() +
    theme(text = element_text(family = "Arial")) +
    xlab("Target gene") +
    ylab("Number of crosstalk pathways") +
    theme_bw() +
    theme(axis.title = element_text(size = 24)) +
    theme(axis.text = element_text(size = 24)) +
    theme(legend.position = "none")
ggsave(filename = "./results/hnscc/figures_new/fig3m.pdf",device = cairo_pdf,width = 4000,height = 3600,unit = "px")

# supp0 the data
SeuratObj <- RunPCA(SeuratObj) %>% RunUMAP(dims = 1:10)
SeuratObj$cluster %>% table()
Idents(SeuratObj)[1:5]
DimPlot(SeuratObj, reduction = "umap", label = F, pt.size = 1.5)
ggsave("./results/hnscc/figures_new/umap.pdf",device = cairo_pdf,width = 4000,height = 3600,unit = "px")

# supp1 CCCpathways
db <- CellChatDB.human$interaction
paths <- LR_original@netP$pathways %>% sort()
LR_Pairprob <- Extract_LR_Prob(LR_original, target_type = target_type, cellchat_use = T)
LR_Pairprob <- LR_Pairprob[which(LR_Pairprob$Weight >= 0.1 * max(LR_Pairprob$Weight)), ]
PlotCCI_CirclePlot(LR_Pairprob, topk = 30)

pathways.show <- c("COLLAGEN")
pdf(file = "./results/hnscc/figures_new/suppa.pdf",width = 10,height = 10)
netVisual_heatmap(LR_original, signaling = pathways.show,font.size = 20, font.size.title = 24, color.heatmap = "Reds")
dev.off()
pathways.show <- c("LAMININ")
pdf(file = "./results/hnscc/figures_new/suppb.pdf",width = 10,height = 10)
netVisual_heatmap(LR_original, signaling = pathways.show,font.size = 20, font.size.title = 24, color.heatmap = "Blues")
dev.off()

# supp2 histogram of No. pathways
df <- table(results_filtered$Target) %>% as.data.frame()
colnames(df) <- c("Target","Number")
df <- df[order(df$Number),] %>% rownames_to_column()
unique_values <- sort(unique(df$Number))
counts <- table(df$Number)
cumulative_probs <- cumsum(counts) / length(df$Number)
cdf_data <- data.frame(
  Value = unique_values,
  CDF = cumulative_probs
)

ggplot(cdf_data, aes(x = Value, y = CDF)) +
    geom_step() +  # Create step plot
    geom_point(size = 3, color = "blue") +  # Add points for clarity
    labs(
    x = "Number of crosstalk pathways",
    y = "Cumulative Probability"
    ) +
    scale_x_continuous(
        breaks = seq(0, 300, by = 50),      # Major ticks with labels
        minor_breaks = seq(0, 300, by = 10) # Minor ticks without labels 
    ) +
    theme(text = element_text(family = "Arial")) +
    theme(axis.text.x = element_text(size = 20, angle = 90, face = "bold", vjust = 0.5)) +
    theme(axis.text.y = element_text(size = 20, face = "bold")) +
    theme(axis.title.x = element_text(size = 24, face = "bold", vjust = 0.5)) +
    theme(axis.title.y = element_text(size = 24, face = "bold", vjust = 1)) +
    theme_bw()
ggsave(filename = "./results/hnscc/figures_new/suppc.pdf",device = cairo_pdf,width = 4000,height = 3600,unit = "px")

# supp3: Heatmap: FOS as the SSC
gene_used <- "FOS"
pdf(file = "./results/hnscc/figures_new/suppd.pdf",width = 10,height = 10)
PlotXT_HeatMap(results_filtered, gene_used, "TF")
dev.off()

# supp4: CDH3 as the signal
gene_used <- "CDH3"
df <- filter(results_filtered, Receptor == gene_used)
df$Spe <- df$Weight/sum(df$Weight)
df <- df[order(df$Spe, decreasing = T), ]
Spe_mat <- Calculate_Specificity_Matrix(results_filtered,KeyRec = gene_used)$all
df <- df[1:20, ]
mat <- tidyfst::df_mat(df, row = SSC, col = Target,value = Spe) %>% t()
ssc_used <- colnames(mat) %>% unique() %>% sort()
plot_order <-  c(sort(rownames(mat)),ssc_used)
cairo_pdf(file = "./results/hnscc/figures_new/suppe.pdf",width = 13,height = 10)
PlotXT_Chord(mat, orders = plot_order)
dev.off()

# supp4 correlation
Fid_all <- Calculate_Pairwise(CC_pair_results, KeyGene = CC_pair_results$Target %>% unique(), type = "Fid")
Spe_all <- Calculate_Pairwise(CC_pair_results, KeyGene = CC_pair_results$Receptor %>% unique(), type = "Spe")
df <- data.frame(Pathway = paste0(CC_pair_results$Receptor,'_',CC_pair_results$Target),Fidelity = Fid_all$Fidelity, Specificity = Spe_all$Specificity)
ggplot(df, aes(x=Fidelity, y=Specificity)) +
    geom_point() +
    theme_bw() +
    xlab("Fidelity") + 
    ylab("Specificity") +
    theme(text = element_text(family = "Arial")) +
    theme(axis.text.x = element_text(size = 20, angle = 90, face = "bold", vjust = 0.5)) +
    theme(axis.text.y = element_text(size = 20, face = "bold")) +
    theme(axis.title.x = element_text(size = 24, face = "bold", vjust = 0.5)) +
    theme(axis.title.y = element_text(size = 24, face = "bold", vjust = 1)) +
    theme(legend.title = element_text(size = 18, face = "bold", vjust = 1)) +
    theme(legend.text = element_text(size = 16)) 
ggsave(filename = "./results/hnscc/figures_new/suppf.pdf",device = cairo_pdf,width = 4000,height = 3600,unit = "px")
cor.test(df$Fidelity,df$Specificity)
# supp5 fid vs spe
KeyTG <- "KRT17"
PlotXT_FidSpe(results_filtered, KeyTG, threshold = 0.01)
ggsave(filename = "./results/hnscc/figures_new/suppg.pdf",device = cairo_pdf,width = 4200,height = 3000,unit = "px")

# suppX subcluster
SeuratObj <- RunPCA(SeuratObj) %>% RunUMAP(dims = 1:10)
SeuratObj$cluster %>% table()
SeuratObj <- FindNeighbors(SeuratObj, graph.name = "Wnn")
SeuratObj <- FindClusters(SeuratObj, graph.name = "Wnn", algorithm = 3, resolution = 1)
Idents(SeuratObj) <- SeuratObj$cluster

SeuratObj <- FindSubCluster(SeuratObj, cluster = "malignant",graph.name = "Wnn",algorithm = 3,resolution = 0.25)
meta <- SeuratObj@meta.data
Idents(SeuratObj) <- SeuratObj$sub.cluster

neworder <- c("myofibroblast", "CAF", "Fibroblast", "Endothelial", "Mast", 
 "T_cell", "Macrophage", "malignant_0", "malignant_1", "malignant_2",
"malignant_3", "malignant_4", "malignant_5")
SeuratObj$sub.cluster <- factor(SeuratObj$sub.cluster,levels = neworder)
Idents(SeuratObj) <- SeuratObj$sub.cluster


table(SeuratObj$sub.cluster)
subobj <- subset(SeuratObj, sub.cluster != "malignant_4")
subobj <- subset(subobj, sub.cluster != "malignant_5")
DimPlot(subobj, reduction = "umap", label = F, pt.size = 1.5)
ggsave("./results/hnscc/figures_new/umap_new.pdf",device = cairo_pdf,width = 4000,height = 3600,unit = "px")

subobj$main_cluster <- subobj$cluster
subobj$cluster <- subobj$sub.cluster
Idents(subobj) <- subobj$cluster
LR_subcluster <- Infer_CCI(subobj, LRDB = NULL, cellchat_output = T, use_spatial = F, db_use = "human")

groupSize <- as.numeric(table(LR_subcluster@idents))
cairo_pdf(file = "./results/hnscc/figures_new/sub_CCC_a.pdf",width = 10,height = 10)
CellChat::netVisual_circle(LR_subcluster@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 vertex.label.cex = 2)
dev.off()
cairo_pdf(file = "./results/hnscc/figures_new/sub_CCC_b.pdf",width = 10,height = 10)
CellChat::netVisual_circle(LR_subcluster@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength",
                 vertex.label.cex = 2)
dev.off()
pathways.show <- c("COLLAGEN")
netVisual_heatmap(LR_subcluster, signaling = pathways.show,font.size = 20, font.size.title = 24, color.heatmap = "Reds")
pathways.show <- c("LAMININ")
netVisual_heatmap(LR_subcluster, signaling = pathways.show,font.size = 20, font.size.title = 24, color.heatmap = "Reds")
