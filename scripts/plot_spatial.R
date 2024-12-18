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
library(CellChat)

setwd("/home/jiawen/myMLnet")
source("./codes/MLNconstruct.R")
source("./codes/Numericals.R")
source("./codes/Utilities.R")
source("./codes/Crosstalk_analysis.R")
source("./codes/Visualization.R")

SeuratObj <- readRDS("./datasets/visium/cortex_new.rds")
SeuratObj@images$anterior1@scale.factors$spot <- 120
cell_anno <- Idents(SeuratObj) %>%
    as.data.frame() %>%
    rownames_to_column()
SeuratObj$cluster <- cell_anno$.

scale_factors <- list(
    spot_diameter_fullres = 89.47789776697327,
    tissue_hires_scalef = 0.17211704,
    fiducial_diameter_fullres = 144.54121946972606,
    tissue_lowres_scalef = 0.051635113
)

LigRec_original <- Infer_CCI(SeuratObj, cell_anno,
    LRDB = NULL,
    cellchat_output = T, use_spatial = T,
    scale_factors = scale_factors, db_use = "mouse"
)

target_type <- "L6 IT"

pathways <- LigRec_original@netP$pathways

pdf("./results/spatial/figure_new/fig6a.pdf")
netVisual_aggregate(LigRec_original,
    signaling = pathways, layout = "spatial",
    edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 8
)
dev.off()

source_type <- "Meis2"
LR_Pairprob <- Extract_LR_Prob(LigRec_original, source_type = source_type, target_type = target_type, cellchat_use = T)
LR_Pairprob$From <- str_to_title(LR_Pairprob$From)
LR_Pairprob$To <- str_to_title(LR_Pairprob$To)
Rec_act <- aggregate(LR_Pairprob$Weight, list(LR_Pairprob$To), sum)
colnames(Rec_act) <- c("Rec", "Weight")
Rec_act <- Rec_act[Rec_act$Weight > 0.1*max(Rec_act$Weight),]
LR_Pairprob <- LR_Pairprob[LR_Pairprob$To %in% Rec_act$Rec,]
pdf("./results/spatial/figure_new/fig6b.pdf")
PlotCCI_CirclePlot(result = LR_Pairprob, topk = 10)
dev.off()

filen <- "/home/jiawen/myMLnet/results/spatial/results_L6_IT.csv"
temp <- read.table(filen,sep = ',',header = T)
temp <- temp[!grepl("^mt-", temp$Target), ]
temp <- temp[!grepl("^Rpl", temp$Target), ]
temp <- temp[!grepl("^Rps", temp$Target), ]
CC_results <- filter(temp, Weight > 0.01*max(temp$Weight))
RTG_pair_results <- Aggregate_Causality(CC_results, sum_type = "sum", data_type = "Target")

table(RTG_pair_results$Receptor)
Receptor_used <- c("Ncl","Sdc3","Ptprz1","Ntrk3")
Spe_all <- Calculate_Pairwise(RTG_pair_results, KeyGene = Receptor_used, type = "Spe")
PlotXT_MultiCircularBar(df = Spe_all[c('Target','Receptor','Specificity')], KeyFactors = Receptor_used)
ggsave("./results/spatial/figure_new/figc.pdf", height = 3000,width = 4200,units = "px")

SeuratObj$celltype <- Idents(SeuratObj)
subobject <- subset(SeuratObj, cluster == target_type)


KeyRecs <- c("Ntrk3")
Subcluster_df <- XT_Subclustering(subobject, RTG_pair_results, KeyRecs, target_type)
PlotXT_Subcluster(SeuratObj, Subcluster_df, source_type, target_type)
ggsave("./results/spatial/figure_new/supp.pdf", height = 3000,width = 4200,units = "px")

PlotXT_distances(SeuratObj, Subcluster_df, source_type, target_type)
ggsave("./results/spatial/figure_new/fige.pdf", device = cairo_pdf, height = 3000,width = 1500,units = "px")

# Overview of the data
pdf(file = "./results/spatial/figure_new/suppa.pdf",width = 15,height = 15)
SpatialDimPlot(SeuratObj, label = TRUE, label.size = 6)
dev.off()

# other receptors
KeyRecs <- c("Ncl")
Subcluster_df <- XT_Subclustering(subobject, RTG_pair_results, KeyRecs, target_type)
PlotXT_Subcluster(SeuratObj, Subcluster_df, source_type, target_type)
ggsave("./results/spatial/figure_new/suppb1.pdf", height = 3000,width = 4200,units = "px")

PlotXT_distances(SeuratObj, Subcluster_df, source_type, target_type)
ggsave("./results/spatial/figure_new/suppb2.pdf", device = cairo_pdf, height = 3000,width = 1500,units = "px")

KeyRecs <- c("Sdc3")
Subcluster_df <- XT_Subclustering(subobject, RTG_pair_results, KeyRecs, target_type)
PlotXT_Subcluster(SeuratObj, Subcluster_df, source_type, target_type)
ggsave("./results/spatial/figure_new/suppc1.pdf", height = 3000,width = 4200,units = "px")

PlotXT_distances(SeuratObj, Subcluster_df, source_type, target_type)
ggsave("./results/spatial/figure_new/suppc2.pdf", device = cairo_pdf, height = 3000,width = 1500,units = "px")

KeyRecs <- c("Ptprz1")
Subcluster_df <- XT_Subclustering(subobject, RTG_pair_results, KeyRecs, target_type)
PlotXT_Subcluster(SeuratObj, Subcluster_df, source_type, target_type)
ggsave("./results/spatial/figure_new/suppd1.pdf", height = 3000,width = 4200,units = "px")

PlotXT_distances(SeuratObj, Subcluster_df, source_type, target_type)
ggsave("./results/spatial/figure_new/suppd2.pdf", device = cairo_pdf, height = 3000,width = 1500,units = "px")
