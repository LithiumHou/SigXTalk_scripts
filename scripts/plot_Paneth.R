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

setwd("/home/jiawen/myMLnet")
source("./codes/MLNconstruct.R")
source("./codes/Numericals.R")
source("./codes/Utilities.R")
source("./codes/Crosstalk_analysis.R")
source("./codes/Visualization.R")

# 0. Load the data and results

filen <- paste0("/home/jiawen/myMLnet/results/Paneth/result_ct7_ISC.csv")
CC_ISC <- read.csv(filen) 
filen <- "/home/jiawen/myMLnet/results/Paneth/result_ct7_Paneth_progenitor.csv"
CC_PC <- read.csv(filen)
CC_pair_ISC <- Aggregate_Causality(CC_ISC, sum_type = "sum",data_type = "Target")
CC_pair_PC <- Aggregate_Causality(CC_PC, sum_type = "sum",data_type = "Target")

KeyGene <- "Lgr5"
PCP_recs <- c("Ror1","Ror2","Ryk","Gpc4","Vangl1") %>% sort()
common_recs <- lapply(1:10, function(x){
  paste0("Fzd",x)
}) %>% unlist() %>% sort()
Canonical_recs <- c("Lrp5","Lrp6")

# ISC chord
Fid_mat <- Calculate_Fidelity_Matrix(CC_ISC,KeyGene)$all %>% t()
custom_colors <- function(from) {
  if (from %in% PCP_recs) {
    return("red")  
  } else if (from %in% common_recs) {
    return("lightskyblue") 
  } else if(from %in% Canonical_recs){
    return("springgreen") 
  }
}
PCP_recs_used <- intersect(PCP_recs, colnames(Fid_mat)) %>% sort()
common_recs_used <- intersect(common_recs, colnames(Fid_mat)) %>% sort()
Canonical_recs_used <- intersect(Canonical_recs, colnames(Fid_mat)) %>% sort()
recs_used <- c(PCP_recs_used,common_recs_used,Canonical_recs_used)

edge_colors <- mapply(
  custom_colors,
  from = rep(colnames(Fid_mat), times = nrow(Fid_mat))
)
plot_order <-  c(sort(rownames(Fid_mat)),recs_used)

cairo_pdf(file = "./results/Paneth/figb1.pdf",width = 12,height = 10)
PlotXT_Chord(Fid_mat, orders = plot_order, edge_colors = edge_colors)
dev.off()

# PC pro cells chord
Fid_mat <- Calculate_Fidelity_Matrix(CC_PC,KeyGene)$all %>% t()
custom_colors <- function(from) {
  if (from %in% PCP_recs) {
    return("red")  
  } else if (from %in% common_recs) {
    return("lightskyblue") 
  } else if(from %in% Canonical_recs){
    return("springgreen") 
  }
}
PCP_recs_used <- intersect(PCP_recs, colnames(Fid_mat)) %>% sort()
common_recs_used <- intersect(common_recs, colnames(Fid_mat)) %>% sort()
Canonical_recs_used <- intersect(Canonical_recs, colnames(Fid_mat)) %>% sort()
recs_used <- c(PCP_recs_used,common_recs_used,Canonical_recs_used)

edge_colors <- mapply(
  custom_colors,
  from = rep(colnames(Fid_mat), times = nrow(Fid_mat))
)
plot_order <-  c(sort(rownames(Fid_mat)),recs_used)

cairo_pdf(file = "./results/Paneth/figb2.pdf",width = 12,height = 10)
PlotXT_Chord(Fid_mat, orders = plot_order, edge_colors = edge_colors)
dev.off()


Fid_ISC <- Calculate_Pairwise(CC_ISC,KeyGene,type = "Fid")
# Fid_ISC <- rbind(Fid_ISC, data.frame(Receptor = "Ryk",Target = "Fzd2", Fidelity = 0) )
Fid_ISC$color <- case_when(
            (Fid_ISC$Receptor %in% common_recs) ~ "Frizzled receptor",  
            (Fid_ISC$Receptor %in% PCP_recs) ~ "PCP co-receptor",  
            (Fid_ISC$Receptor %in% Canonical_recs) ~ "Canonical co-receptor"   
            )
agg <- aggregate(Fid_ISC$Fidelity, list(Fid_ISC$color,Fid_ISC$Target), sum)
colnames(agg) <- c("Signal","Target","Fidelity")
agg$type <- "ISC"

Fid_PC <- Calculate_Pairwise(CC_PC,KeyGene,type = "Fid")
Fid_PC$color <- case_when(
            (Fid_PC$Receptor %in% common_recs) ~ "Frizzled receptor",  
            (Fid_PC$Receptor %in% PCP_recs) ~ "PCP co-receptor",  
            (Fid_PC$Receptor %in% Canonical_recs) ~ "Canonical co-receptor"   
            )
agg2 <- aggregate(Fid_PC$Fidelity, list(Fid_PC$color,Fid_PC$Target), sum)
colnames(agg2) <- c("Signal","Target","Fidelity")
agg2$type <- "PC"

agg <- rbind(agg,agg2)

ggplot(agg, aes(x = type, y = Fidelity, fill = Signal)) +
    geom_bar(stat = "identity", position = "dodge", color = "black") +
    labs(x = "Cell type", y = "Total fidelity", fill = "Signal") +
    scale_fill_manual(
        values = c("PCP co-receptor" = "red", 
        "Frizzled receptor" = "lightskyblue", 
        "Canonical co-receptor" = "springgreen")
    ) +
    theme_bw() +
    theme(text = element_text(family = "Arial")) +
    theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0, hjust = 0.1)) +
    theme(axis.text.y = element_text(size = 16, face = "bold")) +
    theme(axis.title.x = element_text(size = 24, face = "bold", vjust = 0.5)) +
    theme(axis.title.y = element_text(size = 24, face = "bold", vjust = 1)) +
    theme(legend.title = element_text(size = 24, face = "bold")) +
    theme(legend.text = element_text(size = 16)) +
    theme(legend.position = "top")
ggsave(filename = "./results/Paneth/figc.pdf",device = cairo_pdf,width = 2400,height = 3600,unit = "px")

Idents(subobj) <- subobj$proma_cell_type
VlnPlot(subobj,features = "Lgr5")
