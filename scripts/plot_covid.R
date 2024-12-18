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
Seurat_hc <- readRDS("./datasets/COVID/seurat_hc.rds")
Seurat_covid <- readRDS("./datasets/COVID/seurat_covid.rds")

types_used <- c("hc","patient")
target_type <- "Fibroblasts"

results_list <- lapply(types_used, function(x){
    filen <- paste0('./results/COVID/results_',x,'.csv')

    temp <- read.table(filen,sep = ',',header = T)
    temp <- temp[!grepl("^MT", temp$Target), ]
    temp <- temp[!grepl("^RPL", temp$Target), ]
    temp <- temp[!grepl("^RPS", temp$Target), ]
    filter(temp, Weight > 0.01*max(temp$Weight))
})
names(results_list) <- types_used
CC_hc <- results_list$hc
CC_covid <- results_list$patient
CC_pair_hc <- Aggregate_Causality(CC_hc, sum_type = "sum",data_type = "Target")
CC_pair_covid <- Aggregate_Causality(CC_covid, sum_type = "sum",data_type = "Target")

# 1. CCC patterns
cell_anno <- data.frame(cell = names(Idents(Seurat_covid)), cluster = Idents(Seurat_covid) %>% as.character())
CCC_covid <- Infer_CCI(Seurat_covid, LRDB = NULL, cellchat_output = T, use_spatial = F, db_use = "human")
CCC_threshold <- 0.05
LR_Pairprob <- Extract_LR_Prob(CCC_covid, target_type = target_type, cellchat_use = T)
LR_Pairprob <- LR_Pairprob[which(LR_Pairprob$Weight >= CCC_threshold * max(LR_Pairprob$Weight)), ]
LR_covid <- aggregate(LR_Pairprob$Weight,list(LR_Pairprob$To),sum)
colnames(LR_covid) <- c("Receptor","Weight")
LR_covid$Sample <- "COVID"

pdf(file = "./results/COVID/figures_new/fig4a.pdf",width = 10,height = 10)
PlotCCI_CirclePlot(LR_Pairprob, topk = 20)
dev.off()

cell_anno <- data.frame(cell = names(Idents(Seurat_hc)), cluster = Idents(Seurat_hc) %>% as.character())
CCC_hc <- Infer_CCI(Seurat_hc, LRDB = NULL, cellchat_output = T, use_spatial = F, db_use = "human")
CCC_threshold <- 0.05
LR_Pairprob <- Extract_LR_Prob(CCC_hc, target_type = target_type, cellchat_use = T)
LR_Pairprob <- LR_Pairprob[which(LR_Pairprob$Weight >= CCC_threshold * max(LR_Pairprob$Weight)), ]
LR_hc <- aggregate(LR_Pairprob$Weight,list(LR_Pairprob$To),sum)
colnames(LR_hc) <- c("Receptor","Weight")
LR_hc$Sample <- "HC"

pdf(file = "./results/COVID/figures_new/fig4b.pdf",width = 10,height = 10)
PlotCCI_CirclePlot(LR_Pairprob, topk = 20)
dev.off()

# 2. CCC comparison
LR_all <- rbind(LR_hc,LR_covid) %>% as.data.frame()
LR_all$Weight <- as.numeric(LR_all$Weight)

LR_all <- Complete_groups(LR_all)
LR_all$order <- (LR_all$Sample == "COVID")*LR_all$Weight-(LR_all$Sample == "HC")*LR_all$Weight
LR_all$Receptor <- reorder(LR_all$Receptor, -LR_all$order)
ggplot(LR_all) +
  geom_bar(aes(x = Receptor, y = Weight, fill = Sample), stat = "identity",position = "dodge") +
  theme_bw() +
  labs(x = "CCC Signal", y = "CCC strength", fill = "Receptors activated in sample") +
  theme(axis.title = element_text(size = 24)) +
  theme(axis.text = element_text(size = 18)) +
  theme(text = element_text(family = "Arial")) +
  theme(legend.position = "top") +
  theme(
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 18)
  )
ggsave(filename = "./results/COVID/figures_new/fig4c.pdf",device = cairo_pdf,width = 3600,height = 2400,unit = "px")

# 3. common receptor

common_rec <- intersect(CC_hc$Receptor,CC_covid$Receptor)

count_list <- lapply(common_rec, function(x){
  temp_hc <- filter(CC_hc, Receptor == x)
  temp_covid <- filter(CC_covid, Receptor == x)
  
  # count of each target
  temp_hc_count <- aggregate(temp_hc$Weight, list(temp_hc$Target),length)
  temp_hc_count$Sample <- "HC"
  temp_covid_count <- aggregate(temp_covid$Weight, list(temp_covid$Target),length)
  temp_covid_count$Sample <- "COVID"
  df <- rbind(temp_covid_count, temp_hc_count)
  colnames(df) <- c("Target", "Count", "Sample")
  df$Receptor <- x
  df
})
count_df <- do.call(rbind, count_list)
ggplot(count_df,aes(x = Receptor, y = Count, fill = Sample)) +
    geom_boxplot() +
    theme_bw()+
    theme(legend.position = "top") +
    labs(x = "Signal", y = "Crosstalk pathways per target") +
    theme(text = element_text(family = "Arial")) +
    theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0, hjust = 0.1)) +
    theme(axis.text.y = element_text(size = 16, face = "bold")) +
    theme(axis.title.x = element_text(size = 24, face = "bold", vjust = 0.5)) +
    theme(axis.title.y = element_text(size = 24, face = "bold", vjust = 1)) +
    theme(legend.title = element_text(size = 24, face = "bold")) +
    theme(legend.text = element_text(size = 16)) +
    theme(legend.position = "top")
ggsave(filename = "./results/COVID/figures_new/fig4d.pdf",device = cairo_pdf,width = 1600,height = 3600,unit = "px")

# 4. specificity distribution - common Receptors
Spe_list <- lapply(common_rec, function(x){
  spe_hc <- Calculate_Specificity_Matrix(CC_hc, KeyRec = x)$all %>% t()
  spe_hc <- mat_df(spe_hc)
  colnames(spe_hc) <- c("SSC", "Target", "Specificity")
  spe_hc <- filter(spe_hc, Specificity>0)
  spe_hc$Sample <- "HC"

  spe_covid <- Calculate_Specificity_Matrix(CC_covid, KeyRec = x)$all %>% t()
  spe_covid <- mat_df(spe_covid)
  colnames(spe_covid) <- c("SSC", "Target", "Specificity")
  spe_covid <- filter(spe_covid, Specificity>0)
  spe_covid$Sample <- "COVID"
  temp <- rbind(spe_hc[,2:4],spe_covid[,2:4])
  temp$Receptor <- x
  temp
})

df <- do.call(rbind, Spe_list)
ggplot(df, aes(x = Receptor, y = Specificity, fill = Sample)) +
    geom_boxplot(alpha = 0.3) +
    theme_bw()+
    theme(legend.position = "none") +
    labs(x = "Signal", y = "Specificity") +
    theme(text = element_text(family = "Arial")) +
    theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0, hjust = 0.1)) +
    theme(axis.text.y = element_text(size = 16, face = "bold")) +
    theme(axis.title.x = element_text(size = 24, face = "bold", vjust = 0.5)) +
    theme(axis.title.y = element_text(size = 24, face = "bold", vjust = 1)) +
    theme(legend.title = element_text(size = 24, face = "bold")) +
    theme(legend.text = element_text(size = 16)) +
    theme(legend.position = "top")
ggsave(filename = "./results/COVID/figures_new/fig4e.pdf",device = cairo_pdf,width = 1600,height = 3600,unit = "px")

# 5. Specificity distribution - other
rec_hc <- setdiff(CC_hc$Receptor, common_rec)
rec_covid <- setdiff(CC_covid$Receptor, common_rec)

Spe_df_hc <- Calculate_Pairwise(CC_pair_hc, KeyGene = rec_hc, type = "Spe")
Spe_df_covid <- Calculate_Pairwise(CC_pair_covid, KeyGene = rec_covid, type = "Spe")
Spe_df <- rbind(Spe_df_covid, Spe_df_hc)
PlotXT_MultiCircularBar(df = Spe_df[c('Target','Receptor','Specificity')], KeyFactors = c(rec_hc,rec_covid))
ggsave(filename = "./results/COVID/figures_new/fig4f.pdf",device = cairo_pdf,width = 3600,height = 3600,unit = "px")

# 6. common targets
TG_used <- FindMarkers(Seurat_covid, target_type, min.pct = 0.25, only.pos = T, logfc.threshold = 0.25)
TG_used1 <- filter(TG_used, p_val_adj<1e-3) %>% rownames()
TG_used <- FindMarkers(Seurat_hc, target_type, min.pct = 0.25, only.pos = T, logfc.threshold = 0.25)
TG_used2 <- filter(TG_used, p_val_adj<1e-3) %>% rownames()
TGs <- setdiff(TG_used1,TG_used2) %>% sort()
KeyGene <- "STAT3"

# STAT3 as an SSC
df_filtered <- filter(CC_covid, SSC == KeyGene) %>% filter(Target %in% TGs)
df_filtered <- df_filtered[,-2]
temp_mat <- df_mat(df_filtered, row = Receptor, col = Target, value = Weight)
temp_mat[is.na(temp_mat)] <- 0
# temp_mat <- temp_mat/sum(temp_mat)
temp_mat2 <- temp_mat[order(apply(temp_mat/sum(temp_mat),1,mean),decreasing = T),order(apply((temp_mat)/sum(temp_mat),2,mean),decreasing = T)]
temp_mat2 <- temp_mat2[,1:20]
# row (receptor) annotation
ra_values <- sweep(temp_mat2, 2, colSums(temp_mat2), FUN = "/") # mat/colsum
ra = rowAnnotation(Spe = anno_boxplot(ra_values), height = unit(4, "cm"))
# column (Target) annotation
ha_values <- sweep(temp_mat2, 1, rowSums(temp_mat2), FUN = "/")
ha = HeatmapAnnotation(Fid = anno_boxplot(ha_values), which='column',height = unit(4, "cm"))

pdf(file = "./results/COVID/figures_new/fig4g.pdf",width = 10,height = 10)
ComplexHeatmap::Heatmap(temp_mat2,right_annotation = ra, top_annotation = ha,
                        cluster_rows = F, cluster_columns = F,name ="PRS",
                        row_names_gp = gpar(fontsize = 18),
                        column_names_gp = gpar(fontsize = 18),
                        width = unit(18, "cm"), height = unit(4, "cm"))
dev.off()

# STAT3 as a target
df_filtered <- filter(CC_covid, Target == KeyGene) %>% filter(Weight > 0.5)
df_filtered$Receptor <- factor(df_filtered$Receptor, levels = c(rec_covid,common_rec))
df_filtered$SSC <- reorder(df_filtered$SSC, df_filtered$Weight)

p <- PlotXT_Sankey(df_filtered, min_weight = 0)
saveWidget(p, file="./results/COVID/figures_new/fig4h.html")

pdf(file = "./results/COVID/figures_new/fig4h.pdf",width = 12,height = 10)
alluvial::alluvial(df_filtered[,1:3] ,
         freq=df_filtered$Weight, border=NA, alpha = 0.5,
         col = ifelse(df_filtered$Receptor %in% rec_covid, "red", "lightskyblue"),
         cex=1,
         axis_labels = c("Signal", "SSC", "Target")
         )

dev.off()

# compare the common TF
tf_hc <- CC_hc$SSC %>% unique()
tf_covid <- CC_covid$SSC %>% unique()
common_tf <- intersect(tf_hc,tf_covid)
CC_hc_filtered <- CC_hc
CC_covid_filtered <- CC_covid
common_target <- intersect(CC_covid_filtered$Target,CC_hc_filtered$Target)
setdiff(tf_covid,tf_hc)

CC_hc_filtered$TFTG <- paste0(CC_hc_filtered$SSC,"-",CC_hc_filtered$Target)
CC_covid_filtered$TFTG <- paste0(CC_covid_filtered$SSC,"-",CC_covid_filtered$Target)

CC_hc_filtered <- CC_hc_filtered %>%
  group_by(Receptor) %>%
  mutate(Spe = Weight / sum(Weight))
CC_hc_filtered$sample <- "HC"
spe_hc_agg <- aggregate(CC_hc_filtered$Spe, list(CC_hc_filtered$TFTG), sum)
spe_hc_agg$sample <- "HC"
rownames(spe_hc_agg) <- spe_hc_agg$Group.1

CC_covid_filtered <- CC_covid_filtered %>%
  group_by(Receptor) %>%
  mutate(Spe = Weight / sum(Weight))
CC_covid_filtered$sample <- "COVID"
spe_covid_agg <- aggregate(CC_covid_filtered$Spe, list(CC_covid_filtered$TFTG), sum)
spe_covid_agg$sample <- "HC"
rownames(spe_covid_agg) <- spe_covid_agg$Group.1

common_tftg <- intersect(spe_covid_agg$Group.1, spe_hc_agg$Group.1)
df <- data.frame(TFTG = common_tftg, COVID = spe_covid_agg[common_tftg, 2], HC = spe_hc_agg[common_tftg, 2])
sscs <- lapply(common_tftg, function(x){
  temp <- strsplit(x,'-')
  unlist(temp)[1]
}) %>% unlist()
df$SSC <- sscs

ggplot(df, aes(x=COVID, y=HC, color=SSC)) + 
    geom_point(size=4) +
    xlim(0,0.0085)+
    ylim(0,0.0085)+
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    theme_bw()+
    theme(text = element_text(family = "Arial")) +
    theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0, hjust = 0.1)) +
    theme(axis.text.y = element_text(size = 16, face = "bold")) +
    theme(axis.title.x = element_text(size = 24, face = "bold", vjust = 0.5)) +
    theme(axis.title.y = element_text(size = 24, face = "bold", vjust = 1)) +
    theme(legend.title = element_text(size = 24, face = "bold")) +
    theme(legend.text = element_text(size = 16)) +
    theme(legend.position = "right")

# common TG
common_target <- intersect(CC_pair_hc$Target,CC_pair_covid$Target)
CC_hc_filtered <- filter(CC_pair_hc, Target %in% common_target)
CC_covid_filtered <- filter(CC_pair_covid, Target %in% common_target)
CC_hc_filtered <- CC_hc_filtered %>%
  group_by(Target) %>%
  mutate(Fid = Weight / sum(Weight)) %>% filter(Receptor %in% common_rec)
CC_hc_filtered$sample <- "HC"

CC_covid_filtered <- CC_covid_filtered %>%
  group_by(Target) %>%
  mutate(Fid = Weight / sum(Weight)) %>% filter(Receptor %in% common_rec)
CC_covid_filtered$sample <- "COVID"
df <- rbind(CC_hc_filtered[,-3],CC_covid_filtered[,-3])
df1 <- filter(df, Receptor == "CD44")
df2 <- filter(df, Receptor == "PTPRM")

ggplot(data=df1, aes(x=Fid, group=sample, fill=sample)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_ipsum()+
    xlim(0,1) +
    labs(x = "Fidelity", y = "Density")+
    theme(
      axis.title.x = element_blank(),     # Remove x-axis title
      axis.ticks.x = element_blank(),     # Remove x-axis ticks
      axis.text.x = element_blank(),      # Remove x-axis tick labels
      axis.line.x = element_blank()       # Remove x-axis line
    )  +
    theme(text = element_text(family = "Arial")) +
    theme(axis.text.y = element_text(size = 16, face = "bold")) +
    theme(axis.title.y = element_text(size = 24, face = "bold", vjust = 1)) +
    theme(legend.title = element_text(size = 24, face = "bold")) +
    theme(legend.text = element_text(size = 16)) +
    theme(legend.position = "top")
ggsave(filename = "./results/COVID/figures_new/fig4i1.pdf",device = cairo_pdf,width = 3600,height = 1200,unit = "px")

ggplot(data=df2, aes(x=Fid, group=sample, fill=sample)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_ipsum()+
    xlim(0,1) +
    labs(x = "Fidelity", y = "Density")+
    theme(text = element_text(family = "Arial")) +
    theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0, hjust = 0.1)) +
    theme(axis.text.y = element_text(size = 16, face = "bold")) +
    theme(axis.title.x = element_text(size = 24, face = "bold", vjust = 0.5)) +
    theme(axis.title.y = element_text(size = 24, face = "bold", vjust = 1)) +
    theme(legend.title = element_text(size = 24, face = "bold")) +
    theme(legend.text = element_text(size = 16)) +
    theme(legend.position = "top")
ggsave(filename = "./results/COVID/figures_new/fig4i2.pdf",device = cairo_pdf,width = 3600,height = 1200,unit = "px")

# supp 1: overview of the data
DimPlot(Seurat_covid, reduction = "umap",pt.size = 1.5)
ggsave(filename = "./results/COVID/figures_new/suppa_covid.pdf",device = cairo_pdf,width = 3600,height = 3600,unit = "px")
DimPlot(Seurat_hc, reduction = "umap",pt.size = 1.5)
ggsave(filename = "./results/COVID/figures_new/suppa_hc.pdf",device = cairo_pdf,width = 3600,height = 3600,unit = "px")

# supp 2: collagen and laminin signaling
pathways.show <- c("COLLAGEN")
pdf(file = "./results/COVID/figures_new/suppb_covid.pdf",width = 10,height = 10)
netVisual_heatmap(CCC_covid, signaling = pathways.show,font.size = 20, font.size.title = 24, color.heatmap = "Reds")
dev.off()
pdf(file = "./results/COVID/figures_new/suppb_hc.pdf",width = 10,height = 10)
netVisual_heatmap(CCC_hc, signaling = pathways.show,font.size = 20, font.size.title = 24, color.heatmap = "Reds")
dev.off()
pathways.show <- c("LAMININ")
pdf(file = "./results/COVID/figures_new/suppc_covid.pdf",width = 10,height = 10)
netVisual_heatmap(CCC_covid, signaling = pathways.show,font.size = 20, font.size.title = 24, color.heatmap = "Blues")
dev.off()
pdf(file = "./results/COVID/figures_new/suppc_hc.pdf",width = 10,height = 10)
netVisual_heatmap(CCC_hc, signaling = pathways.show,font.size = 20, font.size.title = 24, color.heatmap = "Blues")
dev.off()

# supp 3: PRS comparison

CC_copy_COVID <- data.frame(value = CC_covid$Weight, Signal = "All", Sample = "COVID")
CC_copy_HC <- data.frame(value = CC_hc$Weight, Signal = "All", Sample = "HC")
CC_CD44_COVID <- filter(CC_covid, Receptor == "CD44")
CC_CD44_COVID <- data.frame(value = CC_CD44_COVID$Weight, Signal = "CD44", Sample = "COVID")
CC_CD44_HC <- filter(CC_hc, Receptor == "CD44")
CC_CD44_HC <- data.frame(value = CC_CD44_HC$Weight, Signal = "CD44", Sample = "HC")

CC_PTPRM_COVID <- filter(CC_covid, Receptor == "PTPRM")
CC_PTPRM_COVID <- data.frame(value = CC_PTPRM_COVID$Weight, Signal = "PTPRM", Sample = "COVID")
CC_PTPRM_HC <- filter(CC_hc, Receptor == "PTPRM")
CC_PTPRM_HC <- data.frame(value = CC_PTPRM_HC$Weight, Signal = "PTPRM", Sample = "HC")

df <- rbind(CC_CD44_COVID,CC_CD44_HC,CC_copy_COVID,CC_copy_HC,CC_PTPRM_COVID,CC_PTPRM_HC)
df[which(df$value > 3), "value"] = 3
ggplot(df,aes(x = Signal, y = value, fill = Sample)) +
    geom_boxplot() +
    theme_bw()+
    theme(legend.position = "top") +
    labs(x = "Signal", y = "PRS") +
    # theme(text = element_text(family = "Arial")) +
    theme(axis.text.x = element_text(size = 16, face = "bold", angle = 0, hjust = 0.1)) +
    theme(axis.text.y = element_text(size = 16, face = "bold")) +
    theme(axis.title.x = element_text(size = 24, face = "bold", vjust = 0.5)) +
    theme(axis.title.y = element_text(size = 24, face = "bold", vjust = 1)) +
    theme(legend.title = element_text(size = 24, face = "bold")) +
    theme(legend.text = element_text(size = 16)) +
    theme(legend.position = "top")
ggsave(filename = "./results/COVID/figures_new/suppd.pdf",device = cairo_pdf,width = 2400,height = 3600,unit = "px")

t.test(CC_copy_COVID$value,CC_copy_HC$value)
t.test(CC_CD44_COVID$value,CC_CD44_HC$value)
t.test(CC_PTPRM_COVID$value,CC_PTPRM_HC$value)

# supp 4: ITGB1's main targets:
Keygene <- "ITGB1"
pdf(file = "./results/COVID/figures_new/suppe.pdf",width = 10,height = 10)
PlotXT_HeatMap(CC_covid, Keygene, "Receptor", topk = 20)
dev.off()

# 