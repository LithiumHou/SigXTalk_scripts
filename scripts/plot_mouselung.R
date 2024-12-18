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

current_type <- "Krt8 ADI"
all_timepoints <- list(
  c("day 2", "day 3"),
  c("day 10", "day 11"),
  c("day 15")
)
all_timepoints <- lapply(all_timepoints, function(x) paste(x, collapse = "+")) %>% unlist()

CCC_out <- list()
CCC_in <- list()
CC_list <- list()
CC_pair_list <- list()
Pathway_list <- list()

for (i in 1:length(all_timepoints)) {
  filename <- paste0("./datasets/GSE141259/Seurat_t",i,".rds")
  SeuratObj <- readRDS(filename)

  LigRec_original <- Infer_CCI(SeuratObj, cellchat_output = T, db_use = "mouse")
  LR_Pairprob_in <- Extract_LR_Prob(LigRec_original, target_type = current_type, cellchat_use = T,pv_threshold = 0.1)
  LR_Pairprob_out <- Extract_LR_Prob(LigRec_original, source_type = current_type,target_type = NULL, cellchat_use = T,pv_threshold = 0.1)

  filename2 <- paste0("./results/GSE141259/results_t", i, ".csv")
  # Calculate the causality
  CC_temp <- read.csv(filename2)

  CC_temp <- CC_temp[!grepl("^mt-", CC_temp$Target), ]
  CC_temp <- CC_temp[!grepl("^Rpl", CC_temp$Target), ]
  CC_temp <- CC_temp[!grepl("^Rps", CC_temp$Target), ]
  CC_temp <- CC_temp[!grepl("^Hsp", CC_temp$Target), ]
  CC_temp <- CC_temp[!grepl("^H2-", CC_temp$Target), ]


  CCC_in[[i]] <- LR_Pairprob_in
  CCC_out[[i]] <- LR_Pairprob_out
  CC_list[[i]] <- CC_temp
  CC_pair_list[[i]] <- Aggregate_Causality(CC_temp, sum_type = "sum")
  pathways <- apply(CC_temp, 1,function(x){
    paste0(x[1],'+',x[2],'+',x[3])
  })
  Pathway_list[[i]] <- data.frame(
    pathway = pathways,
    weight = CC_temp$Weight,
    type = all_timepoints[i]
  )
}
names(CC_list) <- names(Pathway_list) <- names(CC_pair_list) <- all_timepoints

# FIGURE A: Seurat object
filename <- paste0("./datasets/GSE141259/Seurat_t1.rds")
SeuratObj <- readRDS(filename)
p <- DimPlot(SeuratObj, reduction = "umap", pt.size = 2.05)
p + theme(text = element_text(family = "Arial"))+
theme(axis.title = element_text(size = 24)) +
  theme(axis.text = element_blank()) +
  theme(legend.text = element_text(size = 24))
ggsave(filename = "./results/GSE141259/figures/fig5a.pdf",device = cairo_pdf,width = 3200,height = 3000,unit = "px")

# FIGURE B: Comparison of numbers of LRIs
CCC_counts <- lapply(1:length(all_timepoints), function(x) {
  timepoint <- all_timepoints[x]
  c1 <- c(timepoint, "Senders", nrow(CCC_out[[x]]))
  c2 <- c(timepoint, "Receivers", nrow(CCC_in[[x]]))
  rbind(c1,c2)
})
CCC_counts <- do.call(rbind, CCC_counts) %>% as.data.frame()
colnames(CCC_counts) <- c("Time", "CCI_type", "Number")

CCC_counts$Number <- as.numeric(CCC_counts$Number)
CCC_counts$Time <- fct_relevel(CCC_counts$Time, all_timepoints)

ggplot(CCC_counts) +
  geom_bar(aes(x = Time, y = Number, fill = CCI_type), stat = "identity", position = "dodge") +
  theme_bw() +
  labs(x = "Days after processing", y = "Number of active LRIs", fill = "Krt8 ADI cells as") +
  theme(axis.title = element_text(size = 24)) +
  theme(axis.text = element_text(size = 18)) +
  theme(text = element_text(family = "Arial")) +
  theme(legend.position = "top") +
  theme(
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 18)
  )
ggsave(filename = "./results/GSE141259/figures/fig5b.pdf",device = cairo_pdf,width = 2400,height = 2000,unit = "px")

# FIGURE C: Comparison of number of crosstalk pathways
topk <- 3
all_results <- list()
Key_counts <- c()
for (i in 1:length(all_timepoints)) {
  CC_results <- CC_list[[i]] %>% as.data.frame()
  Counts_pathway <- Count_Crosstalk(CC_results, verbose = F, datatype = "SSC")
  No_tfs <- length(Counts_pathway)
  Counts_sorted <- sort(Counts_pathway)
  results <- data.frame(orders = 1:No_tfs, gene = names(Counts_sorted), pathways = Counts_sorted)
  results$pathways <- as.numeric(results$pathways)
  results$pathways[results$pathways == 0] <- 1
  results$type <- all_timepoints[i]
  all_results[[i]] <- results
  Key_counts[[i]] <- results[order(results$pathways, decreasing = T)[1:topk], c("gene", "pathways", "type")]
}

all_results <- do.call(rbind, all_results)
all_results$type <- as.factor(all_results$type)
all_results %>%
  mutate(type = fct_reorder(type, orders)) %>%
  ggplot(aes(x = orders, y = pathways, group = type, color = type)) +
  geom_line(linewidth = 2) +
  scale_fill_discrete(breaks = rev(levels(all_results$type))) +
  scale_color_viridis(discrete = TRUE) +
  theme_bw() +
  theme(text = element_text(family = "Arial")) +
  labs(x = "Rank of SSC", y = "Crosstalk pathways", color = "Days after \nprocessing") +
  theme(axis.title = element_text(size = 24)) +
  theme(axis.text = element_text(size = 16)) +
  theme(
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 16)
  ) +
  theme(
    legend.position = "right",
    legend.direction = "vertical"
  )
ggsave(filename = "./results/GSE141259/figures/fig5c1.pdf",device = cairo_pdf,width = 3400,height = 2800,unit = "px")

Key_counts <- do.call(rbind, Key_counts)
Key_counts$type <- as.factor(Key_counts$type)
Key_counts$type <- fct_relevel(Key_counts$type, all_timepoints)
Key_counts$order <- 1:nrow(Key_counts)
Key_counts %>%
  ggplot(aes(x = order, y = pathways, fill = type)) +
  geom_bar(stat = "identity", alpha = .6, width = .4) +
  scale_x_continuous(breaks = Key_counts$order, labels = Key_counts$gene) +
  labs(x = "SSC", y = "Crosstalk pathways", fill = "Days after \nprocessing") +
  theme_bw() +
  theme(text = element_text(family = "Arial")) +
  theme(axis.title = element_text(size = 24)) +
  theme(axis.text = element_text(size = 18)) +
  theme(
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 16)
  ) +
  theme(
    legend.position = "right",
    legend.direction = "vertical"
  )
ggsave(filename = "./results/GSE141259/figures/fig5c2.pdf",device = cairo_pdf,width = 3000,height = 2000,unit = "px")

# Density plot
names(CC_pair_list) <- all_timepoints
Fid_list <- lapply(all_timepoints,function(x){
  temp <- CC_pair_list[[x]]
  temp2 <- Calculate_Pairwise(temp, type = "Fid")
  temp2$Time <- x
  temp2
})
Fid_df <- do.call(rbind, Fid_list)
Pathway_df <- do.call(rbind,Pathway_list)
average_Fid <- aggregate(Fid_df$Fidelity, list(Fid_df$Time),mean)
colnames(average_Fid) <- c("Days_after_processing", "Average_Fidelity")
average_Fid$order <- c(2,3,1)

Fid_df$order <- match(Fid_df$Time,all_timepoints) %>% as.character()

Fid_df %>%
    mutate(Time = fct_reorder(Time, order)) %>%
    ggplot(aes(y=order, x=Fidelity,  fill=Time)) +
      geom_density_ridges(alpha=0.6) +
      theme_ridges() +
      theme_bw() +
      theme(
        legend.position="none",
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 8)
      ) +
      xlab("Fidelity") +
      ylab("Days after processing") +
      theme(text = element_text(family = "Arial")) +
      theme(axis.text.x = element_text(size = 20, angle = 90, face = "bold", vjust = 0.5)) +
      theme(axis.text.y = element_text(size = 20, face = "bold")) +
      theme(axis.title.x = element_text(size = 24, face = "bold", vjust = 0.5)) +
      theme(axis.title.y = element_text(size = 24, face = "bold", vjust = 1)) +
      theme(legend.title = element_text(size = 18, face = "bold", vjust = 1)) +
      theme(legend.text = element_text(size = 16)) 

ggsave(filename = "./results/GSE141259/figures/fig5d_outer.pdf",device = cairo_pdf,width = 3000,height = 3000,unit = "px")

average_Fid %>% 
  mutate(Days_after_processing = fct_reorder(Days_after_processing, order)) %>%
  ggplot(aes(x = Days_after_processing, y = Average_Fidelity, fill = Days_after_processing))+
  geom_bar(stat = "identity", alpha = .6, width = .4) +
  theme_bw()+
  theme(legend.position="none")+
  theme(text = element_text(family = "Arial")) +
  theme(axis.title = element_text(size = 24)) +
  theme(axis.text = element_text(size = 18)) +
  theme(
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )+
  labs(x = "Days after processing", y = "Average Fidelity", fill = NULL)
ggsave(filename = "./results/GSE141259/figures/fig5d_inner.pdf",device = cairo_pdf,width = 3000,height = 3000,unit = "px")

# FIGURE E: Rec-TFTG bundles in three stages
topk <- 20
plot_list <- list()
for (i in 1:length(all_timepoints)) {
  CC_temp <- CC_list[[i]] %>% as.data.frame()
  CC_temp <- CC_temp[order(CC_temp$Weight, decreasing = T), ]
  topk_new <- min(topk, nrow(CC_temp))
  CC_df <- data.frame(
    Rec = CC_temp$Receptor[1:topk_new], TFTG = paste0(CC_temp$SSC[1:topk_new], "--", CC_temp$Target[1:topk_new]),
    Weight_all = CC_temp$Weight[1:topk_new]
  )
  CC_df$Weight_all <- as.numeric(CC_df$Weight_all)

  filename <- paste0("./datasets/GSE141259/Seurat_t",i,".rds")
  SeuratObj <- readRDS(filename)
  Exp_clu <- Get_Exp_Clu(SeuratObj, clusterID = current_type, datatype = "counts")
  vertices_1 <- data.frame(name = unique(CC_df$Rec), type = "Signal")
  vertices_1$value <- sapply(vertices_1$name, function(x) {
    temp <- Exp_clu[x,]
    mean(temp)
  })

  vertices_2 <- data.frame(name = unique(CC_df$TFTG), type = "SSC-Target")
  vertices_2$value <- sapply(vertices_2$name, function(x) {
    TF <- strsplit(x,split = "--")[[1]][1]
    TG <- strsplit(x,split = "--")[[1]][2]
    temp1 <- Exp_clu[TF,]
    temp2 <- Exp_clu[TG,]

    mean(sqrt(temp1*temp2))
  })

  all_vertices <- rbind(vertices_1, vertices_2) %>% as.data.frame()
  edges_1 <- rbind(c("Origin", "Recs"), c("Origin", "TFTGs"))
  edges_2 <- rbind(
    cbind(rep("Recs", vertices_1 %>% nrow()), vertices_1$name),
    cbind(rep("TFTGs", vertices_2 %>% nrow()), vertices_2$name)
  )
  all_edge <- rbind(edges_1, edges_2)
  vertices_3 <- data.frame(name = c("Origin", "Recs", "TFTGs"), type = c(NA, "Origin", "Origin"), value = runif(3))
  all_vertices <- rbind(vertices_3, all_vertices)
  all_vertices$value <- as.numeric(all_vertices$value)*20
  mygraph <- graph_from_data_frame(all_edge, vertices = all_vertices)
  from <- match(CC_df$Rec, all_vertices$name)
  to <- match(CC_df$TFTG, all_vertices$name)

  p <- do.call(ggraph, list(mygraph, layout = "dendrogram", circular = TRUE)) +
    geom_node_point(aes(filter = leaf, x = x * 1.05, y = y * 1.05, colour = type, size = value)) +
    geom_node_text(aes(
      x = x * 1.3, y = y * 1.3, filter = leaf, label = name,
      colour = type,size = 20
    )) +
    geom_conn_bundle(
      data = get_con(from = from, to = to, Strength = CC_df$Weight),
      aes(color = Strength, edge_width = 1.),
      arrow = arrow(length = unit(6, "mm")), start_cap = circle(4, "mm"), end_cap = circle(6, "mm")
    ) +
    theme(text = element_text(family = "Arial")) +
    scale_edge_color_continuous(low = "yellow", high = "blue") +
    theme_void() +
    theme(legend.position = "right")

  p
  plotname <- paste0("./results/GSE141259/figures/fig5e_",i,".pdf")
  ggsave(filename = plotname,device = cairo_pdf,width = 5000,height = 4000,unit = "px")
}

TG_list <- lapply(CC_list, function(x){
  x$Target %>% unique() %>% sort()
})
Reduce(intersect, TG_list) %>% sort()

# Figure F
KeyTG <- c("Actb")
p1 <- PlotXT_FidSpe(CC_list[[1]], KeyTG, threshold = 0.01)
p2 <- PlotXT_FidSpe(CC_list[[2]], KeyTG, threshold = 0.5)
p3 <- PlotXT_FidSpe(CC_list[[3]], KeyTG, threshold = 0.5)
p1 | p2 | p3

ggsave(filename = "./results/GSE141259/figures/fig5f.pdf",device = cairo_pdf,width = 10800,height = 2400,unit = "px")

# Figure G
SSC_list <- lapply(all_timepoints, function(x){
  temp <- CC_list[[x]]
  temp2 <- aggregate(temp$Weight, list(temp$SSC),mean)
  colnames(temp2) <- c("SSC","value")
  temp2$Time <- x
  temp2
})
SSC_df <- do.call(rbind,SSC_list)
SSC_df$order <- match(SSC_df$Time,all_timepoints) %>% as.character()
SSC_df$SSC <- reorder(SSC_df$SSC, -SSC_df$order)

shared_SSC <- intersect(intersect(CC_list[[1]]$SSC,CC_list[[2]]$SSC),CC_list[[3]]$SSC)
SSC_filtered <- filter(SSC_df, SSC %in% shared_SSC)
SSC_filtered %>%
  mutate(Time = fct_reorder(Time,order)) %>% 
  ggplot() +
  geom_bar(aes(x = SSC, y = value, fill = Time), stat = "identity",position = "dodge") +
  theme_bw() +
  labs(x = " ", y = "Average PRS", fill = "Days after processing") +
  theme(axis.title = element_text(size = 24)) +
  theme(axis.text = element_text(size = 18)) +
  theme(text = element_text(family = "Arial")) +
  theme(legend.position = "top") +
  theme(
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 18)
  )
ggsave(filename = "./results/GSE141259/figures/fig5g.pdf",device = cairo_pdf,width = 3000,height = 2400,unit = "px")

# supp a: overview of dataset
filename <- paste0("./datasets/GSE141259/Seurat_t1.rds")
SeuratObj <- readRDS(filename)
p <- DimPlot(SeuratObj, reduction = "umap", pt.size = 2.05)
p + theme(text = element_text(family = "Arial"))+
theme(axis.title = element_text(size = 24)) +
  theme(axis.text = element_blank()) +
  theme(legend.text = element_text(size = 24))
ggsave(filename = "./results/GSE141259/figures/suppa1.pdf",device = cairo_pdf,width = 3200,height = 3000,unit = "px")

filename <- paste0("./datasets/GSE141259/Seurat_t2.rds")
SeuratObj <- readRDS(filename)
p <- DimPlot(SeuratObj, reduction = "umap", pt.size = 2.05)
p + theme(text = element_text(family = "Arial"))+
theme(axis.title = element_text(size = 24)) +
  theme(axis.text = element_blank()) +
  theme(legend.text = element_text(size = 24))
ggsave(filename = "./results/GSE141259/figures/suppa2.pdf",device = cairo_pdf,width = 3200,height = 3000,unit = "px")

filename <- paste0("./datasets/GSE141259/Seurat_t3.rds")
SeuratObj <- readRDS(filename)
p <- DimPlot(SeuratObj, reduction = "umap", pt.size = 2.05)
p + theme(text = element_text(family = "Arial"))+
theme(axis.title = element_text(size = 24)) +
  theme(axis.text = element_blank()) +
  theme(legend.text = element_text(size = 24))
ggsave(filename = "./results/GSE141259/figures/suppa3.pdf",device = cairo_pdf,width = 3200,height = 3000,unit = "px")

# supp b: CCC patterns at three times
pdf(file = "./results/GSE141259/figures/suppb1.pdf",width = 10,height = 10)
PlotCCI_CirclePlot(CCC_in[[1]], topk = 20)
dev.off()
pdf(file = "./results/GSE141259/figures/suppb2.pdf",width = 10,height = 10)
PlotCCI_CirclePlot(CCC_in[[2]], topk = 20)
dev.off()
pdf(file = "./results/GSE141259/figures/suppb3.pdf",width = 10,height = 10)
PlotCCI_CirclePlot(CCC_in[[3]], topk = 20)
dev.off()

# supp c: Jun as the SSC
KeyTF <- "Jun"
topk <- 20
df1 <- filter(CC_list[[1]], SSC==KeyTF)
df1 <- df1[order(df1$Weight,decreasing = T),]
df1 <- df1[1:topk,]
df1$pathway <- paste0(df1$Receptor,'-',df1$SSC,'-',df1$Target)
df1 <- mutate(df1, pathway = fct_reorder(pathway, Weight))
ggplot(df1, aes(x=pathway, y=Weight)) +
    geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
    coord_flip() +
    xlab("Pathway") +
    ylab("PRS")+
    theme_bw()+
    theme(axis.title = element_text(size = 24)) +
    theme(axis.text = element_text(size = 18)) +
    theme(legend.position = "none")
ggsave(filename = "./results/GSE141259/figures/suppc1.pdf",device = cairo_pdf,width = 3500,height = 3000,unit = "px")

df1 <- filter(CC_list[[2]], SSC==KeyTF)
df1 <- df1[order(df1$Weight,decreasing = T),]
df1 <- df1[1:topk,]
df1$pathway <- paste0(df1$Receptor,'-',df1$SSC,'-',df1$Target)
df1 <- mutate(df1, pathway = fct_reorder(pathway, Weight))
ggplot(df1, aes(x=pathway, y=Weight)) +
    geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
    coord_flip() +
    xlab("Pathway") +
    ylab("PRS")+
    theme_bw()+
    theme(axis.title = element_text(size = 24)) +
    theme(axis.text = element_text(size = 18)) +
    theme(legend.position = "none")
ggsave(filename = "./results/GSE141259/figures/suppc2.pdf",device = cairo_pdf,width = 3500,height = 3000,unit = "px")

df1 <- filter(CC_list[[3]], SSC==KeyTF)
df1 <- df1[order(df1$Weight,decreasing = T),]
df1 <- df1[1:topk,]
df1$pathway <- paste0(df1$Receptor,'-',df1$SSC,'-',df1$Target)
df1 <- mutate(df1, pathway = fct_reorder(pathway, Weight))
ggplot(df1, aes(x=pathway, y=Weight)) +
    geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
    coord_flip() +
    xlab("Pathway") +
    ylab("PRS")+
    theme_bw()+
    theme(axis.title = element_text(size = 24)) +
    theme(axis.text = element_text(size = 18)) +
    theme(legend.position = "none")
ggsave(filename = "./results/GSE141259/figures/suppc3.pdf",device = cairo_pdf,width = 3200,height = 3000,unit = "px")


# supp d: