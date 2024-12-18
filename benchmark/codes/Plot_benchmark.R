library(ggplot2)
library(tidyverse)
library(patchwork) # To display 2 charts together
library(hrbrthemes)

# Plot the comparison
result_list <- lapply(1:10, function(x){
    filen <- paste0('/home/jiawen/myMLnet/benchmark/Beeline/outputs/comparison/comparison_',x,'.csv')
    temp <- read.csv(filen)
    temp
})

df <- do.call(rbind, result_list)
write.table(df,file = "/home/jiawen/myMLnet/benchmark/Beeline/outputs/comparison_all.csv", quote = F, sep = ",",row.names = F)
df$method %>% unique()
allm <- c("SigXTalk", "GeneLink","PIDC", "GRNVBEM", "GRNBOOST2", "PPCOR", "SCODE", "SINCERITIES", 
"LEAP", "GRISLI", "SINGE", "SCRIBE", "SCSGL")

df %>% 
    arrange(AUC) %>%
    mutate(name = factor(method, levels=allm)) %>%
    ggplot( aes(x=name, y=AUC, fill=name)) + 
    geom_boxplot() +
    xlab("Method") +
    ylab("AUROC") +
    scale_y_continuous(limits = c(0.4,0.9))+
    theme_bw() +
    theme(text = element_text(family = "Arial")) +
    theme(axis.title = element_text(size = 24)) +
    theme(axis.text = element_text(size = 24)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position="none") 

ggsave("/home/jiawen/myMLnet/benchmark/Beeline/outputs/results_comparison_AUROC.pdf",device = cairo_pdf,width = 3600,height = 2000,unit = "px")

df %>% 
    arrange(AUPR) %>%
    mutate(name = factor(method, levels=allm)) %>%
    ggplot( aes(x=name, y=AUPR, fill=name)) + 
    geom_boxplot() +
    xlab("Method") +
    ylab("AUPRC") +
    scale_y_continuous(limits = c(0.4,0.9))+
    theme_bw() +
    theme(text = element_text(family = "Arial")) +
    theme(axis.title = element_text(size = 24)) +
    theme(axis.text = element_text(size = 24)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position="none") 
ggsave("/home/jiawen/myMLnet/benchmark/Beeline/outputs/results_comparison_AUPRC.pdf",device = cairo_pdf,width = 3600,height = 2000,unit = "px")

# Cells
df <- read.csv("/home/jiawen/myMLnet/benchmark/Beeline/outputs/results_cells.csv")
df1 <- data.frame(Round = df$Round, Metric = "AUROC", Cells = df$Cells, value = df$AUROC)
df2 <- data.frame(Round = df$Round, Metric = "AUPRC", Cells = df$Cells, value = df$AUPRC)
df <- rbind(df1,df2)
allm <- df$Cells %>% unique() %>% sort()
df %>% 
    mutate(name = factor(Cells, levels=allm)) %>%
    mutate(Metric = factor(Metric, levels = c("AUROC","AUPRC"))) %>% 
    ggplot( aes(x=name, y=value, fill=Metric)) + 
    geom_boxplot() +
    xlab("Number of cells") +
    ylab("AUROC/AUPRC") +
    scale_y_continuous(limits = c(0.6,0.9))+
    theme_bw() +
    theme(text = element_text(family = "A")) +
    theme(axis.title = element_text(size = 24)) +
    theme(axis.text = element_text(size = 24)) +
    theme(legend.text = element_text(size = 24))+
    theme(legend.title = element_text(size = 24)) +
    theme(legend.position="top") 
ggsave("/home/jiawen/myMLnet/benchmark/Beeline/outputs/results_cells.pdf",device = cairo_pdf,width = 3000,height = 2000,unit = "px")

# Genes
df <- read.csv("/home/jiawen/myMLnet/benchmark/Beeline/outputs/results_genes.csv")
df1 <- data.frame(Round = df$Round, Metric = "AUROC", Genes = df$Genes, value = df$AUROC)
df2 <- data.frame(Round = df$Round, Metric = "AUPRC", Genes = df$Genes, value = df$AUPRC)
df <- rbind(df1,df2)
allm <- df$Genes %>% unique() %>% sort()
df %>% 
    mutate(name = factor(Genes, levels=allm)) %>%
    mutate(Metric = factor(Metric, levels = c("AUROC","AUPRC"))) %>% 
    ggplot( aes(x=name, y=value, fill=Metric)) + 
    geom_boxplot() +
    xlab("Number of genes") +
    ylab("AUROC/AUPRC") +
    scale_y_continuous(limits = c(0.6,0.9))+
    theme_bw() +
    theme(text = element_text(family = "A")) +
    theme(axis.title = element_text(size = 24)) +
    theme(axis.text = element_text(size = 24)) +
    theme(legend.text = element_text(size = 24))+
    theme(legend.title = element_text(size = 24)) +
    theme(legend.position="top") 
ggsave("/home/jiawen/myMLnet/benchmark/Beeline/outputs/results_genes.pdf",device = cairo_pdf,width = 3000,height = 2000,unit = "px")

# Learning rate
df <- read.csv("/home/jiawen/myMLnet/benchmark/Beeline/outputs/results_lr.csv")
df <- filter(df,batch==64)
df1 <- data.frame(Round = df$Round, Metric = "AUROC", LR = df$LR, value = df$AUROC)
df2 <- data.frame(Round = df$Round, Metric = "AUPRC", LR = df$LR, value = df$AUPRC)
df <- rbind(df1,df2)
allm <- df$LR %>% unique() %>% sort()
df %>% 
    mutate(name = factor(LR, levels=allm)) %>%
    mutate(Metric = factor(Metric, levels = c("AUROC","AUPRC"))) %>% 
    ggplot( aes(x=name, y=value, fill=Metric)) + 
    geom_boxplot() +
    xlab("Initial learning rate") +
    ylab("AUROC/AUPRC") +
    scale_y_continuous(limits = c(0.6,0.9))+
    theme_bw() +
    theme(text = element_text(family = "A")) +
    theme(axis.title = element_text(size = 24)) +
    theme(axis.text = element_text(size = 24)) +
    theme(legend.text = element_text(size = 24))+
    theme(legend.title = element_text(size = 24)) +
    theme(legend.position="top") 
ggsave("/home/jiawen/myMLnet/benchmark/Beeline/outputs/results_lr.pdf",device = cairo_pdf,width = 3000,height = 2000,unit = "px")

df <- read.csv("/home/jiawen/myMLnet/benchmark/Beeline/outputs/results_lr.csv")
df <- filter(df,LR == 0.01)

df1 <- data.frame(Round = df$Round, Metric = "AUROC", BS = df$batch, value = df$AUROC)
df2 <- data.frame(Round = df$Round, Metric = "AUPRC", BS = df$batch, value = df$AUPRC)
df <- rbind(df1,df2)
allm <- df$BS %>% unique() %>% sort()
df %>% 
    mutate(name = factor(BS, levels=allm)) %>%
    mutate(Metric = factor(Metric, levels = c("AUROC","AUPRC"))) %>% 
    ggplot( aes(x=name, y=value, fill=Metric)) + 
    geom_boxplot() +
    xlab("Batch size") +
    ylab("AUROC/AUPRC") +
    scale_y_continuous(limits = c(0.6,0.9))+
    theme_bw() +
    theme(text = element_text(family = "A")) +
    theme(axis.title = element_text(size = 24)) +
    theme(axis.text = element_text(size = 24)) +
    theme(legend.text = element_text(size = 24))+
    theme(legend.title = element_text(size = 24)) +
    theme(legend.position="top") 
ggsave("/home/jiawen/myMLnet/benchmark/Beeline/outputs/results_bs.pdf",device = cairo_pdf,width = 3000,height = 2000,unit = "px")

# Plot netdensity
df <- read.csv("/home/jiawen/myMLnet/benchmark/Beeline/outputs/results_netdense_1.csv")
df1 <- data.frame(Round = df$Round, Metric = "AUROC", Density = df$Density, value = df$AUROC)
df2 <- data.frame(Round = df$Round, Metric = "AUPRC", Density = df$Density, value = df$AUPRC)
df <- rbind(df1,df2)
allm <- df$Density %>% unique() %>% sort()
df %>% 
    mutate(name = factor(Density, levels=allm)) %>%
    mutate(Metric = factor(Metric, levels = c("AUROC","AUPRC"))) %>% 
    ggplot( aes(x=name, y=value, fill=Metric)) + 
    geom_boxplot() +
    xlab("Network wiring probability") +
    ylab("AUROC/AUPRC") +
    scale_y_continuous(limits = c(0.6,0.9))+
    theme_bw() +
    theme(text = element_text(family = "A")) +
    theme(axis.title = element_text(size = 24)) +
    theme(axis.text = element_text(size = 24)) +
    theme(legend.text = element_text(size = 24))+
    theme(legend.title = element_text(size = 24)) +
    theme(legend.position="top") 
ggsave("/home/jiawen/myMLnet/benchmark/Beeline/outputs/results_Density.pdf",device = cairo_pdf,width = 3000,height = 2000,unit = "px")

# Plot the train size
df <- read.csv("/home/jiawen/myMLnet/benchmark/Beeline/outputs/results_trainsize.csv")
df1 <- data.frame(Round = df$Round, Metric = "AUROC", ts = df$ts, value = df$AUROC)
df2 <- data.frame(Round = df$Round, Metric = "AUPRC", ts = df$ts, value = df$AUPRC)
df <- rbind(df1,df2)
allm <- df$ts %>% unique() %>% sort()
df %>% 
    mutate(name = factor(ts, levels=allm)) %>%
    mutate(Metric = factor(Metric, levels = c("AUROC","AUPRC"))) %>% 
    ggplot( aes(x=name, y=value, fill=Metric)) + 
    geom_boxplot() +
    xlab("Fraction of training sets") +
    ylab("AUROC/AUPRC") +
    scale_y_continuous(limits = c(0.6,0.9))+
    theme_bw() +
    theme(text = element_text(family = "A")) +
    theme(axis.title = element_text(size = 24)) +
    theme(axis.text = element_text(size = 24)) +
    theme(legend.text = element_text(size = 24))+
    theme(legend.title = element_text(size = 24)) +
    theme(legend.position="top") 
ggsave("/home/jiawen/myMLnet/benchmark/Beeline/outputs/results_trainsize.pdf",device = cairo_pdf,width = 3000,height = 2000,unit = "px")

# Plot the curve
df <- read.csv("/home/jiawen/myMLnet/benchmark/Beeline/outputs/results_curve.csv")
temperatureColor <- "blue"
priceColor <- "red"
ggplot(df, aes(x=Epoch)) +
    geom_line( aes(y=train_loss), size=2, color=temperatureColor) + 
    geom_line( aes(y=validation_loss*50), size=2, color=priceColor) +
    # Custom the Y scales:
    scale_y_continuous(
    # Features of the first axis
    name = "Training loss",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./50, name="Validation loss")
    ) +
    theme_bw() +
    theme(text = element_text(family = "A")) +
    theme(axis.title = element_text(size = 24)) +
    theme(axis.text = element_text(size = 24)) +
    theme(legend.text = element_text(size = 24))+
    theme(legend.title = element_text(size = 24)) +
    theme(legend.position="top") +
    geom_rect(aes(xmin=140,xmax=148,ymin=59.75,ymax=60.25),
            fill=temperatureColor,color=temperatureColor)+
    geom_rect(aes(xmin=140,xmax=148,ymin=57.75,ymax=58.25),
            fill=priceColor,color=priceColor)+
    annotate(geom='text',x=168,y=60,label='Training loss',size=8)+
    annotate(geom='text',x=170,y=58,label='Validation loss',size=8)

ggsave("/home/jiawen/myMLnet/benchmark/Beeline/outputs/results_curve.pdf",device = cairo_pdf,width = 6000,height = 2000,unit = "px")

# The random seed
df <- read.csv("/home/jiawen/myMLnet/benchmark/Beeline/outputs/results_randomseed.csv")
df1 <- data.frame(Round = df$Round, Metric = "AUROC", value = df$AUROC)
df2 <- data.frame(Round = df$Round, Metric = "AUPRC", value = df$AUPRC)
df <- rbind(df1,df2)
df %>% 
    ggplot( aes(x=Metric, y=value, fill=Metric)) + 
    geom_boxplot() +
    xlab(" ") +
    ylab("AUROC/AUPRC") +
    scale_y_continuous(limits = c(0.6,0.9))+
    theme_bw() +
    theme(text = element_text(family = "A")) +
    theme(axis.title = element_text(size = 24)) +
    theme(axis.text = element_text(size = 24)) +
    theme(legend.text = element_text(size = 24))+
    theme(legend.title = element_text(size = 24)) +
    theme(legend.position="top") 
ggsave("/home/jiawen/myMLnet/benchmark/Beeline/outputs/results_seed.pdf",device = cairo_pdf,width = 750,height = 2000,unit = "px")

# self-supervised
df <- read.csv("/home/jiawen/myMLnet/benchmark/Beeline/outputs/results_PPR.csv")
df$ts <- round(df$ts,digits = 1)
df1 <- data.frame(Round = df$Round, Metric = "AUROC", ts = df$ts, value = df$AUROC)
df2 <- data.frame(Round = df$Round, Metric = "AUPRC", ts = df$ts, value = df$AUPRC)
df <- rbind(df1,df2)
df <- filter(df, ts>0.1)
df <- filter(df, ts<0.8)

allm <- df$ts %>% unique() %>% sort()
df %>% 
    mutate(name = factor(ts, levels=allm)) %>%
    mutate(Metric = factor(Metric, levels = c("AUROC","AUPRC"))) %>% 
    ggplot( aes(x=name, y=value, fill=Metric)) + 
    geom_boxplot() +
    xlab("Fraction of training sets") +
    ylab("AUROC/AUPRC") +
    scale_y_continuous(limits = c(0.6,0.9))+
    theme_bw() +
    theme(text = element_text(family = "A")) +
    theme(axis.title = element_text(size = 24)) +
    theme(axis.text = element_text(size = 24)) +
    theme(legend.text = element_text(size = 24))+
    theme(legend.title = element_text(size = 24)) +
    theme(legend.position="top") 
ggsave("/home/jiawen/myMLnet/benchmark/Beeline/outputs/results_PPR.pdf",device = cairo_pdf,width = 2000,height = 2000,unit = "px")
