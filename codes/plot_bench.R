library(ggplot2)
library(dplyr)
library(forcats)
results <- readxl::read_xlsx("./results.xlsx")

aurocs <- results[,c(1,seq(2,20,2))]
auprcs <- results[,c(1,seq(3,21,2))]
auroc_long <- tidyr::gather(aurocs, key = "sample",value = "AUROC",-'Model')
auprc_long <- tidyr::gather(auprcs, key = "sample",value = "AUPRC",-'Model')


auroc_long %>% 
  ggplot( aes(x=Model, y=AUROC, fill=Model)) + 
  geom_boxplot() +
  xlab("Method") +
  ylab("AUROC") +
  scale_y_continuous(limits = c(0.4,0.9))+
  theme_bw() +
  theme(text = element_text(family = "A")) +
  theme(axis.title = element_text(size = 24)) +
  theme(axis.text = element_text(size = 24)) +
  theme(legend.position="none") 

auprc_long %>% 
  ggplot( aes(x=Model, y=AUPRC, fill=Model)) + 
  geom_boxplot() +
  xlab("Method") +
  ylab("AUPRC") +
  scale_y_continuous(limits = c(0.4,0.9))+
  theme_bw() +
  theme(text = element_text(family = "A")) +
  theme(axis.title = element_text(size = 24)) +
  theme(axis.text = element_text(size = 24)) +
  theme(legend.position="none") 

