setwd('c:/Users/Lenovo/Desktop/lr')

library(tidyverse)
library(corrplot)
library(psych)
library(plyr)
library(reshape2)
library(tidyr)
library(patchwork)


#macro_UMAP <- readRDS('../Cluster/second/macro_UMAP.rds')
#macro_m2 <- macro_UMAP[,macro_UMAP$seurat_clusters2%in%c(2,4,8,9)]
macro_m2$seurat_clusters4 <- factor(macro_m2$seurat_clusters4,levels = c('C1QC+ TAM','M2b like TAM','M2c like TAM','M2a like TAM'))
macro_m2$TumorStage <- mapvalues(macro_m2$TumorStage,from = c( 'pT1','pT2','pT3','pt4a','pT4a','pT4b'),to = c('T1','T2','T3','T4','T4','T4'))
macro_m2$NodeStatus_detailed <- mapvalues(macro_m2$NodeStatus_detailed,
                                          from = c('N0','N0 (i+)* (isolated tumor cells)','N1','N1a','N1b','N1c','N2a','N2b'),
                                          to = c('N0','N0','N1','N1','N1','N1','N2','N2'))

m2_dmmr <- macro_m2[,macro_m2$MMRStatus=='MMRd']
m2_pmmr <- macro_m2[,macro_m2$MMRStatus=='MMRp']
colorboard <- c('#9A3671','#25BCCD','#95A238','#E1C392')


#免疫浸润差异
datBar <- matrix(as.numeric(table(macro_UMAP$MMRStatus,macro_UMAP$seurat_clusters4)),nrow = 2,dimnames = list(c('MMRd','MMRp'),1:10)) %>% as.data.frame() %>% t() 

##卡方检验
p_chi <- c()
for (i in c(1:2,4:10)) {
  datchi <- matrix(c(datBar[i,],sum(datBar[i,])*(6270 /13704),sum(datBar[i,])*(7434 /13704)),nrow = 2,ncol = 2)
  chi_result <- chisq.test(datchi)
  sub1 <- chi_result$p.value
  p_chi <- c(p_chi,sub1)
}

##Fisher精确检验
datfisher <- matrix(c(datBar[3,],sum(datBar[3,])*(6270 /13704),sum(datBar[3,])*(7434 /13704)),nrow = 2,ncol = 2)
fisher_result <- fisher.test(datfisher)
p_fisher <- fisher_result$p.value

p_all <- c(p_chi[1:2],p_fisher,p_chi[3:9]) %>% p.adjust(.,method = 'fdr')
log_p_all <- log(p_all,base = 20)
log_p_all <- ifelse(log_p_all<(-10),-10,log_p_all)
log_p_all[c(1,7,10)] <- -log_p_all[c(1,7,10)]

##绘制棒棒糖图
datBar <- data.frame(cluster=levels(macro_UMAP$seurat_clusters4),logp=log_p_all)
plot1 <- ggplot(datBar,aes(x=cluster,y=logp))+
  geom_hline(yintercept = c(-1,0,1),color=c('grey','black','grey'),linetype=c('dashed'))+
  geom_segment(aes(x=cluster,xend=cluster,y=0,yend=logp),size=1,color='darkgrey',linetype="solid")+
  geom_point(size=8,color=ifelse(log_p_all>0,'#6768AB','#D5837B'),
             fill=ifelse(log_p_all>1,'#6768AB',ifelse(log_p_all<(-1),'#D5837B','white')),
             shape=21)+
  ylab('log20p')+xlab('')+
  theme_test()
plot1
ggsave('cell_ratio_inf.pdf',plot1,width = 6,height = 4)


#细胞比例随肿瘤进展变化比例图
##M2巨噬细胞
##dMMR tumorstage
{cell_pro1 <- as.data.frame(table(m2_dmmr@meta.data[["seurat_clusters4"]],
                                  m2_dmmr$TumorStage)) %>% spread(.,key = Var1,value = Freq)

cell_pro1[,-1] <- (cell_pro1[,-1]/rowSums(cell_pro1[,-1])) 
cell_pro1 <- melt(cell_pro1,id.vars = 'Var2')
colnames(cell_pro1) <- c('group','cell','ratio')}

plot2 <- ggplot(cell_pro1 , mapping = aes(x = group, y = ratio,group=cell,color=cell)) + geom_line(cex=1)+
  geom_point(pch=20,cex=5)+
  theme_classic()+
  theme(axis.text = element_text(face='bold',size = 10,angle = 0,family = 'sans',colour = 'black'),
        axis.text.x =element_text(vjust = 0.5),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = '')+
  scale_color_manual(values = colorboard)+
  xlab('')+ylab('dMMR')+
  ggtitle('Changes in cell type composition of M2 like TAMs 
    at each tumor stage')
plot2 
 
##pMMR tumorstage
{cell_pro2 <- as.data.frame(table(m2_pmmr@meta.data[["seurat_clusters4"]],
                                  m2_pmmr$TumorStage)) %>% spread(.,key = Var1,value = Freq)
  
  cell_pro2[,-1] <- (cell_pro2[,-1]/rowSums(cell_pro2[,-1])) 
  cell_pro2 <- melt(cell_pro2,id.vars = 'Var2')
  colnames(cell_pro2) <- c('group','cell','ratio')}

plot3 <- ggplot(cell_pro2 , mapping = aes(x = group, y = ratio,group=cell,color=cell)) + geom_line(cex=1)+
  geom_point(pch=20,cex=5)+
  theme_classic()+
  theme(axis.text = element_text(face='bold',size = 10,angle = 0,family = 'sans',colour = 'black'),
        axis.text.x =element_text(vjust = 0.5),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = '')+
  scale_color_manual(values = colorboard)+
  xlab('')+ylab('pMMR')
  
plot3

##dMMR nodestage
{cell_pro3 <- as.data.frame(table(m2_dmmr@meta.data[["seurat_clusters4"]],
                                  m2_dmmr$NodeStatus_detailed)) %>% spread(.,key = Var1,value = Freq)
  cell_pro3[,-1] <- (cell_pro3[,-1]/rowSums(cell_pro3[,-1])) 
  cell_pro3 <- melt(cell_pro3,id.vars = 'Var2')
  colnames(cell_pro3) <- c('group','cell','ratio')}

plot4 <- ggplot(cell_pro3 , mapping = aes(x = group, y = ratio,group=cell,color=cell)) + geom_line(cex=1)+
  geom_point(pch=20,cex=5)+
  theme_classic()+
  theme(axis.text = element_text(face='bold',size = 10,angle = 0,family = 'sans',colour = 'black'),
        axis.text.x =element_text(vjust = 0.5),
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_color_manual(values = colorboard)+
  xlab('')+ylab('')+
  ggtitle('Changes in cell type composition of M2 like TAMs 
    at each lymph nodes stage')
plot4

##pMMR nodestage
{cell_pro4 <- as.data.frame(table(m2_pmmr@meta.data[["seurat_clusters4"]],
                                  m2_pmmr$NodeStatus_detailed)) %>% spread(.,key = Var1,value = Freq)
  cell_pro4[,-1] <- (cell_pro4[,-1]/rowSums(cell_pro4[,-1])) 
  cell_pro4 <- melt(cell_pro4,id.vars = 'Var2')
  colnames(cell_pro4) <- c('group','cell','ratio')}

plot5 <- ggplot(cell_pro4 , mapping = aes(x = group, y = ratio,group=cell,color=cell)) + geom_line(cex=1)+
  geom_point(pch=20,cex=5)+
  theme_classic()+
  theme(axis.text = element_text(face='bold',size = 10,angle = 0,family = 'sans',colour = 'black'),
        axis.text.x =element_text(vjust = 0.5),
        plot.title = element_text(hjust = 0.5,size = 15))+
  scale_color_manual(values = colorboard)+
  xlab('')+ylab('')
plot5
ggsave('cell_ratio_all.pdf',(plot2|plot4)/(plot3|plot5),height = 5,width = 12)
