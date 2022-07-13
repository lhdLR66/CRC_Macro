setwd('E:/work/Hub_gene_futher/ICB_response/')

library(easier)
library(maftools)
library(ggpubr)
library(ggstatsplot)
library(reshape2)
library(tidyverse)
library(psych)
library(plyr)
library(GSVA)

colorboard <- c('#6768AB','#D5837B','#F9959D','#FFBEF2')

#数据读取
crc_counts <- read.csv('E:/work/data/bulk/crc_datexpr_counts.csv',row.names = 1,header = T,check.names = F)
crc_tpm <- read.csv('E:/work/data/bulk/crc_datexpr_tpm.csv',row.names = 1,header = T,check.names = F)
crc_status <- read.csv('E:/work/data/bulk/crc_datstatus.csv',row.names = 2,header = T,check.names = F)
crc_logfpkm <- read.csv('E:/work/data/bulk/crc_datexpr.csv',header = T,row.names = 1,check.names = F)
hub_gene <- read.csv('../Choose_gene/Hub_gene.csv',row.names = 1)

coad_maf <- read.maf('E:/work/data/bulk/rawdata/TCGA.COAD.mutect.03652df4-6090-4f5a-a2ff-ee28a37f9301.DR-10.0.somatic.maf.gz',isTCGA = T)
read_maf <- read.maf('E:/work/data/bulk/rawdata/TCGA.READ.mutect.faa5f62a-2731-4867-a264-0e85b7074e87.DR-10.0.somatic.maf.gz',isTCGA = T)

markers_4<- read.csv('E:/work/Cluster/total/macro_markers_all_logfc0.25.csv',row.names = 1,header = T,check.names = F) %>% filter(.,cluster=='M2c like TAM')
markers_4 <- list(as.vector(markers_4[1:30,]$gene))


#计算TMB
coad_tmb <- tmb(coad_maf)
read_tmb <- tmb(read_maf)
crc_tmb <- rbind(coad_tmb,read_tmb)

TMB <- crc_tmb$total_perMB
names(TMB) <- crc_tmb$Tumor_Sample_Barcode


#数据整合与构建
TMB <- TMB[intersect(names(TMB),colnames(crc_counts))]
crc_counts <- crc_counts[,names(TMB)]
crc_tpm <- crc_tpm[,names(TMB)]
crc_status_sub <- crc_status[names(TMB),]
crc_logfpkm_sub <- crc_logfpkm[,names(TMB)]

RNA_counts <- crc_counts
RNA_tpm <- crc_tpm


#免疫治疗预测
##计算免疫反应相关分数
hallmarks_of_immune_response <- c("CYT", "Roh_IS", "chemokines", "Davoli_IS", "IFNy", "Ayers_expIS", "Tcell_inflamed", "RIR", "TLS")
immune_response_scores <- compute_scores_immune_response(RNA_tpm = RNA_tpm, 
                                                         selected_scores = hallmarks_of_immune_response)#tpm数据

##计算TME相关定量描述（共五个）
cell_fractions <- compute_cell_fractions(RNA_tpm = RNA_tpm)
pathway_activities <- compute_pathway_activity(RNA_counts = RNA_counts,
                                               remove_sig_genes_immune_response = TRUE)
tf_activities <- compute_TF_activity(RNA_tpm = RNA_tpm)
lrpair_weights <- compute_LR_pairs(RNA_tpm = RNA_tpm,
                                   cancer_type = "CRC")

ccpair_scores <- compute_CC_pairs(lrpairs = lrpair_weights, 
                                  cancer_type = "CRC")

##免疫反应预测
predictions <- predict_immune_response(pathways = pathway_activities,
                                       immunecells = cell_fractions,
                                       tfs = tf_activities,
                                       lrpairs = lrpair_weights,
                                       ccpairs = ccpair_scores,
                                       cancer_type = 'CRC', 
                                       verbose = TRUE)

output_eval_no_resp <- assess_immune_response(predictions_immune_response = predictions,
                                              RNA_tpm = RNA_tpm,
                                              TMB_values = TMB,
                                              easier_with_TMB = "weighted_average",
                                              weight_penalty = 0.5)

##提取easier评分
easier_derived_scores <- retrieve_easier_score(predictions_immune_response = predictions,
                                               TMB_values = TMB,
                                               easier_with_TMB = c("weighted_average", 
                                                                   "penalized_score"),
                                               weight_penalty = 0.5)

##生物标志物体现
output_biomarkers <- explore_biomarkers(pathways = pathway_activities,
                                        immunecells = cell_fractions,
                                        lrpairs = lrpair_weights,
                                        tfs = tf_activities,
                                        ccpairs = ccpair_scores,
                                        cancer_type = 'CRC')


#d/p组免疫治疗likelihood及S100A6表达箱线图
##将pmmr分成三组(low/median/high)
{datbox <- data.frame(row.names = names(TMB),likelihood=easier_derived_scores$w_avg_score,expression=t(crc_logfpkm_sub['S100A6',])[,1],group=crc_status_sub$dMMR.or.pMMR)
datbox_sub <- datbox[datbox$group=='pMMR',]
group_kmeans <- kmeans(datbox_sub$likelihood,3)
datbox_sub$group <- group_kmeans[["cluster"]]
datbox_sub$group <- mapvalues(datbox_sub$group,from = c(3,1,2),to = c('pMMR high','pMMR median','pMMR low'))
datbox <- left_join(datbox,datbox_sub,'likelihood')
datbox$group.y <- mapvalues(datbox$group.y,NA,'dMMR')
datbox <- datbox[,-c(3,4)]
colnames(datbox) <- c('score','expression','group')
datbox$group <- factor(datbox$group,levels = c('dMMR','pMMR high','pMMR median','pMMR low'))}

##likelihood箱线图
plot1 <- ggplot(aes(x = group, y = score,fill = group),data = datbox) +
          stat_boxplot(geom = 'errorbar',width=0.3,cex=1)+
          geom_boxplot(outlier.size = -1,width=0.5,
                       fill=colorboard,alpha=c(1,1,0.8,0.6)) +
          scale_y_continuous(name = "likelihood score")+
          scale_x_discrete(name = "") +
          ggtitle("Likehood of immune therapy") +
          theme_bw()+
          theme(plot.title = element_text(size = 10, face =  "bold",family = 'sans',hjust = 0.5),
                text = element_text(size = 12),
                axis.title = element_text(face="bold"),
                axis.text.x=element_text(size = 10,angle = 45,family = 'sans',vjust = 1,hjust = 1,face = 'bold'),
                panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ 
          geom_signif(comparisons = list(c(1,2)),test="wilcox.test",
                      tip_length = 0.02,size = 0.8,textsize = 4,y_position = 0.45)+
  geom_signif(comparisons = list(c(2,3)),test="wilcox.test",
              tip_length = 0.02,size = 0.8,textsize = 4,y_position = 0.4)+
  geom_signif(comparisons = list(c(4,3)),test="wilcox.test",
              tip_length = 0.02,size = 0.8,textsize = 4,y_position = 0.25)
plot1
ggsave('likelihoodofICB.pdf',plot1,height = 6,width = 3)

##S100A6表达箱线图
datbox2 <- data.frame(group=crc_status$dMMR.or.pMMR,expression=t(crc_logfpkm['S100A6',])[,1])
plot2 <- ggplot(aes(x = group, y = expression,fill = group),data = datbox2) +
  stat_boxplot(geom = 'errorbar',width=0.3,cex=1)+
  geom_boxplot(outlier.size = -1,width=0.5,fill=c('#6768AB','#D5837B')) +
  scale_y_continuous(name = "Expression",)+
  scale_x_discrete(name = "") +
  ggtitle("Expression of S100A6") +
  theme_bw()+
  theme(plot.title = element_text(size = 10, face =  "bold",family = 'sans',hjust = 0.5),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 10,angle = 45,family = 'sans',vjust = 1,hjust = 1,face = 'bold'),
        panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ 
  geom_signif(comparisons = list(c(1,2)),test="wilcox.test",
              tip_length = 0.02,size = 0.8,textsize = 4,y_position = 13.5)
plot2
ggsave('expressionofS100A6_bulk.pdf',plot2,width = 3,height = 6)


#cluster4 免疫浸润与免疫治疗关系
##ssgsea免疫浸润
ssgsea <- gsva(as.matrix(crc_logfpkm_sub), markers_4,
               verbose=TRUE,
               method="ssgsea",
               kcdf="Gaussian",
               min.sz=1,
               max.sz=1000)

datbox3 <- datbox
datbox3$expression <- t(ssgsea)
colnames(datbox3)[2] <- 'infscore'

plot3 <- ggscatter(datbox3,x = "score", y = 'infscore',size = 2,
                   xlab = 'Immunotherapy response score',ylab = 'Infiltration of M2c like TAMs',color = '#E69947',
                   add = "reg.line", conf.int = TRUE,    
                   add.params = list(fill = "darkgray"))+
         stat_cor(method = "pearson")+theme_test()+
         facet_wrap(~group)
plot3
ggsave('cor_therapy_inf.pdf',plot3,width = 6,height = 6)

##密度图
plot4 <- ggplot(datbox,aes(x=score,fill=group))+
         geom_density(alpha=0.8,aes(x=score))+
         scale_fill_manual(values = colorboard)+
         theme_bw()+
         xlab('Response score')+ylab('Density')
plot4
ggsave('density.pdf',plot4,width = 8,height = 5)
