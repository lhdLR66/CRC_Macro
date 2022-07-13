setwd('E:/work/Apo_gene_choose/')

library(Seurat)
library(tidyverse)
library(reshape2)
library(psych)
library(ggpubr)
library(ggthemes)
library(clusterProfiler)
library(org.Hs.eg.db)
library(VennDiagram)
library(corrplot)
library(ggdendro)
library(survival)
library(survminer)
library(glmnet)
library(pheatmap)


#数据读取
##单细胞数据
macro_sub <- readRDS('macro_4_umap.rds')
macro_dmmr <- macro_sub[,macro_sub$MMRStatus=='MMRd']
macro_pmmr <- macro_sub[,macro_sub$MMRStatus=='MMRp']
apo_score <- read.csv('../Hallmark/HALLMARK_APOPTOSIS_score.csv',header = T,check.names = F,na.strings = '',row.names = 1)
apo_gene <- read.csv('Apoptosis_gene.csv',row.names = 1)

##bulk数据
dmmr_datsur <- read.csv('../data/bulk/dmmr_datsur.csv',header = T,row.names = 1,check.names = F)
pmmr_datsur <- read.csv('../data/bulk/pmmr_datsur.csv',header = T,row.names = 1,check.names = F)
crc_datsur <- read.csv('../data/bulk/crc_datsur.csv',header = T,row.names = 1,check.names = F)


#单细胞差异分析
markers_sc <- FindMarkers(macro_sub,ident.1 = 'MMRd',group.by = 'MMRStatus',
                          slot = 'data',
                          test.use = 'wilcox',
                          min.pct = 0.3,
                          random.seed = 1,
                          logfc.threshold = 0) 

##绘制火山图
DESeq2_result <- markers_sc
DESeq2_result$logP <- -log(DESeq2_result$p_val_adj,base = 20)
DESeq2_result$Group <- 'non-significant'
DESeq2_result$Group[which((DESeq2_result$p_val_adj<0.05) & (DESeq2_result$avg_log2FC>0.25))] <- 'highly expressed in dMMR'
DESeq2_result$Group[which((DESeq2_result$p_val_adj<0.05) & (DESeq2_result$avg_log2FC< (-0.25)))] <- 'highly expressed in pMMR'
table(DESeq2_result$Group)

DESeq2_result$Label <- ''
DESeq2_result <- DESeq2_result[order(DESeq2_result$avg_log2FC,decreasing = TRUE),]
up_genes <- head(rownames(DESeq2_result)[which(DESeq2_result$Group=='highly expressed in dMMR')],10)
down_genes <- tail(rownames(DESeq2_result)[which(DESeq2_result$Group=='highly expressed in pMMR')],10)
DEGs_top10 <- c(up_genes,down_genes)
DESeq2_result$Label[match(DEGs_top10,rownames(DESeq2_result))] <- DEGs_top10

plot1 <- ggscatter(DESeq2_result,
                   x='avg_log2FC',font.family = 'sans',
                   y='logP',
                   color = 'Group',
                   palette = c('#6768AB','#D5837B','#BBBBBB'),
                   size = 3,
                   label = 'Label',
                   repel = T,
                   font.label = 8)+theme_base()+
  geom_hline(yintercept = 1,linetype='dashed')+
  geom_vline(xintercept = c(-0.25,0.25),linetype='dashed')+
  theme(legend.position = 'top')+
  geom_point(size = 3, alpha = 0)+
  xlab('Log fold change')+ylab('-Log20(adj P-value)')
plot1
ggsave('sc_MMR_volcano.pdf',plot1,width = 6,height = 6)
DEG_MMR <-  filter(DESeq2_result,Group!='non-significant')
write.csv(DEG_MMR,'DEG_sig_0.05.csv')
DEGs <- rownames(DEG_MMR)#284个差异基因


#go富集分析
ENSEMBL <- bitr(DEGs,fromType="SYMBOL",toType="ENTREZID", OrgDb="org.Hs.eg.db",drop = T) 
go<- enrichGO(ENSEMBL$ENTREZID, 
                OrgDb = 'org.Hs.eg.db',
                ont = 'BP',readable = T,
                pvalueCutoff = 0.05,
                pAdjustMethod = 'BH', 
                minGSSize = 5,
                maxGSSize = 500,
                qvalueCutoff = 0.05)
go_fil <- go_fil %>% subset(.,FDR<0.05)
go_fil <- go_fil[order(go_fil$Count,decreasing = T),]
write.csv(go_fil,'kegg_result.csv')

##go结果可视化
go_fil2 <- go_fil[1:20,]
datbar <- go_fil2[,c(2,3,13)]
datbar$Term <- factor(datbar$Term,levels = rev(datbar$Term))
plot2 <- ggplot(datbar,aes(x=Count,y=Term,fill=FDR))+
         geom_bar(stat = 'identity',mapping = aes(fill=FDR))+
         theme_test()+
         scale_colour_gradient(low="#B5BFC9",high="#536887",aesthetics = "fill")+
  theme(axis.text = element_text(face='bold',size = 8,angle = 0,family = 'sans',colour = 'black'),
        legend.position = c(0.9,0.5))+
  xlab('Gene Counts')+ylab('')
plot2      
ggsave('go_enrich.pdf',plot1,width = 8,height = 6)

#lasso回归筛选基因
##表达矩阵与凋亡分数构建
dmmr_deg_expr_sc <- macro_dmmr@assays[["RNA"]]@data[DEGs,] %>% as.matrix()
dmmr_aposcore_sc <- data.frame(apo_score[colnames(dmmr_deg_expr_sc),],row.names = colnames(dmmr_deg_expr_sc))

pmmr_deg_expr_sc <- macro_pmmr@assays[["RNA"]]@data[DEGs,] %>% as.matrix()
pmmr_aposcore_sc <- data.frame(apo_score[colnames(pmmr_deg_expr_sc),],row.names = colnames(pmmr_deg_expr_sc))

##dmmr
set.seed(9)
{fit <- glmnet(x = t(dmmr_deg_expr_sc), y = as.matrix(dmmr_aposcore_sc),family = "gaussian",alpha = 1)
  plot(fit,xvar = "lambda")
  cv_fit <- cv.glmnet(x = t(dmmr_deg_expr_sc), y = as.matrix(dmmr_aposcore_sc),family="gaussian",type.measure="deviance",alpha = 1,nfolds = 3)
  plot(cv_fit)
  
  coefficient1 <- coef(cv_fit, s = cv_fit$lambda.min)
  Active.Index <- which(as.numeric(coefficient1) != 0)
  active.coefficient <- as.numeric(coefficient1)[Active.Index]
  sig_gene_apo_dmmr <- rownames(coefficient1)[Active.Index][-1]
  sig_gene_apo_dmmr}

##pmmr
set.seed(3)
{fit <- glmnet(x = t(pmmr_deg_expr_sc), y = as.matrix(pmmr_aposcore_sc),family = "gaussian",alpha = 1)
  plot(fit,xvar = "lambda")
  cv_fit <- cv.glmnet(x = t(pmmr_deg_expr_sc), y = as.matrix(pmmr_aposcore_sc),family="gaussian",type.measure="deviance",alpha = 1,nfolds = 3)
  plot(cv_fit)
  
  coefficient2 <- coef(cv_fit, s = cv_fit$lambda.min)
  Active.Index <- which(as.numeric(coefficient2) != 0)
  active.coefficient <- as.numeric(coefficient2)[Active.Index]
  sig_gene_apo_pmmr <- rownames(coefficient2)[Active.Index][-1]
  sig_gene_apo_pmmr}

inter_gene <- intersect(sig_gene_apo_dmmr,sig_gene_apo_pmmr)
inter_gene
write.csv(inter_gene,'Apo_related_DEG.csv')

##lasso系数热图可视化
coeff <- data.frame(dMMR=coefficient1[-1,],pMMR=coefficient2[-1,])
coeff <- coeff[union(sig_gene_apo_pmmr,sig_gene_apo_dmmr),] %>% round(.,5)
coeff_fil <- coeff
coeff_fil[coeff_fil<(-0.01)] <- (-0.01)
coeff_fil[coeff_fil>0.01] <- 0.01

coeff[coeff==0] <- NA
coeff_fil[coeff_fil==0] <- NA

coeff <- rbind(coeff[!is.na(rowSums(coeff)),],coeff[is.na(coeff$dMMR),],coeff[is.na(coeff$pMMR),])
coeff_fil <- coeff_fil[rownames(coeff),]

annotation_row <- data.frame(row.names = rownames(coeff_fil),Group=c(rep('Related to apoptosis in both dMMR and pmmr',25),
                                                              rep('Related to apoptosis in pMMR',24),
                                                              rep('Related to apoptosis in dMMR',42)))

plot3 <- pheatmap(coeff_fil,border_color = 'white',na_col = 'white',
                  annotation_names_row = FALSE,legend_labels = seq(-0.01,0.01,0.005),
                  show_rownames = T,show_colnames = T,
                  angle_col = 0,
                  display_numbers = coeff,
                  fontsize_row = 16,fontsize_col = 16,fontsize_number = 16,
                  cellwidth = 100,
                  cellheight = 19,
                  cluster_cols = F,cluster_rows = F,
                  color = colorRampPalette(c("#4E6CE4","white","#D69CE9"))(100),
                  use_raster = TRUE)

ggsave('lasso_coeff.pdf',plot3,width = 6,height = 25)




##交集基因Venn图展示
venn.diagram(list('dMMR'=sig_gene_apo_dmmr,'pMMR'=sig_gene_apo_pmmr),
             resolution = 800, imagetype = "png",
             filename = "Apogene_intersection.tiff",##韦恩图的名字
             lty = 1,
             lwd = 1,
             col = "black",  ##圈的颜色
             alpha = 0.60,
             fill=c("#6768AB","#D5837B"),##对应每个圈的颜色，有几个数据集，就需要有相应数量的颜色
             cat.col = "black",##此处设置每个数据集的名称颜色，也可以使用c（）函数输入三种颜色
             cat.cex = 0.8,
             cat.fontface = "bold",
             margin = 0.1,
             cex = 1,main.fontface = 'bold',main.cex = 0.8,main.fontfamily = 'sans',main.pos = c(0.5,1),
             main = "Intersection of apoptosis related gene in dMMR and pMMR")


#DEG与凋亡基因相关性(推断DEG可能的功能)
##单细胞DEG、apoptosis表达矩阵构建
crc_deg_expr <- macro_sub@assays[["RNA"]]@data[inter_gene,] %>% as.matrix()
crc_apo_expr <- macro_sub@assays[["RNA"]]@data[apo_gene$x,] %>% as.matrix()

##相关性分析
cor_interdeg_apo_crc <- corr.test(t(crc_deg_expr),t(crc_apo_expr),method = 'pearson',adjust = 'fdr') 
r_crc <- cor_interdeg_apo_crc$r 
p_crc <- cor_interdeg_apo_crc$p.adj 

##区分促凋亡和促生存的基因
table(rowSums(r_crc[,1:8]>0)>4&rowSums(r_crc[,9:16]>0)>4)
r_crc_fil <- r_crc[!(rowSums(r_crc[,1:8]>0)>4&rowSums(r_crc[,9:16]>0)>4),]#去除正相关较多的DEG，这类基因和两类apo基因均正相关较多，效应不明显
table(rowSums(r_crc_fil[,1:8]<0)>4&rowSums(r_crc_fil[,9:16]<0)>4)
r_crc_fil <- r_crc_fil[!(rowSums(r_crc_fil[,1:8]<0)>4&rowSums(r_crc_fil[,9:16]<0)>4),]#去除负相关较多的DEG，这类基因和两类apo基因均负相关较多，效应不明显

pos_relatedtoapo_deg <- rownames(r_crc_fil)[rowSums(r_crc_fil[,1:8]>0)>rowSums(r_crc_fil[,9:16]>0)]#与促凋亡基因正相关更多，认为这些DEG与促凋亡功能有关
neg_relatedtoapo_deg <- rownames(r_crc_fil)[rowSums(r_crc_fil[,1:8]>0)<rowSums(r_crc_fil[,9:16]>0)]#与促存活基因正相关更多，认为这些DEG与促存活功能有关

deg_final <- c(pos_relatedtoapo_deg,neg_relatedtoapo_deg)
write.csv(deg_final,'Apo_related_DEG_final.csv')

deg_final_anno <- DEG_MMR[deg_final,]
deg_final_anno <- mutate(deg_final_anno,gene=rownames(deg_final_anno),relevance=c(rep('positive related to apoptosis',2),rep('positive related to survival',19)))
deg_final_anno <- deg_final_anno[,c(2,5,7,9,10)]
deg_final_anno <- deg_final_anno[order(deg_final_anno$relevance,deg_final_anno$Group),]
write.csv(deg_final_anno,'Apo_related_DEG_final_anno.csv')

##相关性热图
r_crc_fil <- r_crc_fil[deg_final,]

apo_gene2 <- c('BAX','AIFM2','BID','HRK','CASP8','BBC3','FASLG','PIDD1',
               'CFLAR','NFKB1','BIRC2','BIRC3','XIAP','RELA','GADD45A','PTPN13')
r_crc_fil <- r_crc_fil[,apo_gene2]
p_crc_fil <- p_crc[rownames(r_crc_fil),colnames(r_crc_fil)]

pdf('corplot_deg_apo.pdf',width = 8,height = 8)
corrplot(r_crc_fil, p.mat = p_crc_fil,method = 'color',is.corr = FALSE,
         tl.cex = 0.5,tl.col = 'black',tl.offset = 0.5,tl.srt = 90,cl.length = 11,
         col.lim = c(-0.25,0.25),cl.cex = 1,cl.align.text = 'l',cl.ratio = 0.1,
         col = c(colorRampPalette(colors = c("#005B43","white","#7E155F"))(100)),
         sig.level = c(.001, .01, .05),outline="white",
         insig = "label_sig",pch.cex = 0.7, pch.col = "black")
dev.off()




#筛选核心基因
hub_gene <- subset(deg_final_anno,(deg_final_anno$relevance=='positive related to survival'&deg_final_anno$Group=='highly expressed in pMMR')|
                                  (deg_final_anno$relevance=='positive related to apoptosis'&deg_final_anno$Group=='highly expressed in dMMR'))
write.csv(hub_gene,'Hub_gene.csv')
