setwd('E:/work/Function/')
library(GSVA)
library(pheatmap)
library(tidyverse)
library(ggsci)
library(reshape2)
library(ggpubr)
library(Seurat)
library(fmsb)


#细胞注释(以功能确定细胞类型和)
##GSVA分析
macro_UMAP <- readRDS('../Cluster/second/macro_UMAP.rds')
datExpr <- macro_UMAP@assays[["RNA"]]@data %>% as.matrix()
#markers_paper <- read.csv('../Cluster/first/markers_paper.csv',header = T,check.names = F)

##gmt文件构建
#geneset <- read.csv('../data/scRNA/178341/rawdata/macro_geneset_v3.csv',header = FALSE,row.names = 1,na.strings = '') %>% t()
#geneset <- geneset[-1,]
#GSVA_gs <- melt(geneset,measure.vars = colnames(geneset),variable.name = 'Description',value.name = 'gene') %>% na.omit()
#GSVA_gs <- split(GSVA_gs$gene, GSVA_gs$Var2)

## all cells GSVA
#gsva.es <- gsva(datExpr, GSVA_gs,
#                verbose=TRUE,
#                method="gsva",
#                kcdf="Gaussian",
#                min.sz=1,
#                max.sz=1000)
#write.csv(gsva.es,paste('gsva.es_',length(GSVA_gs),'.csv',sep = ''))
gsva.es <- read.csv('gsva.es_39.csv',header = T,row.names = 1,check.names = F)
colnames(gsva.es) <- macro_UMAP@meta.data[["seurat_clusters4"]]
table(colnames(gsva.es))

gsva.es <- gsva.es[-c(7,10,12,13,31,36,39),]
rownames(gsva.es)[c(10,11)] <- c('M1 Macrophage Polarization','M2 Macrophage Polarization')

##不分组取平均
gsva.es_avg <- data.frame(row.names = rownames(gsva.es))
for (i in levels(macro_UMAP$seurat_clusters4)) {
  gsva.es_sub <- gsva.es[,colnames(gsva.es)==i]
  gsva.es_sub <- rowMeans(gsva.es_sub)
  gsva.es_avg <- cbind(gsva.es_avg,gsva.es_sub)
}
colnames(gsva.es_avg) <- levels(macro_UMAP$seurat_clusters4)

gsva.es_avg <- gsva.es_avg %>% t() %>% scale() %>% t()
summary(as.numeric(gsva.es_avg))

gsva_hclust <- hclust(dist(gsva.es_avg))
function_group <- data.frame(cutree(gsva_hclust,k=5))
gsva.es_avg <- gsva.es_avg[order(function_group$cutree.gsva_hclust..k...5.),]

gsva_hclust2 <- hclust(dist(t(gsva.es_avg)))
cluster_group <- data.frame(cutree(gsva_hclust2,k=4))
gsva.es_avg <- gsva.es_avg[,order(cluster_group$cutree.gsva_hclust2..k...4.)]

##d/p组富集分数wilcox差异分析
p <- data.frame(row.names = rownames(gsva.es_avg))
for (i in rownames(p)) {
  p_sub2 <- c()
  for (j in colnames(gsva.es_avg)[c(1:6,8:10)]) {
    gsva.es_sub <- gsva.es[i,colnames(gsva.es)==j]
    group <- macro_UMAP$MMRStatus[macro_UMAP$seurat_clusters4==j]
    gsva.es_sub <- gsva.es_sub %>% t() %>% as.data.frame() %>%  mutate(.,group=group)
    colnames(gsva.es_sub)[1] <- 'score'
    p_sub1 <- wilcox.test(score~group,gsva.es_sub,paired=F,alternative='two.sided')$p.value
    p_sub2 <- c(p_sub2,p_sub1)
  }
  p <- rbind(p,p_sub2)
}
p_adj <- p.adjust(melt(p)$value,method = "fdr") %>% matrix(nrow = 32,ncol = 9)
rownames(p_adj) <- rownames(gsva.es_avg)
colnames(p_adj) <- colnames(gsva.es_avg)[c(1:6,8:10)]
write.csv(p_adj,'cluster_group_function.csv')
p_adj <- ifelse(p_adj<0.05,'*','')
p_adj <- cbind(p_adj[,1:6],rep('',32),p_adj[,7:9])
colnames(p_adj)[7] <- 'M2a like TAM'

##分d/p组取平均
gsva.es_group <- data.frame(row.names = rownames(gsva.es_avg))
for (i in colnames(gsva.es_avg)) {
  sub1 <- gsva.es[,macro_UMAP$MMRStatus=='MMRd'&colnames(gsva.es)==i]
  sub2 <- gsva.es[,macro_UMAP$MMRStatus=='MMRp'&colnames(gsva.es)==i]
  sub3 <- ifelse(rowMeans(sub1)>rowMeans(sub2),'dMMR','pMMR')
  gsva.es_group <- cbind(gsva.es_group,sub3)
}
colnames(gsva.es_group) <- colnames(gsva.es_avg)
gsva.es_group[p_adj==''] <- ''
write.csv(gsva.es_group,'cluster_group_function2.csv')

##热图绘制
#annotation_color <- list(cluster=c('1'='#800020','2'='#9A3671','3'='#008C8C','4'='#95A238','5'='#FADA5E',
#                                   '6'='#E85827','7'='#B05923','8'='#25BCCD','9'='#E1C392','10'='#FF0000'))
number_color <- ifelse(gsva.es_group=='dMMR','black',ifelse(gsva.es_group=='pMMR','orange',NA))

plot1 <- pheatmap(gsva.es_avg,border_color = 'white',
                  annotation_names_row = FALSE,
                  show_rownames = T,show_colnames = T,
                  gaps_row = c(18,23,28,31),gaps_col = c(3,6,8),
                  angle_col = 45,
                  annotation_colors = annotation_color,
                  display_numbers = p_adj,
                  fontsize_row = 5,fontsize_col = 5,fontsize_number = 15,number_color = number_color,
                  cellwidth = 8,
                  cellheight = 8,
                  cluster_cols = F,cluster_rows = F,
                  color = colorRampPalette(c("#005B43","white","#7E155F"))(100),
                  use_raster = TRUE)
cairo_pdf('cluster_gsva.pdf',width = 5,height = 5)
plot1
dev.off()


#M2亚群功能分析
##gmt文件构建
#geneset_m2 <- read.csv('../data/scRNA/178341/rawdata/macro_geneset_v4.csv',header = FALSE,row.names = 1,na.strings = '') %>% t()
#geneset_m2 <- geneset_m2[-1,]
#GSVA_gs_m2 <- melt(geneset_m2,measure.vars = colnames(geneset_m2),variable.name = 'Description',value.name = 'gene') %>% na.omit()
#GSVA_gs_m2 <- split(GSVA_gs_m2$gene, GSVA_gs_m2$Var2)

## M2 GSVA
#datExpr_m2 <- datExpr[,macro_UMAP$seurat_clusters2 %in% c(2,4,8,9)]
#gsva.es_m2 <- gsva(datExpr_m2, GSVA_gs_m2,
#                verbose=TRUE,
#                method="gsva",
#                kcdf="Gaussian",
#                min.sz=1,
#                max.sz=1000)
#write.csv(gsva.es_m2,paste('gsva.es_',length(GSVA_gs_m2),'.csv',sep = ''))
gsva.es_m2 <- read.csv('gsva.es_8.csv',header = T,row.names = 1,check.names = F)
colnames(gsva.es_m2) <- macro_UMAP$seurat_clusters4[macro_UMAP$seurat_clusters2 %in% c(2,4,8,9)]
table(colnames(gsva.es_m2))

##取平均
gsva.es_avg_m2 <- data.frame(row.names = rownames(gsva.es_m2))
for (i in unique(colnames(gsva.es_m2))) {
  gsva.es_sub <- gsva.es_m2[,colnames(gsva.es_m2)==i]
  gsva.es_sub <- rowMeans(gsva.es_sub)
  gsva.es_avg_m2 <- cbind(gsva.es_avg_m2,gsva.es_sub)
}
colnames(gsva.es_avg_m2) <- unique(colnames(gsva.es_m2))

gsva.es_avg_m2 <- gsva.es_avg_m2 %>% t() %>% scale() %>% t()
summary(as.numeric(gsva.es_avg_m2))

##绘制雷达图
datradar <- t(gsva.es_avg_m2) %>% as.data.frame()
datradar <- rbind(rep(1.5,8),rep(-1.5,8),datradar)
rownames(datradar) <- c('max','min',unique(colnames(gsva.es_m2)))

pdf('M2_function_radar.pdf',width = 6,height = 5)
radarchart(datradar,
           axistype= 1,
           seg = 2,
           pty = 16,
           pcol = c('#9A3671','#25BCCD','#95A238','#E1C392'), pfcol = NULL,
           plwd = 3,plty = 1,cglty = 2,cglwd = 1,cglcol = 'black',axislabcol = 'black',
           vlcex = 0.5,caxislabels = c('-1.5','0','1.5'),calcex = 0.8,title = '')
dev.off()

##绘制分组柱状图
datBar <- melt(gsva.es_avg_m2)
colnames(datBar) <- c('pathway','cell','score')

plot2 <- ggplot(data=datBar,mapping=aes(x=pathway,y=score,fill=cell,width=0.6))+
          geom_bar(stat = 'identity',position = position_dodge(0.6))+
          scale_fill_manual(values = c('#9A3671','#25BCCD','#95A238','#E1C392'))+
          geom_hline(yintercept = 0,linetype="dashed")+
          theme_test()+
          ggtitle('')+
          xlab('')+ylab('GSVA z-score')+
          theme(axis.text = element_text(face='bold',size = 5,angle = 0,family = 'sans',colour = 'black'),
                axis.text.x =element_text(vjust = 1,angle = 45,hjust = 1,size = 5),
                plot.title = element_text(hjust = 0.5,size = 10),
                legend.position = 'right')
        
plot2
ggsave('M2_function_bar.pdf',plot2,width = 8,height = 4)

##M2巨噬细胞功能(完整版)
gsva.es_avg_m2_total <- gsva.es_avg[c(1,2,10,11,15,
                                      7,8,9,19,26,27,28,
                                      4,5,
                                      20,30,31),c(4,2,8,7)]
number_color <- ifelse(gsva.es_group=='dMMR','black',ifelse(gsva.es_group=='pMMR','orange',NA))
p_adj_m2 <- p_adj[c(1,2,10,11,15,
                    7,8,9,19,26,27,28,
                    4,5,
                    20,30,31),c(4,2,8,7)]
number_color_m2 <- number_color[c(1,2,10,11,15,
                                  7,8,9,19,26,27,28,
                                  4,5,
                                  20,30,31),c(4,2,8,7)]

plot3 <- pheatmap(gsva.es_avg_m2_total,border_color = 'white',
                  annotation_names_row = FALSE,
                  show_rownames = T,show_colnames = T,
                  gaps_row = c(5,12,14),
                  angle_col = 45,
                  display_numbers = p_adj_m2,
                  fontsize_row = 5,fontsize_col = 5,fontsize_number = 15,number_color = number_color_m2,
                  cellwidth = 20,
                  cellheight = 20,
                  cluster_cols = F,cluster_rows = F,
                  color = colorRampPalette(c("#005B43","white","#7E155F"))(100),
                  use_raster = TRUE)
cairo_pdf('cluster_gsva_m2.pdf',width = 7,height = 7)
plot3
dev.off()
