setwd('E:/work/Cluster/second/')

library(Seurat)
library(clustree)
library(tidyverse)
library(plyr)
library(ggstatsplot)
library(reshape2)
library(pheatmap)
library(ggplot2)
library(ggforce)
library(schex)


#亚群细分聚类
#macro_all <- readRDS('../../Process/macro_seu_fil.rds')
#macro_UMAP <- readRDS('../first/macro_UMAP.rds')
#macro_0_UMAP <- readRDS('macro_0_UMAP.rds')
#markers_paper <- read.csv('../first/markers_paper.csv',header = TRUE)
colorboard <- c('#9A3671','#95A238','#25BCCD','#E1C392')


#cluster0(M2)二次聚类
macro_0 <- macro_all[,macro_UMAP@meta.data[["seurat_clusters"]]==0]

##标准化、log转换
macro_0 <- NormalizeData(macro_0, normalization.method = "LogNormalize", scale.factor = 10000)

##鉴定高变基因
i0 <- 2000
macro_0 <- FindVariableFeatures(macro_0, selection.method = "vst", nfeatures = i0)
hvg0<- VariableFeatures(macro_0)
write.csv(hvg0,'hvg0_second.csv')

##scale
macro_0 <- ScaleData(macro_0,features = rownames(macro_0))

##主成分
macro_0 <- RunPCA(macro_0)

##选择最佳pca成分
macro_0 <- JackStraw(macro_0, num.replicate = 100)
macro_0 <- ScoreJackStraw(macro_0,dims = 1:20)

plot1 <- JackStrawPlot(macro_0, dims = 1:20)#选20个主成分
plot2 <- ElbowPlot(object = macro_0,ndims = 30,reduction = 'pca')

pdf(paste('JackStrawPlot+ElbowPlot_',i0,'.pdf',sep = ''),width = 12,height = 5)
plot1+plot2
dev.off()


#降维
macro_0 <- FindNeighbors(macro_0, dims = 1:12,reduction = 'pca')#调整PCs

##确定resolution 0.1-1.5
#res.used <- seq(0.1,1.5,by=0.1)
#sce <- FindClusters(object = macro_0, verbose = T, resolution = res.used)
#clus.tree.out <- clustree(sce,prefix = 'RNA_snn_res.') +
#  theme(legend.position = "bottom") + 
#  scale_color_brewer(palette = "Set1") +
#  scale_edge_color_continuous(low = "grey80", high = "red")

#pdf(paste('clustertree_0_',i0,'.pdf',sep = ''),width = 30,height = 30)
#clus.tree.out
#dev.off()

##聚类
macro_0 <- FindClusters(macro_0, resolution = 0.2)#调整cluster
macro_0_cluster <- macro_0
table(macro_0_cluster@meta.data[["seurat_clusters"]])

##确定最佳n.neighbors
#map(seq(5,50,5),function(x) { macro_0_cluster %>%  RunUMAP(n.neighbors = x,n.epochs=500,dims = 1:12) %>% DimPlot()}) %>% cowplot::plot_grid(plotlist = .)

macro_0_UMAP<- RunUMAP(macro_0, dims = 1:12,reduction = 'pca',seed.use = 1,n.neighbors = 40)

##细胞命名
cellnames <- mapvalues(macro_0_UMAP$seurat_clusters,from = c(0,1,2,3),to = c('C1QC+ TAM','M2c like TAM','M2b like TAM','M2a like TAM'))
macro_0_UMAP$seurat_clusters3 <- cellnames


##可视化
plot3 <- DimPlot(macro_0_UMAP, reduction = "umap",
                   label = F,cols = colorboard,
                   pt.size = 0.7,group.by = 'seurat_clusters3',
                   shuffle = TRUE,repel = TRUE)+
           ggtitle("M2 like TAM") +theme_void()+
           theme(text=element_text(size=12, family="sans",face = 'bold'),
                 plot.title = element_text(hjust = 0.5))
           

points <-data.frame(macro_0_UMAP@reductions$umap@cell.embeddings, cluster=macro_0_UMAP$seurat_clusters3)
theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
circle <- cbind(cos(theta), sin(theta))

aux <- function(x, one, two, prob=0.8) {
  if(nrow(x) <= 2) {
    return(NULL)
  }
  sigma <- var(cbind(x[,one], x[,two]))
  mu <- c(mean(x[,one]), mean(x[,two]))
  ed <- sqrt(qchisq(prob, df = 2))
  data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'))
}
ell <- ddply(points, "cluster", aux, one="UMAP_1", two="UMAP_2",prob=0.1)
plot3 <- plot3+geom_mark_ellipse(data=ell, aes(x=X1, y=X2, label=cluster,group=cluster),inherit.aes = F,color=NA) +NoLegend()
plot3

pdf(paste('cluster_0_',i0,'_0.2_UMAP','.pdf',sep = ''),width = 6,height = 6)
plot3
dev.off()

##差异分析FindAllMarkers
macro_0_markers <- FindAllMarkers(object = macro_0_UMAP, test.use = 'wilcox',slot = 'data',
                                  random.seed = 2,
                                  return.thresh = 0.05, 
                                  logfc.threshold = 0.25,only.pos = TRUE) %>% filter(.,p_val_adj<0.05)
macro_0_markers <- arrange(macro_0_markers,cluster,-avg_log2FC)
table(macro_0_markers$cluster)

macro_0_markers_fil <- data.frame(NULL)
for (i in unique(macro_0_markers$gene)) {
  sub1 <- macro_0_markers[macro_0_markers$gene==i,]
  sub2 <- macro_0_markers[macro_0_markers$avg_log2FC==max(sub1$avg_log2FC),]
  macro_0_markers_fil <- rbind(macro_0_markers_fil,sub2)
}
macro_0_markers_fil <- arrange(macro_0_markers_fil,cluster,-avg_log2FC)
table(macro_0_markers_fil$cluster)

datExpr_heatmap <- AverageExpression(macro_0_UMAP,group.by = 'seurat_clusters3',
                                     features = macro_0_markers_fil$gene,
                                     slot = 'data')$RNA %>% as.data.frame()
datExpr_heatmap$cluster <- macro_0_markers_fil$cluster[match(rownames(datExpr_heatmap),macro_0_markers_fil$gene)]
table(datExpr_heatmap$cluster)

datExpr_heatmap_sub <- data.frame(NULL)
for (i in as.numeric(levels(datExpr_heatmap$cluster))) {
  sub1 <- datExpr_heatmap[datExpr_heatmap$cluster==i,]
  sub2 <- apply(sub1[,1:4], 1, max)
  sub3 <- datExpr_heatmap[datExpr_heatmap$cluster==i,i+1]==sub2
  sub4 <- sub1[sub3,]
  datExpr_heatmap_sub <- rbind(datExpr_heatmap_sub,sub4)
}
datExpr_heatmap <- datExpr_heatmap_sub[,1:4]
datExpr_heatmap <- datExpr_heatmap %>% t() %>% scale() %>% t()

macro_0_markers_fil <- macro_0_markers_fil[match(rownames(datExpr_heatmap),macro_0_markers_fil$gene),]
table(macro_0_markers_fil$cluster)
write.csv(macro_0_markers_fil,'macro_markers_second_logfc0.25.csv')

##第二次分群平均表达量热图
plot4 <- pheatmap(datExpr_heatmap,
                    #annotation_col = colanno,
                    show_colnames = TRUE,show_rownames = FALSE,
                    #annotation_colors = coloranno,
                    cellheight = 0.3,cellwidth = 20,
                    fontsize = 6,
                    cluster_rows = FALSE,
                    cluster_cols = FALSE,
                    color = colorRampPalette(c("#3A54A5", "white", "#ED2225"))(100),
                    use_raster = TRUE)
pdf('Heatmap.pdf',height = 10)
plot4
dev.off()

##基因表达featureplot
macro_0_UMAP <- make_hexbin(macro_0_UMAP,50, dimension_reduction = "umap")
##M2a/M2b/M2c
plot5 <- (plot_hexbin_feature(macro_0_UMAP,type = 'scale.data',feature = 'IGF1',action = 'mean',title = 'IGF1',xlab = '',ylab = '')/
          plot_hexbin_feature(macro_0_UMAP,type = 'scale.data',feature = 'KDM6B',action = 'mean',title = 'KDM6B',xlab = '',ylab = ''))|
          (plot_hexbin_feature(macro_0_UMAP,type = 'scale.data',feature = 'CXCL2',action = 'mean',title = 'CXCL2',xlab = '',ylab = '')/
          plot_hexbin_feature(macro_0_UMAP,type = 'scale.data',feature = 'CXCL3',action = 'mean',title = 'CXCL3',xlab = '',ylab = ''))|
          (plot_hexbin_feature(macro_0_UMAP,type = 'scale.data',feature = 'MRC1',action = 'mean',title = 'MRC1',xlab = '',ylab = '')/
          plot_hexbin_feature(macro_0_UMAP,type = 'scale.data',feature = 'TLR8',action = 'mean',title = 'TLR8',xlab = '',ylab = ''))
ggsave('Featureplot_M2_all.pdf',plot5,width = 15,height = 10)

##基因表达vlnplot
##M2a/M2b/M2c
plot6<- VlnPlot(macro_0_UMAP,
                   features = c('IGF1','KDM6B','CXCL2','CXCL3','MRC1','TLR8','C1QC','C1QA'),
                   stack = TRUE,flip = T,
                   group.by = 'seurat_clusters3',
                   split.by = 'seurat_clusters3',
                   cols = colorboard,
                   same.y.lims = FALSE)+
           theme(axis.text.x = element_text(hjust = 1,size = 8,family = "sans",face = 'bold'),
                  axis.text.y = element_text(hjust = 1,vjust = 0,size = 9,family = "sans",face = 'bold'),
                 legend.position = 'none')+
           xlab('')+ylab('')
plot6
ggsave('Vlnplot_M2_all.pdf',plot6,width = 6,height = 8)


#提取cluster编号
##导出编号
clustername_0 <- data.frame(cell=names(macro_0_UMAP$seurat_clusters),cluster=macro_0_UMAP$seurat_clusters)#4群

##替换字符
clustername_0$cluster <- mapvalues(x = clustername_0$cluster,from = 0:3,to = LETTERS[1:4])

seurat_clusters <- data.frame(row.names = colnames(macro_UMAP),cell=colnames(macro_UMAP),cluster=as.character(macro_UMAP@meta.data[["seurat_clusters"]]))
table(seurat_clusters$cluster)
seurat_clusters$cluster[seurat_clusters$cluster==0] <- as.character(clustername_0$cluster)
table(seurat_clusters$cluster)

##替换命名
seurat_clusters$cluster <- mapvalues(seurat_clusters$cluster,names(sort(table(seurat_clusters$cluster),decreasing = TRUE)),1:10)
seurat_clusters$cluster<- as.numeric(seurat_clusters$cluster)
seurat_clusters$cluster<- as.factor(seurat_clusters$cluster)
table(seurat_clusters$cluster)

seurat_clusters$name <- mapvalues(seurat_clusters$cluster,from = 1:10,to = c('SPP1+ TAM','C1QC+ TAM','FTL+ Macro','M2c like TAM','M1 like Macro',
                                                                             'MKI67+ Macro','GNLY+ Monolike Macro','M2b like TAM','M2a like TAM','STMN1+ Macro'))
seurat_clusters$name <- factor(seurat_clusters$name,levels = c('SPP1+ TAM','C1QC+ TAM','M2a like TAM','M2b like TAM','M2c like TAM',
                                                               'M1 like Macro','GNLY+ Monolike Macro',
                                                               'MKI67+ Macro','STMN1+ Macro',
                                                               'FTL+ Macro'))

##添加新分群至seurat对象
macro_UMAP <- AddMetaData(macro_UMAP,metadata = seurat_clusters$cluster,col.name = 'seurat_clusters2')
macro_UMAP <- AddMetaData(macro_UMAP,metadata = seurat_clusters$name,col.name = 'seurat_clusters4')
saveRDS(macro_UMAP,'macro_UMAP.rds')
