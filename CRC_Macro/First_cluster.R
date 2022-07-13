setwd('E:/work/Cluster/first/')
library(Seurat)
library(clustree)
library(tidyverse)
library(plyr)
library(ggstatsplot)
library(reshape2)
library(ggsci)
library(pheatmap)
library(harmony)
library(ggforce)
library(schex)


#macro_all <- readRDS('../../Process/macro_seu_fil.rds')
#macro_UMAP <- readRDS('macro_UMAP.rds')
colorboard <- c('#1A5599','#800020','#008C8C','#FADA5E','#E85827','#B05923','#FF0000')


#seurat流程
##标准化、log转换
macro_all <- NormalizeData(macro_all, normalization.method = "LogNormalize", scale.factor = 10000)

##鉴定高变基因
i <- 2500
macro_all <- FindVariableFeatures(macro_all, selection.method = "vst", nfeatures = i)
hvg <- VariableFeatures(macro_all)
write.csv(hvg,'hvg_first.csv')

##scale
macro_all <- ScaleData(macro_all,features = rownames(macro_all))

##主成分
macro_all <- RunPCA(macro_all)

##harmony
plot1 <-  DimPlot(macro_all,group.by = 'orig.ident')/(DimPlot(macro_all,group.by = 'SOURCE_HOSPITAL')|DimPlot(macro_all,group.by = 'PROCESSING_TYPE'))/
         (DimPlot(macro_all,group.by = 'SINGLECELL_TYPE')|DimPlot(macro_all,group.by = 'TISSUE_PROCESSING_TEAM'))
ggsave('harmony_before.pdf',plot1,width = 15,height = 15)

macro_all <- RunHarmony(macro_all,
                        c('orig.ident','SOURCE_HOSPITAL','PROCESSING_TYPE','SINGLECELL_TYPE','TISSUE_PROCESSING_TEAM'),
                        plot_convergence = TRUE)

plot2 <- DimPlot(macro_all,group.by = 'orig.ident',reduction = 'harmony')/(DimPlot(macro_all,group.by = 'SOURCE_HOSPITAL',reduction = 'harmony')|DimPlot(macro_all,group.by = 'PROCESSING_TYPE',reduction = 'harmony'))/
         (DimPlot(macro_all,group.by = 'SINGLECELL_TYPE',reduction = 'harmony')|DimPlot(macro_all,group.by = 'TISSUE_PROCESSING_TEAM',reduction = 'harmony'))
ggsave('harmony_after.pdf',plot2,width = 15,height = 15)

##选择最佳harmony成分
macro_all <- JackStraw(macro_all, num.replicate = 100)
macro_all <- ScoreJackStraw(macro_all,dims = 1:20)

plot3 <- JackStrawPlot(macro_all, dims = 1:20)#选20个主成分
plot4 <- ElbowPlot(object = macro_all,ndims = 30,reduction = 'harmony')

pdf(paste('JackStrawPlot+ElbowPlot_',i,'.pdf',sep = ''),width = 12,height = 5)
plot3+plot4
dev.off()


#降维
macro_all <- FindNeighbors(macro_all, dims = 1:14,reduction = 'harmony')#调整PCs

##确定resolution 0.1-1.5
#res.used <- seq(0.1,1.5,by=0.1)
#sce <- FindClusters(object = macro_all, verbose = T, resolution = res.used)
#clus.tree.out <- clustree(sce,prefix = 'RNA_snn_res.') +
#  theme(legend.position = "bottom") + 
#  scale_color_brewer(palette = "Set1") +
#  scale_edge_color_continuous(low = "grey80", high = "red")

#pdf(paste('clustertree_',i,'.pdf',sep = ''),width = 30,height = 30)
#clus.tree.out
#dev.off()

##聚类
macro_all<- FindClusters(macro_all, resolution = 0.1)#调整cluster
table(macro_all@meta.data[["seurat_clusters"]])

##确定最佳n.neighbors
#map(seq(5,50,5),function(x) { macro_cluster %>%  RunUMAP(n.neighbors = x,n.epochs=500,dims = 1:14,reduction = 'harmony') %>% DimPlot()}) %>% cowplot::plot_grid(plotlist = .)

##细胞命名
cellnames <- mapvalues(macro_all$seurat_clusters,from = c(0,1,2,3,4,5,6),to = c('M2 like TAM','SPP1+ TAM','FTL+ Macro',
                                                                                'M1 like Macro','MKI67+ Macro','GNLY+ Monolike Macro','STMN1+ Macro'))
macro_all$seurat_clusters3 <- factor(cellnames,levels = c('M2 like TAM','SPP1+ TAM','FTL+ Macro',
                                                          'M1 like Macro','MKI67+ Macro','GNLY+ Monolike Macro','STMN1+ Macro'))

##可视化
macro_UMAP<- RunUMAP(macro_all, dims = 1:14,seed.use = 1,reduction = 'harmony',n.neighbors = 35)
plot5 <- DimPlot(macro_UMAP, reduction = "umap",
                 label = F,
                 pt.size = 0.5,group.by = 'seurat_clusters3',
                 shuffle = TRUE,cols = colorboard,
                 repel = TRUE)+
         ggtitle("Macrophages") +
         theme(text=element_text(size=7,face = 'bold'))

points <-data.frame(macro_UMAP@reductions$umap@cell.embeddings, cluster=macro_UMAP$seurat_clusters3)
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
plot5 <- plot5+geom_mark_ellipse(data=ell, aes(x=X1, y=X2, label=cluster,group=cluster),inherit.aes = F,color=NA,label.fontsize = 5) +NoLegend()
plot5


pdf(paste('cluster_',i,'_0.1_UMAP','.pdf',sep = ''),width = 6,height = 6)
plot5
dev.off()

##差异分析
macro_markers <- FindAllMarkers(object = macro_UMAP, 
                                random.seed = 2,
                                return.thresh = 0.05,
                                logfc.threshold = 0.25,only.pos = TRUE) %>% filter(.,p_val_adj<0.05)
macro_markers <- arrange(macro_markers,cluster,-avg_log2FC)
table(macro_markers$cluster)

macro_markers_fil <- data.frame(NULL)
for (i in unique(macro_markers$gene)) {
  sub1 <- macro_markers[macro_markers$gene==i,]
  sub2 <- macro_markers[macro_markers$avg_log2FC==max(sub1$avg_log2FC),]
  macro_markers_fil <- rbind(macro_markers_fil,sub2)
}
macro_markers_fil <- arrange(macro_markers_fil,cluster,-avg_log2FC)
table(macro_markers_fil$cluster)

datExpr_heatmap <- AverageExpression(macro_UMAP,group.by = 'seurat_clusters3',
                                     features = macro_markers_fil$gene,
                                     slot = 'data')$RNA %>% as.data.frame()
datExpr_heatmap$cluster <- macro_markers_fil$cluster[match(rownames(datExpr_heatmap),macro_markers_fil$gene)]
table(datExpr_heatmap$cluster)

datExpr_heatmap_sub <- data.frame(NULL)
for (i in as.numeric(levels(datExpr_heatmap$cluster))) {
  sub1 <- datExpr_heatmap[datExpr_heatmap$cluster==i,]
  sub2 <- apply(sub1[,1:7], 1, max)
  sub3 <- datExpr_heatmap[datExpr_heatmap$cluster==i,i+1]==sub2
  sub4 <- sub1[sub3,]
  datExpr_heatmap_sub <- rbind(datExpr_heatmap_sub,sub4)
}
datExpr_heatmap <- datExpr_heatmap_sub[,1:7]
datExpr_heatmap <- datExpr_heatmap %>% t() %>% scale() %>% t()

macro_markers_fil <- macro_markers_fil[match(rownames(datExpr_heatmap),macro_markers_fil$gene),]
table(macro_markers_fil$cluster)
write.csv(macro_markers_fil,'macro_markers_first_logfc0.25.csv')

##第一次分群平均表达量热图
plot6 <- pheatmap(datExpr_heatmap,
                  #annotation_col = colanno,
                  show_colnames = TRUE,show_rownames = FALSE,
                  #annotation_colors = coloranno,
                  cellheight = 0.3,cellwidth = 20,
                  fontsize = 6,
                  cluster_rows = FALSE,
                  cluster_cols = FALSE,
                  color = colorRampPalette(c("#3A54A5", "white", "#ED2225"))(100),
                  use_raster = TRUE)
pdf('Heatmap.pdf',height = 50)
plot6
dev.off()

##不同group UMAP图
plot7 <- DimPlot(macro_UMAP, reduction = "umap",label.size = 3,pt.size = 0.5,
                 shuffle = TRUE,cols =c('#6768AB','#D5837B'),group.by = 'MMRStatus')+
         ggtitle("Macrophages")+
         theme(text=element_text(size=7,face = 'bold'))
plot7
ggsave('cluster_group.pdf',plot7,width = 6,height = 6)


#基因表达featureplot
macro_UMAP <- make_hexbin(macro_UMAP,80, dimension_reduction = "umap")

##M1和M2marker基因
plot8 <- plot_hexbin_feature(macro_UMAP,type = 'scale.data',feature = 'F13A1',action = 'mean',title = 'F13A1',xlab = '',ylab = '')+
         scale_fill_viridis_c(option = "inferno")+
         plot_hexbin_feature(macro_UMAP,type = 'scale.data',feature = 'FOLR2',action = 'mean',title = 'FOLR2',xlab = '',ylab = '')+
         scale_fill_viridis_c(option = "inferno")#M2

plot9 <- plot_hexbin_feature(macro_UMAP,type = 'scale.data',feature = 'CXCL10',action = 'mean',title = 'CXCL10',xlab = '',ylab = '')+
         scale_fill_viridis_c(option = "inferno")+
         plot_hexbin_feature(macro_UMAP,type = 'scale.data',feature = 'IL15RA',action = 'mean',title = 'IL15RA',xlab = '',ylab = '')+
         scale_fill_viridis_c(option = "inferno")#M1

ggsave('Featureplot_M1+M2.pdf',plot8/plot9,width = 15,height = 12)

##所有细胞marker基因
plot10 <- plot_hexbin_feature(macro_UMAP,type = 'scale.data',feature = 'CXCL10',action = 'mean',title = 'CXCL10',xlab = '',ylab = '',upper_cutoff = 0.8)+
          theme_void()+theme(legend.position = 'none')+
          scale_fill_viridis_c(option = "inferno")+
          plot_hexbin_feature(macro_UMAP,type = 'scale.data',feature = 'IL15RA',action = 'mean',title = 'IL15RA',xlab = '',ylab = '',upper_cutoff = 0.7)+
          theme_void()+theme(legend.position = 'none')+
          scale_fill_viridis_c(option = "inferno")+#M1
          plot_hexbin_feature(macro_UMAP,type = 'scale.data',feature = 'F13A1',action = 'mean',title = 'F13A1',xlab = '',ylab = '',upper_cutoff = 0.8)+
          theme_void()+theme(legend.position = 'none')+
          scale_fill_viridis_c(option = "inferno")+
          plot_hexbin_feature(macro_UMAP,type = 'scale.data',feature = 'FOLR2',action = 'mean',title = 'FOLR2',xlab = '',ylab = '',upper_cutoff = 0.9)+
          theme_void()+theme(legend.position = 'none')+
          scale_fill_viridis_c(option = "inferno")+#M2
          plot_hexbin_feature(macro_UMAP,type = 'scale.data',feature = 'SPP1',action = 'mean',title = 'SPP1',xlab = '',ylab = '')+
          theme_void()+theme(legend.position = 'none')+
          scale_fill_viridis_c(option = "inferno")+
          plot_hexbin_feature(macro_UMAP,type = 'scale.data',feature = 'FTL',action = 'mean',title = 'FTL',xlab = '',ylab = '',upper_cutoff = 0.9)+
          theme_void()+theme(legend.position = 'none')+
          scale_fill_viridis_c(option = "inferno")+
          plot_hexbin_feature(macro_UMAP,type = 'scale.data',feature = 'MKI67',action = 'mean',title = 'MKI67',xlab = '',ylab = '',upper_cutoff = 0.5)+
          theme_void()+theme(legend.position = 'none')+
          scale_fill_viridis_c(option = "inferno")+
          plot_hexbin_feature(macro_UMAP,type = 'scale.data',feature = 'GNLY',action = 'mean',title = 'GNLY',xlab = '',ylab = '',upper_cutoff = 0.3)+
          theme_void()+theme(legend.position = 'none')+
          scale_fill_viridis_c(option = "inferno")+
          plot_hexbin_feature(macro_UMAP,type = 'scale.data',feature = 'STMN1',action = 'mean',title = 'STMN1',xlab = '',ylab = '')+
          theme_void()+theme(legend.position = 'none')+
          scale_fill_viridis_c(option = "inferno")
plot10
ggsave('Featureplot_all.pdf',plot10,width = 4,height = 4)


#基因表达vlnplot
##M1和M2marker基因
plot11 <- VlnPlot(macro_UMAP,
                 features = c('F13A1','FOLR2','CXCL10','IL15RA'),
                 stack = TRUE,
                 group.by = 'seurat_clusters3',
                 split.by = 'seurat_clusters3',
                 cols = colorboard,
                 same.y.lims = FALSE)+
  theme(axis.text.x = element_text(hjust = 0,size = 8,family = "sans",face = 'bold'),
        axis.text.y = element_text(hjust = 1,vjust = 0,size = 9,family = "sans",face = 'bold'))+
  theme(legend.position = 'none')+
  xlab('')+ylab('')
plot11
ggsave('Vlnplot_M1+M2.pdf',plot11,width = 8,height = 6)#最终认为cluster3是M1

##所有细胞marker基因
plot12 <- VlnPlot(macro_UMAP,
                 features = c('SPP1','FTL','MKI67','GNLY','STMN1'),
                 stack = TRUE,
                 group.by = 'seurat_clusters3',
                 split.by = 'seurat_clusters3',
                 cols = colorboard,
                 same.y.lims = FALSE)+
  theme(axis.text.x = element_text(hjust = 0,size = 8,family = "sans",face = 'bold'),
        axis.text.y = element_text(hjust = 1,vjust = 0,size = 9,family = "sans",face = 'bold'))+
  theme(legend.position = 'none')+
  xlab('')+ylab('')
plot12
ggsave('Vlnplot_all.pdf',plot12,width = 8,height = 6)


#基因表达dotplot
macro_markers_top5 <- c()
for (i in levels(macro_markers_fil$cluster)) {
  sub1 <- macro_markers_fil[macro_markers_fil$cluster==i,]
  sub2 <- slice_max(sub1,order_by = avg_log2FC,n = 5)
  macro_markers_top5 <- rbind(macro_markers_top5,sub2)
}
plot13 <- DotPlot(macro_UMAP,
                  features = macro_markers_top5$gene,
                  group.by  = 'seurat_clusters3',
                  cols = c('#FFFFFF','#AD152A'),
                  col.min = -1,dot.scale = 5)+
          theme(axis.text.x = element_text(hjust = 1,angle = 45,size = 10,family = "sans",face = 'bold'),
                axis.text.y = element_text(hjust = 1,vjust = 0,size = 10,family = "sans",face = 'bold'))+
          coord_flip()+
          xlab('')+ylab('')
plot13
ggsave('Dotplot_all.pdf',plot13,width = 8,height = 8)


#M2(cluster0)提取继续细分（M2a,M2b,M2c,M2d）
saveRDS(macro_UMAP,'macro_UMAP.rds')
