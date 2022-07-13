setwd('E:/work/Cluster/total/')
library(Seurat)
library(tidyverse)
library(plyr)
library(ggsci)
library(reshape2)
library(ggstatsplot)
library(psych)
library(pheatmap)
library(ggcorrplot)

colorboard <- c('#800020','#9A3671','#E1C392','#25BCCD','#95A238','#FADA5E','#B05923','#E85827','#FF0000','#008C8C')
#macro_UMAP <- readRDS('../second/macro_UMAP.rds')

#总体情况
##总体分群差异分析
macro_markers_all <- data.frame(NULL)
for (i in levels(macro_UMAP$seurat_clusters4)) {
  sub1 <- FindMarkers(macro_UMAP,ident.1 = i,random.seed = 1,
                      only.pos = TRUE,group.by = 'seurat_clusters4',
                      logfc.threshold = 0.25) %>% filter(.,p_val_adj<0.05)
  sub1$cluster <- i
  sub1$gene <- rownames(sub1)
  macro_markers_all <- rbind(macro_markers_all,sub1)
}
macro_markers_all <- arrange(macro_markers_all,cluster,-avg_log2FC)
table(macro_markers_all$cluster)

macro_markers_fil <- data.frame(NULL)
for (i in unique(macro_markers_all$gene)) {
  sub1 <- macro_markers_all[macro_markers_all$gene==i,]
  sub2 <- macro_markers_all[macro_markers_all$avg_log2FC==max(sub1$avg_log2FC),]
  macro_markers_fil <- rbind(macro_markers_fil,sub2)
}
macro_markers_fil <- arrange(macro_markers_fil,cluster,-avg_log2FC)

datExpr_heatmap_all <- AverageExpression(macro_UMAP,group.by = 'seurat_clusters4',
                                         features = macro_markers_fil$gene,
                                         slot = 'data')$RNA %>% as.data.frame()

datExpr_heatmap_all$cluster <- macro_markers_fil$cluster[match(rownames(datExpr_heatmap_all),macro_markers_fil$gene)]
table(datExpr_heatmap_all$cluster)

datExpr_heatmap_sub_all <- data.frame(NULL)
for (i in levels(macro_UMAP$seurat_clusters4)) {
  sub1 <- datExpr_heatmap_all[datExpr_heatmap_all$cluster==i,]
  sub2 <- apply(sub1[,1:10], 1, max)
  sub3 <- datExpr_heatmap_all[datExpr_heatmap_all$cluster==i,i]==sub2
  sub4 <- sub1[sub3,]
  datExpr_heatmap_sub_all <- rbind(datExpr_heatmap_sub_all,sub4)
}
datExpr_heatmap_all <- datExpr_heatmap_sub_all[,1:10]
datExpr_heatmap_all <- datExpr_heatmap_all %>% t() %>% scale() %>% t()

macro_markers_fil <- macro_markers_fil[match(rownames(datExpr_heatmap_all),macro_markers_fil$gene),]
table(macro_markers_fil$cluster)
write.csv(macro_markers_fil,'macro_markers_all_logfc0.25.csv')

##平均表达量(top35)热图
markers_all_top35 <- data.frame(NULL)
for (i in levels(as.factor(macro_markers_fil$cluster))) {
  sub1 <- macro_markers_fil[macro_markers_fil$cluster==i,]
  sub1 <- slice_max(sub1,n = 35,order_by = avg_log2FC)
  markers_all_top35 <- rbind(markers_all_top35,sub1)
}
table(markers_all_top35$cluster)
datExpr_heatmap_top35 <- datExpr_heatmap_all[markers_all_top35$gene,]

plot1 <- pheatmap(datExpr_heatmap_top35,
                  #annotation_col = colanno,
                  show_colnames = TRUE,show_rownames = FALSE,
                  #annotation_colors = coloranno,
                  cellheight = 0.5,cellwidth = 10,
                  fontsize = 6,
                  cluster_rows = FALSE,
                  cluster_cols = FALSE,
                  color = colorRampPalette(c("#3A54A5", "white", "#ED2225"))(100),
                  use_raster = TRUE)
ggsave('Heatmap_all.pdf',plot1,height = 10,width = 8)


#全局marker差异基因可视化
##all cell Dotplot
macro_markers_top5 <- c()
for (i in unique(macro_markers_fil$cluster)) {
  sub1 <- macro_markers_fil[macro_markers_fil$cluster==i,]
  sub2 <- slice_max(sub1,order_by = avg_log2FC,n = 5)
  macro_markers_top5 <- rbind(macro_markers_top5,sub2)
}

plot4 <- DotPlot(macro_UMAP,
                 features = macro_markers_top5$gene,
                 group.by  = 'seurat_clusters4',
                 cols = c('#FFFFFF','#AD152A'),
                 col.min = 0,col.max = 3,dot.scale = 0.5,dot.min = 0.1,)+
          theme(axis.text.x = element_text(hjust = 1,vjust = 1,angle = 45,size = 3.5,family = "sans",face = 'bold'),
                axis.text.y = element_text(hjust = 1,vjust = 0.5,angle = 0,size = 3.5,family = "sans",face = 'bold'),
                legend.key.size = unit(1,'mm'),legend.position = 'top',
                legend.text = element_text(size = 3,family = "sans",face = 'bold'),
                legend.title = element_text(vjust = 1,size = 3,family = "sans",face = 'bold'))+
          xlab('')+ylab('')
plot4
ggsave('Dotplot_all.pdf',plot4,width = 3.2,height = 2)

plot5 <- DimPlot(macro_UMAP, reduction = "umap",
                 label = F,
                 pt.size = 0.5,group.by = 'seurat_clusters4',
                 shuffle = TRUE,cols = colorboard,
                 repel = TRUE)+
         ggtitle("Macrophages") +
         theme(text=element_text(size=7,face = 'bold'))
plot5
ggsave('all_macro_UMAP.pdf',plot5,width = 7,height = 6)
