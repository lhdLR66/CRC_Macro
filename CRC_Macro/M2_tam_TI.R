setwd('E:/work/TI/')

library(dyno)
library(Seurat)
library(tidyverse)
library(Matrix)
library(reshape2)
library(monocle)


#数据读取
macro_umap <- readRDS('../Cluster/second//macro_UMAP.rds')
macro_m2 <- macro_umap[,macro_umap$seurat_clusters2 %in% c(2,4,8,9)]
markers_all <- read.csv('../Cluster/total/macro_markers_all_logfc0.25.csv',header = T,row.names = 1,check.names = F)
markers_m2 <- markers_all[markers_all$cluster %in% c(2,4,8,9),]#所有M2 marker gene并集用于拟时序


#dyno轨迹分析
##构建dataset
dataset_m2 <- wrap_expression(counts = t(macro_m2@assays$RNA@counts[markers_m2$gene,]),
                              expression = t(macro_m2@assays$RNA@data[markers_m2$gene,])) ##wrap_expression要求raw counts和normalised表达为 sparse matrix (dgCMatrix）(列为genes/features，行为细胞)

dataset_m2_all <- wrap_expression(counts = t(macro_m2@assays$RNA@counts),
                                  expression = t(macro_m2@assays$RNA@data))


#添加先验信息(细胞id)，后期可视化可以根据具体的轨迹推断结果进行调整
start_cell <- names(macro_m2$seurat_clusters2[macro_m2$seurat_clusters2==2][1])
groups_id <- data.frame(cell_id=names(macro_m2$seurat_clusters2),group_id=macro_m2$seurat_clusters2)
dataset_m2 <- add_prior_information(dataset_m2,
                                    start_id = start_cell,
                                    start_n = 1,end_n = 1)##有些方法需要先验信息
dataset_m2 <- add_grouping(dataset_m2,grouping = groups_id)

groups <- data.frame(cell_id=names(macro_m2$seurat_clusters2),group_id=macro_m2$seurat_clusters4)

##选择最佳TI方法（shiny）
guidelines_m2 <- guidelines_shiny(dataset_m2)

##选定方法
methods_m2 <- guidelines_m2$methods_selected

##轨迹推断
m2_TI <- infer_trajectory(dataset_m2,method = 'slingshot',give_priors = c('start_id','start_n','end_n'),
                                 verbose = TRUE,seed = 1)
#dimred_m2 <- dyndimred::dimred_umap(dataset_m2$expression)

##添加root
m2_TI <- add_root(m2_TI,root_cell_id = start_cell)

##可视化
plot1 <- 
      plot_dimred(m2_TI,
              expression_source = dataset_m2,
              grouping = groups,
              color_density = 'grouping',
              alpha_cells = 1,size_cells = 1,
  )+ggtitle("Cell cluster")+
    scale_color_manual(values = c('#9A3671','#E1C392','#25BCCD','#95A238'))+
    scale_fill_manual(values = c('#9A3671','#E1C392','#25BCCD','#95A238'))
plot2 <- 
  plot_dimred(m2_TI,alpha_cells = 1,size_cells = 1,
              "pseudotime", pseudotime = calculate_pseudotime(m2_TI),
  )+ggtitle("Pseudotime")+
  theme(legend.position = 'right',legend.key.size = unit(4,'mm'))

ggsave('TI_m2.pdf',plot1,width = 6,height = 4)
ggsave('TI_m2_pseu.pdf',plot2,width = 7,height = 4)

#monocle进行相关基因可视化
##构建monocle对象
phenoData <- macro_m2@meta.data
featureData <- data.frame(row.names = rownames(macro_m2),gene_short_name=rownames(macro_m2))
pd <- new("AnnotatedDataFrame", data = phenoData)
fd <- new("AnnotatedDataFrame", data = featureData)

m2_monocle <- newCellDataSet(macro_m2@assays[["RNA"]]@data,
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit = 0.1,
                        expressionFamily = tobit(Lower = 0.1))


m2_monocle_sub@phenoData@data[["seurat_clusters4"]] <- factor(m2_monocle_sub@phenoData@data[["seurat_clusters4"]],
                                                              levels = c('C1QC+ TAM','M2b like TAM','M2c like TAM','M2a like TAM'))
genetocheck <- c('MSR1','IGF1',
                 'VEGFA','HIF1A')
for (i in 1:4) {
  gene <- genetocheck[i]
  sub1 <- m2_monocle_sub[gene,]
  plot_genes_jitter(sub1,grouping = 'seurat_clusters4',color_by = 'seurat_clusters4',plot_trend = T)+
                           scale_color_manual(values = c('#9A3671','#25BCCD','#95A238','#E1C392'))+
  theme(axis.text.x = element_text(hjust = 1,size = 8,angle = 45,family = "sans",face = 'bold'),
        axis.text.y = element_text(hjust = 1,vjust = 0,size = 9,family = "sans",face = 'bold'))+
    theme(legend.position = 'none')+
    xlab('')+ylab('Expression')
  ggsave(paste(gene,'_polarization_marker.pdf'),width = 8,height = 6)
}






