setwd('E:/work/Process/')
library(Seurat)
library(tidyverse)
library(Matrix)
library(MASS)
library(cowplot)

#数据处理
##矩阵读取
macro <- readRDS("E:/work/data/scRNA/178341/data_counts/macro_T.rds")

##临床信息
macro_meta <- read.csv('../data/scRNA/178341/data_counts/macro_meta_T.csv',header = T,row.names = 1)

##基因注释信息
gene_anno <- read.csv('../data/bulk/gene_annofile.csv',header = T,row.names = 1)

#构建seurat对象以及数据过滤
##保留protein coding基因
macro <- macro[intersect(rownames(macro),gene_anno$gene_name),]

##删除ERCC基因
erccs <-  grep('^ERCC-', x= rownames(macro),value = T) #并没有
percent_ercc <-   Matrix::colSums(macro[erccs, ])/Matrix::colSums(macro)
fivenum(percent_ercc) 
ercc_index <-  grep(pattern = "^ERCC-", x = rownames(macro), value = FALSE) #获取index
###macro <-  macro[-ercc_index,]
###dim(macro)


#构建seurat对象
macro_seu<- CreateSeuratObject(counts = macro, project = "macrophage",meta.data = macro_meta,min.cells =10 ,min.features = 200)

##统计核糖体基因和线粒体基因
table(grepl("^MT-",rownames(macro_seu)))
table(grepl("^RP[SL][[:digit:]]",rownames(macro_seu)))

##向seurat添加新的meta信息
macro_seu <-  AddMetaData(object= macro_seu, percent_ercc, col.name = "percent.ercc") 
macro_seu[["percent.mt"]] <- PercentageFeatureSet(macro_seu, pattern = "^MT-")
macro_seu[["percent.rp"]] <- PercentageFeatureSet(macro_seu, pattern = "^RP[SL][[:digit:]]")
summary(macro_seu[["nCount_RNA"]])
summary(macro_seu[["nFeature_RNA"]])
summary(macro_seu[["percent.mt"]])
summary(macro_seu[["percent.ercc"]])
summary(macro_seu[["percent.rp"]])

pdf('VlnPlot.pdf',width = 20)
VlnPlot(macro_seu, features = c('nCount_RNA','nFeature_RNA','percent.mt','percent.ercc','percent.rp'))
dev.off()

plot1 <- FeatureScatter(macro_seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(macro_seu, feature1 = "nCount_RNA", feature2 = "percent.rp")
plot3 <- FeatureScatter(macro_seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf('feature_counts.pdf',width = 20)
plot1+plot2+plot3
dev.off()
"
##过滤线粒体和核糖体基因异常表达的细胞
filt_mt_rp <- function(object,data) {
  ###计算 MAD
  mad_value <- mad(object, constant = 1.4826)#为了达到渐进正态
  ###低于median - 3*MAD的值认为是离群值
  ###高于median + 3*MAD的值认为是离群值
  upper_bound <- median(object) + 3*mad_value
  macro_seu_fil <<- subset(x = data, cells =which(object <= upper_bound))
  print(ncol(macro_seu_fil))
}

filt_mt_rp(object = macro_seu$percent.mt,data = macro_seu)
filt_mt_rp(object = macro_seu_fil$percent.rp,data = macro_seu_fil)

##过滤nFeature和nCount异常的细胞
filt_nc_nf <- function(object,data) {
  ###计算 MAD
  mad_value <- mad(object, constant = 1.4826)#为了达到渐进正态
  ###低于median - 3*MAD的值认为是离群值
  ###高于median + 3*MAD的值认为是离群值
  lower_bound <- median(object) - 3*mad_value
  upper_bound <- median(object) + 3*mad_value
  macro_seu_fil <<- subset(x = data, cells =which(object >= lower_bound & object <= upper_bound))
  print(ncol(macro_seu_fil))
}
filt_nc_nf(object = macro_seu_fil$nCount_RNA,data = macro_seu_fil)
filt_nc_nf(object = macro_seu_fil$nFeature_RNA,data = macro_seu_fil)
"

#数据过滤
##查看分布情况
VlnPlot(macro_seu, features = 'nCount_RNA')
nCount_RNA <- as.data.frame(macro_seu$nCount_RNA)
ggplot(nCount_RNA, aes(x = macro_seu$nCount_RNA))+geom_density(color ='black', fill = 'gray')
table(macro_seu@meta.data[["nCount_RNA"]]<40000)

VlnPlot(macro_seu, features = 'nFeature_RNA')
nFeature_RNA <- as.data.frame(macro_seu$nFeature_RNA)
ggplot(nFeature_RNA, aes(x = macro_seu$nFeature_RNA))+geom_density(color ='black', fill = 'gray')
table(macro_seu@meta.data[["nFeature_RNA"]]<5000)

VlnPlot(macro_seu, features = 'percent.mt')
percent.mt <- as.data.frame(macro_seu$percent.mt)
ggplot(percent.mt, aes(x = macro_seu$percent.mt))+geom_density(color ='black', fill = 'gray')
table(macro_seu@meta.data[["percent.mt"]]<20)

VlnPlot(macro_seu, features = 'percent.rp')
percent.rp <- as.data.frame(macro_seu$percent.rp)
ggplot(percent.rp, aes(x = macro_seu$percent.rp))+geom_density(color ='black', fill = 'gray')
table(macro_seu@meta.data[["percent.rp"]]<25)

##过滤矩阵
macro_seu_fil <- subset(macro_seu,nCount_RNA<25000&nFeature_RNA<5000&percent.mt<15&percent.rp<20)

##删除线粒体、核糖体基因
macro_seu_fil <- macro_seu_fil[!grepl("^MT-",rownames(macro_seu_fil)),]
macro_seu_fil <- macro_seu_fil[!grepl("^RP[SL][[:digit:]]",rownames(macro_seu_fil)),]


#效果展示

pdf('vlnplot_fil.pdf',width = 20)
VlnPlot(macro_seu_fil, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.ercc','percent.rp'))
dev.off()

plot4 <- FeatureScatter(macro_seu_fil, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot5 <- FeatureScatter(macro_seu_fil, feature1 = "nCount_RNA", feature2 = "percent.rp")
plot6 <- FeatureScatter(macro_seu_fil, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf('feature_counts_fil.pdf',width = 20)
plot4+plot5+plot6
dev.off()


metadata1 <- macro_seu@meta.data
metadata1$version <- 'untreated'
metadata2 <- macro_seu_fil@meta.data
metadata2$version <- 'treated'
metadata <- rbind(metadata1,metadata2)



version_plot <- function(object) {
  ggplot(aes(color=version, x=object, fill=version),data = metadata) +     
  geom_density(alpha = 0.2) +     
  scale_x_continuous(breaks = seq(0, max(object), by=5)) +     
  theme_classic() 
}
plot7 <- version_plot(metadata$nFeature_RNA)
plot8 <- version_plot(metadata$nCount_RNA)
plot9 <- version_plot(metadata$percent.mt)
plot10 <- version_plot(metadata$percent.rp)

pdf('version.pdf',width = 20,height = 10)
plot_grid(plot7, plot8, plot9, plot10,labels=c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.rp'))
dev.off()

#结果输出
saveRDS(macro_seu_fil,file = 'macro_seu_fil.rds')
