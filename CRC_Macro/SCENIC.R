setwd('E:/work/SCENIC/')

library(SCENIC)
library(Seurat)
library(tidyverse)
library(foreach)
library(pheatmap)
library(stringr)
library(reshape2)
library(limma)


#数据载入
macro_umap <- readRDS('../Cluster/second/macro_UMAP.rds')
macro_sub <- macro_umap[,macro_umap@meta.data[["seurat_clusters2"]]==4]

##细胞信息和表达矩阵构建
cellInfo <- data.frame(macro_sub@meta.data)
datExpr <- GetAssayData(object = macro_sub,slot = "data") %>% as.matrix()


#初始化设置
scenicOptions <- initializeScenic(org = 'hgnc',
                                  dbDir = 'cisTarget_databases',
                                  nCores=1) 

scenicOptions@inputDatasetInfo$cellInfo <- cellInfo
#saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

##将表达矩阵中未在cisTarget database中收录的基因去除
dbFilePath <- getDatabases(scenicOptions)[[1]]
motifRankings <- importRankings(dbFilePath)
genesInDatabase <- colnames(getRanking(motifRankings))
genesLeft_minCells<-rownames(datExpr)
length(genesLeft_minCells)
genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)]
length(genesLeft_minCells_inDatabases)
genesKept <- genesLeft_minCells_inDatabases
datExpr_filtered <- datExpr[genesKept, ]
dim(datExpr_filtered)

#共表达网络构建（GENIE3）
##相关性计算
runCorrelation(datExpr_filtered,scenicOptions)

##鉴定潜在靶点
runGenie3(datExpr_filtered, scenicOptions)#相当费时

##推断共表达模块
scenic_result <- runSCENIC_1_coexNetwork2modules(scenicOptions)

##推断转录调控网络（regulon）
coexMethod <- c("w001","w005","top50","top5perTarget","top10perTarget","top50perTarget")
scenic_result2 <- runSCENIC_2_createRegulons(scenic_result,coexMethods = coexMethod)

##regulon活性评分与可视化##
#3.4_regulonAUC.rds最为重要，有每个细胞的auc值
scenic_result3 <- runSCENIC_3_scoreCells(scenic_result2,exprMat = datExpr_filtered)

##手动调整AUC阈值(可选)
#aucellApp <- plotTsne_AUCellApp(scenic_result3,datExpr_filtered)
#savedSelections <- shiny::runApp(aucellApp)
#scenic_result3@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
#saveRDS(newThresholds, file=getIntName(scenic_result3, "aucell_thresholds"))#会覆盖，有需要重跑result3
#saveRDS(scenic_result3, file="int/scenicOptions.Rds")#会覆盖，有需要重跑result3

##AUC结果二进制转换
scenic_result4 <- runSCENIC_4_aucell_binarize(scenic_result3,exprMat = datExpr_filtered)


#regulon AUC去除extended TF
AUCmatrix <- readRDS("history/int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
  
anno_col <- data.frame(cell=colnames(AUCmatrix),group=macro_sub$MMRStatus)
anno_col <- anno_col[order(anno_col$group),][-1]
AUCmatrix <- AUCmatrix[,rownames(anno_col)]

AUCmatrix_fil <- AUCmatrix[!str_detect(rownames(AUCmatrix),'[extended]'),]
AUCmatrix_fil <- AUCmatrix_fil[,colnames(macro_sub)]


#limma差异分析
##分组信息构建
datgroup <- data.frame(macro_sub$seurat_clusters2,macro_sub$MMRStatus)
colnames(datgroup) <- c('cluster','MMRstatus')

##分组矩阵
design.matrix <- model.matrix(~0+factor(datgroup$MMRstatus))
colnames(design.matrix) <- c('dMMR','pMMR')
rownames(design.matrix) <- rownames(datgroup)

##比较矩阵
contracts.matrix <- makeContrasts(dMMR-pMMR,levels=design.matrix)

##差异分析
fit1 <- lmFit(AUCmatrix_fil, design.matrix)
fit2 <- contrasts.fit(fit1, contracts.matrix)
fit3 <- eBayes(fit2)
summary(decideTests(fit3,adj.P.Val = 0.05))

##结果输出
auc_limma <- topTable(fit3,
                       number = Inf,
                       sort.by = 'logFC') %>% filter(.,adj.P.Val<0.05)
auc_limma <- auc_limma[order(auc_limma$logFC,decreasing = T),]
tf_to_show <- rownames(rbind(head(auc_limma,5),tail(auc_limma,5)))

AUCmatrix_filnal <- AUCmatrix_fil[tf_to_show,]
AUCmatrix_filnal <- AUCmatrix_filnal %>% t() %>% scale() %>% t()
summary(as.numeric(AUCmatrix_filnal))
AUCmatrix_filnal <- AUCmatrix_filnal[,rownames(anno_col)]

AUCmatrix_filnal[AUCmatrix_filnal<(-3)] <- (-3)
AUCmatrix_filnal[AUCmatrix_filnal>3] <- 3

##对患者排序以优化可视化结果
dmmr_auc <- AUCmatrix_filnal[1:5,1:754]
dmmr_auc <- dmmr_auc[,order(colMeans(dmmr_auc),decreasing = T)]

pmmr_auc <- AUCmatrix_filnal[6:10,755:1583]
pmmr_auc <- pmmr_auc[,order(colMeans(pmmr_auc),decreasing = T)]

auc_index <- c(colnames(dmmr_auc),colnames(pmmr_auc))
AUCmatrix_filnal <- AUCmatrix_filnal[,auc_index]

##热图可视化
anno_color <- list(group=c('MMRd'='#6768AB','MMRp'='#D5837B'))

plot1 <- pheatmap(AUCmatrix_filnal,border_color = 'white',
                  annotation_names_row = FALSE,legend = T,annotation_legend = T,
                  main = '1583 dMMR and pMMR M2c like TAMs',
                  show_rownames = T,show_colnames = F,annotation_col = anno_col,
                  annotation_colors = anno_color,
                  fontsize_row = 12,fontsize_number = 18,fontsize_col = 12,
                  cellwidth = 0.2,
                  cellheight = 15,
                  cluster_cols = F,cluster_rows = F,
                  color = colorRampPalette(c("#2771B1", "white", "#BC2A33"))(100),
                  use_raster = TRUE,)
ggsave('TF_regulon_top5.pdf',plot1,width = 8,height = 6)

