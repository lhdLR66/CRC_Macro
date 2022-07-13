setwd('E:/work/CellChat/')

library(CellChat)
library(dplyr)
library(Seurat)
library(NMF)
library(ggalluvial)
library(openxlsx)


#数据载入
macro_umap <- readRDS('../Cluster/second/macro_UMAP.rds')
color_10 <- c('#800020','#9A3671','#E1C392','#25BCCD','#95A238','#FADA5E','#B05923','#E85827','#FF0000','#008C8C') 
color_9 <- c('#800020','#9A3671','#25BCCD','#95A238','#FADA5E','#B05923','#E85827','#FF0000','#008C8C') 
source <- c('SPP1+ TAM','C1QC+ TAM','M2a like TAM','M2b like TAM','M1 like Macro','GNLY+ Monolike Macro','MKI67+ Macro','STMN1+ Macro','FTL+ Macro')
target <- c('M2c like TAM')


#导入配体受体相互作用数据库
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)# Show the structure of the database
# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
future::plan("multisession", workers = 4)



####dmmr患者中cellchat分析####
#Part1 CellChat对象的数据输入和处理以及初始化
##创建cellchat对象
macro_dmmr <- macro_umap[,macro_umap@meta.data[["MMRStatus"]]=='MMRd']
dmmr_data.input <- macro_dmmr@assays$RNA@data
dmmr_meta <- macro_dmmr@meta.data
dmmr_cellchat <- createCellChat(object = dmmr_data.input, meta = dmmr_meta, group.by = "seurat_clusters4")

dmmr_cellchat@DB <- CellChatDB

##数据预处理
dmmr_cellchat <- subsetData(dmmr_cellchat)# subset the expression data of signaling genes for saving computation cost.This step is necessary even if using the whole database
dmmr_cellchat <- identifyOverExpressedGenes(dmmr_cellchat)
dmmr_cellchat <- identifyOverExpressedInteractions(dmmr_cellchat)
dmmr_cellchat <- projectData(dmmr_cellchat, PPI.human)# project gene expression data onto PPI network (optional)


#Part2细胞间通讯网络推断
##计算通讯概率并推断网络
dmmr_cellchat@idents <-  droplevels(dmmr_cellchat@idents, exclude = setdiff(levels(dmmr_cellchat@idents),unique(dmmr_cellchat@idents)))##dmmr中没有cluster9，因此要drop 9
dmmr_cellchat <- computeCommunProb(dmmr_cellchat)
dmmr_cellchat <- filterCommunication(dmmr_cellchat, min.cells = 10)

##提取推断得到的细胞通讯网络
dmmr_df.net <- subsetCommunication(dmmr_cellchat)
#df.net <- subsetCommunication(all_cellchat, signaling = 'WNT')#netp返回通路级别数据/sources.use返回信号起点终点/signaling 返回特定信号通路

##在信号通路水平上推断细胞间通讯
dmmr_cellchat <- computeCommunProbPathway(dmmr_cellchat)

##计算整合网络
dmmr_cellchat <- aggregateNet(dmmr_cellchat)
dmmr_groupSize <- as.numeric(table(dmmr_cellchat@idents))

#mat <- all_cellchat@net$weight
#par(mfrow = c(3,4), xpd=TRUE)
#for (i in 1:nrow(mat)) {
#  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#  mat2[i, ] <- mat[i, ]
#  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
#}



####pmmr患者中cellchat分析####
#Part1 CellChat对象的数据输入和处理以及初始化
##创建cellchat对象
macro_pmmr <- macro_umap[,macro_umap@meta.data[["MMRStatus"]]=='MMRp']
pmmr_data.input <- macro_pmmr@assays$RNA@data
pmmr_meta <- macro_pmmr@meta.data
pmmr_cellchat <- createCellChat(object = pmmr_data.input, meta = pmmr_meta, group.by = "seurat_clusters4")

##导入配体受体相互作用数据库
#CellChatDB <- CellChatDB.human
#showDatabaseCategory(CellChatDB)
#dplyr::glimpse(CellChatDB$interaction)# Show the structure of the database
# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
pmmr_cellchat@DB <- CellChatDB

##数据预处理
pmmr_cellchat <- subsetData(pmmr_cellchat)# subset the expression data of signaling genes for saving computation cost
#future::plan("multisession", workers = 4)
pmmr_cellchat <- identifyOverExpressedGenes(pmmr_cellchat)
pmmr_cellchat <- identifyOverExpressedInteractions(pmmr_cellchat)
pmmr_cellchat <- projectData(pmmr_cellchat, PPI.human)# project gene expression data onto PPI network (optional)


#Part2细胞间通讯网络推断
##计算通讯概率并推断网络
pmmr_cellchat <- computeCommunProb(pmmr_cellchat)
pmmr_cellchat <- filterCommunication(pmmr_cellchat, min.cells = 10)

##提取推断得到的细胞通讯网络
pmmr_df.net <- subsetCommunication(pmmr_cellchat)
#df.net <- subsetCommunication(all_cellchat, signaling = 'WNT')#netp返回通路级别数据/sources.use返回信号起点终点/signaling 返回特定信号通路

##在信号通路水平上推断细胞间通讯
pmmr_cellchat <- computeCommunProbPathway(pmmr_cellchat)

##计算整合网络
pmmr_cellchat <- aggregateNet(pmmr_cellchat)
pmmr_groupSize <- as.numeric(table(pmmr_cellchat@idents))

#mat <- all_cellchat@net$weight
#par(mfrow = c(3,4), xpd=TRUE)
#for (i in 1:nrow(mat)) {
#  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
#  mat2[i, ] <- mat[i, ]
#  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
#}



####cellchat对象整合####
#数据整合
group.new <-  levels(pmmr_cellchat@idents)
dmmr_cellchat_lift <- liftCellChat(dmmr_cellchat, group.new)
object.list <- list(dmmr = dmmr_cellchat_lift, pmmr = pmmr_cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
weight.max <- getMaxWeight(object.list, slot.name = c("idents",'net'), attribute = c('idents','counts'))


# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets",
                                       pos.dataset = 'dmmr', features.name = 'dmmr',
                                       only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 0.05)
net <- netMappingDEG(cellchat, features.name = 'dmmr')
net.up <- subsetCommunication(cellchat, net = net, datasets = "dmmr",ligand.logFC = 0.2, receptor.logFC = NULL)
net.down <- subsetCommunication(cellchat, net = net, datasets = "pmmr",ligand.logFC = -0.1, receptor.logFC = -0.1)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

tnf.use <- net.up %>% subset(.$ligand=='TNF'&.$target=='M2c like TAM') %>% .[,'interaction_name'] %>% as.data.frame()
gas.use <- net.down %>% subset(.$ligand=='GAS6'&.$target=='M2c like TAM') %>% .[,'interaction_name'] %>% as.data.frame()
colnames(tnf.use) <- 'interaction_name'
colnames(gas.use) <- 'interaction_name'



####全部结果可视化####
pathways.pro <- c("GAS") 
pathways.apo <- c("TNF") 

##弦图  
plot1 <- netVisual_chord_gene(dmmr_cellchat, big.gap = 20,small.gap = 5,
                              targets.use = target,transparency = 0.5,
                              signaling = pathways.pro,show.legend = T,
                              color.use = color_9)
plot2 <- netVisual_chord_gene(pmmr_cellchat, big.gap = 20,small.gap = 5,
                             targets.use = target,transparency = 0.5,
                             signaling = pathways.pro,show.legend = T,
                             color.use = color_10)
plot3 <- netVisual_chord_gene(dmmr_cellchat, big.gap = 20,small.gap = 5,
                             targets.use = target,transparency = 0.5,
                             signaling = pathways.apo,show.legend = T,
                             color.use = color_9)
plot4 <- netVisual_chord_gene(pmmr_cellchat, big.gap = 20,small.gap = 5,
                             sources.use = source[-3],targets.use = target,transparency = 0.5,
                             signaling = pathways.apo,show.legend = T,
                             color.use = color_10)
pdf('chord.pdf',height = 8,width = 8)
plot1
plot2
plot3
plot4
dev.off()

##热图
plot5 <- netVisual_heatmap(dmmr_cellchat, color.heatmap = 'Reds',color.use = color_9,
                           signaling = pathways.pro,
                           targets.use = target,
                           title.name = 'GAS signaling pathway interaction in dMMR')
plot6 <- netVisual_heatmap(pmmr_cellchat, color.heatmap = 'Reds',color.use = color_10,
                           signaling = pathways.pro,
                           targets.use = target,
                           title.name = 'GAS signaling pathway interaction in pMMR')
plot7 <- netVisual_heatmap(dmmr_cellchat, color.heatmap = 'Reds',color.use = color_9,
                           signaling = pathways.apo,
                           targets.use = target,
                           title.name = 'TNF signaling pathway interaction in dMMR')
plot8 <- netVisual_heatmap(pmmr_cellchat, color.heatmap = 'Reds',color.use = color_10,
                           signaling = pathways.apo,
                           sources.use = source[-3],targets.use = target,
                           title.name = 'TNF signaling pathway interaction in pMMR')
pdf('heatmap.pdf',height = 8,width = 8)
plot5
plot6
plot7
plot8
dev.off()


##点图
plot9 <- netVisual_bubble(cellchat,
                 pairLR.use = tnf.use,
                 sources.use = c(source,target)[-3], targets.use =  target, 
                 comparison = c(1,2),  angle.x = 45, remove.isolate =T,
                 title.name = 'Up-regulated signaling in dMMR')
plot10 <- netVisual_bubble(cellchat,
                           pairLR.use = gas.use,
                           targets.use =  target, 
                           comparison = c(1,2),  angle.x = 45, remove.isolate =T,
                           title.name = 'Up-regulated signaling in pMMR')
ggsave('dot.pdf',plot9+plot10,width = 8,height = 4)


#导出prob数据用于标注在热图上
dmmr_tnf <- dmmr_cellchat@netP[["prob"]][,,'TNF'] 
dmmr_gas <- dmmr_cellchat@netP[["prob"]][,,'GAS'] 
pmmr_tnf <- pmmr_cellchat@netP[["prob"]][,,'TNF'] 
pmmr_gas <- pmmr_cellchat@netP[["prob"]][,,'GAS'] 

list_of_datasets <- list('dmmr_tnf' = dmmr_tnf, 
                         'dmmr_gas' = dmmr_gas,
                         'pmmr_tnf' = pmmr_tnf,
                         'pmmr_gas' = pmmr_gas)
write.xlsx(list_of_datasets, "pathways_prob.xlsx",rowNames=TRUE)
