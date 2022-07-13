setwd('E:/work/Hallmark//')
library(GSVA)
library(reshape2)
library(Seurat)
library(tidyverse)
library(limma)
library(clusterProfiler)


#数据读入与处理
##表达矩阵
#macro_UMAP <- readRDS('../Cluster/second/macro_UMAP.rds')
macro_sub <- macro_UMAP[,macro_UMAP$seurat_clusters2==4]
datExpr <- macro_sub@assays[["RNA"]]@data %>% as.matrix()

##gmt文件构建
hallmark <- read.gmt('../data/scRNA/178341/rawdata/h.all.v7.5.1.symbols.gmt')
hallmark <-  split(hallmark$gene, hallmark$term) 


#GSVA
gsva.es <- gsva(datExpr, hallmark,
                verbose=TRUE,
                method="gsva",
                kcdf="Gaussian",
                min.sz=1,
                max.sz=500)


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
fit1 <- lmFit(gsva.es, design.matrix)
fit2 <- contrasts.fit(fit1, contracts.matrix)
fit3 <- eBayes(fit2)
summary(decideTests(fit3,p.value = 0.05))

##结果输出
GSVA_limma <- topTable(fit3,
                       number = Inf,
                       sort.by = 'logFC')


##构建画图数据
datBar <- data.frame(GSVA_limma$t,GSVA_limma$P.Value,rownames(GSVA_limma))
colnames(datBar) <- c('t','padj','pathway')
rownames(datBar) <- rownames(GSVA_limma)
datBar <- datBar %>% mutate(.,group=as.factor(ifelse(padj>=0.05,'ns',ifelse(t>0,'enriched in dMMR','enriched in pMMR'))))
datBar <- datBar %>% mutate(.,pathway=factor(pathway,levels = pathway[order(t,decreasing = F)]))

plot1 <- ggplot(data=datBar,mapping=aes(x=t,y=pathway,fill=group))+
         geom_bar(stat = 'identity')+
         scale_fill_manual(values = c('#6768AB','#D5837B','#CCCCCC'))+
         geom_vline(xintercept = c(-2,2),linetype="dashed",color='white')+
         theme_test()+
         ggtitle('M2c like TAM, dMMR vs pMMR')+
         xlab('t value of GSVA score,dMMR versus pMMR')+ylab('')+
        theme(axis.text = element_text(face='bold',size = 10,angle = 0,family = 'sans',colour = 'black'),
              axis.text.x =element_text(vjust = 0.5,angle = 0),
              plot.title = element_text(hjust = 0.5,size = 6),
              legend.position = '')
plot1
ggsave('hallmark_cluster4.pdf',plot1,width = 10,height = 8)
write.csv(gsva.es['HALLMARK_APOPTOSIS',],'HALLMARK_APOPTOSIS_score.csv')
