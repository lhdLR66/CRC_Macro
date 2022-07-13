setwd('E:/work/Hub_gene_futher/Survival/')

library(survival)
library(survminer)
library(tidyverse)


#数据读取
crc_datsur <- read.csv('E:/work/data/bulk/crc_datsur.csv',row.names = 1,header = T,check.names = F)
crc_datsur_hub <- crc_datsur[,c('Patient','OS','OS.time','S100A6','TP53')]
crc_logfpkm <- read.csv('E:/work/data/bulk/crc_datexpr.csv',header = T,row.names = 1,check.names = F)
markers_4<- read.csv('E:/work/Cluster/total/macro_markers_all_logfc0.25.csv',row.names = 1,header = T,check.names = F) %>% filter(.,cluster=='M2c like TAM')
markers_4 <- list(as.vector(markers_4[1:50,]$gene))

icb_group <- read.csv('crc_icbgroup.csv',header = T,row.names = 1,check.names = F)


##TLR8生存分析
cutpoint1<-surv_cutpoint(crc_datsur_hub,time="OS.time",event = "OS",variables = "S100A6")
crc_datsur_sub1 <- surv_categorize(cutpoint1,labels = c('Low','High'))
sfit1<-survfit(Surv(OS.time,OS)~S100A6,data = crc_datsur_sub1)
plot1 <- ggsurvplot(sfit1,
                    legend.title = 'Group',
                    title='S100A6',
                    conf.int = F,
                    pval = T,
                    risk.table = F,
                    palette = c("#5CACEE","#000000"),
                    legend.labs = c("High","Low"),
                    risk.table.col = 'strata',
                    tables.height = 0.2,
                    xlab = "OS Time In Days",
                    ylab = "Overall Survival",
                    pval.method = TRUE,
                    surv.median.line = "hv",
                    ncensor.plot = F)

##MRC1生存分析
cutpoint2<-surv_cutpoint(crc_datsur_hub,time="OS.time",event = "OS",variables = "TP53")
crc_datsur_sub2 <- surv_categorize(cutpoint2,labels = c('Low','High'))
sfit2<-survfit(Surv(OS.time,OS)~TP53,data = crc_datsur_sub2)
plot2 <- ggsurvplot(sfit2,
                    legend.title = 'Group',
                    title='P53',
                    conf.int = F,
                    pval = T,
                    risk.table = F,
                    palette = c("#5CACEE","#000000"),
                    legend.labs = c("High","Low"),
                    risk.table.col = 'strata',
                    tables.height = 0.2,
                    xlab = "OS Time In Days",
                    ylab = "Overall Survival",
                    pval.method = TRUE,
                    surv.median.line = "hv",
                    ncensor.plot = F)

pdf('hub_survplot.pdf',width = 4,height = 4)
plot1
plot2
dev.off()

##浸润生存分析
ssgsea <- gsva(as.matrix(crc_logfpkm), markers_4,
               verbose=TRUE,
               method="ssgsea",
               kcdf="Gaussian",
               min.sz=1,
               max.sz=1000)
crc_datsur_inf <- cbind(crc_datsur_sub1[,1:2],t(ssgsea))
colnames(crc_datsur_inf) <- c('OS.time','OS','infscore')

cutpoint3<-surv_cutpoint(crc_datsur_inf,time="OS.time",event = "OS",variables = "infscore")
crc_datsur_inf <- surv_categorize(cutpoint3,labels = c('Low','High'))
sfit3<-survfit(Surv(OS.time,OS)~infscore,data = crc_datsur_inf)
plot3 <- ggsurvplot(sfit3,
                    legend.title = 'Group',
                    title='Infiltration of M2c like TAMs',
                    conf.int = F,
                    pval = T,
                    risk.table = F,
                    palette = c("#5CACEE","#000000"),
                    legend.labs = c("High","Low"),
                    risk.table.col = 'strata',
                    tables.height = 0.2,
                    xlab = "OS Time In Days",
                    ylab = "Overall Survival",
                    pval.method = TRUE,
                    surv.median.line = "hv",
                    ncensor.plot = F)
pdf('inf_survplot.pdf',width = 4,height = 4)
plot3
dev.off()
