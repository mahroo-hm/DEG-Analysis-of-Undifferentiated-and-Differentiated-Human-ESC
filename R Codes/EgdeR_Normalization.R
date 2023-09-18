#https://web.stanford.edu/class/bios221/labs/rnaseq/lab_4_rnaseq.html#:~:text=edgeR%20is%20concerned%20with%20differential,with%20estimating%20absolute%20expression%20levels.


cat('\14')
library("readxl")
library("edgeR")
#library(dplyr)
############################################
#
#     RNASeq Normalization (edgeR)
#
############################################
setwd('F:\\Bio_WorkShop\\LBB\\4-NGS\\Gene')
Dat = read_excel("Read_Count_Full.xlsx");

#Dat1=distinct(Dat,Gene_Symbol, .keep_all= TRUE) #Remove repeated Gene Symbol
Dat1=Dat[!duplicated(Dat$Gene_Symbol), ]

col1<-Dat1$Gene_Symbol
Dat2<-Dat1[,-1]
Dat2=as.data.frame(Dat2)
row.names(Dat2) <- col1
Expression_Raw=as.matrix(Dat2)

means <- rowMeans(Expression_Raw)
filter <- means >= 1
table(filter)
keepSamples = (filter==TRUE)
geneCountHigh_M <- Expression_Raw[keepSamples,]
dim(geneCountHigh_M)

####################Load Trait
Trait= read_excel("Trait_COAD.xlsx");


col1=colnames(geneCountHigh_M)
indx <-match(Trait$Sample , col1)
Group=Trait[indx,2]

#--------------------------------------------------------#
#-----------------TMM--Normalization---------------------#
#--------------------------------------------------------#
dgellist <- DGEList(counts=geneCountHigh_M,group=factor(Group$Type))
dgellist <- calcNormFactors(dgellist,method = "TMM") #method ="TMM","TMMwsp","RLE","upperquartile","none"
dgellist <- estimateCommonDisp(dgellist)
dgellist <-estimateTagwiseDisp(dgellist,trend="movingave")
Normalexpr <-dgellist$pseudo.counts


lNormexpr <- log(Normalexpr+1)


#######################################
#  QC & Plots
#######################################
pdf(file='Heatmap.pdf',width = 10,height=10)
heatmap(cor(lNormexpr))
dev.off()


# pdf(file='Heatmap.pdf',width = 10,height=10)
# heatmap(cor(lNormexpr[1:10,1:10]))
# dev.off()


# Boxplot Before Normalization
pdf(file='B_boxplot.pdf',width = 10,height=10)
par(mar=c(14,5,1,1))
boxplot(Expression_Raw,las = 2,ylim = c(0, 20),labels=FALSE,col="green")
mtext("Normal samples \n Before normalization",side=2,line = 2)
dev.off()

# Boxplot After Normalization
pdf(file='A_boxplot.pdf',width = 30,height=10)
par(mar=c(14,5,1,1))
boxplot(lNormexpr,las = 2,col="green",ylim = c(0, 20))
mtext("After normalization",side=2,line = 2)
dev.off()


pdf(file='MDS.pdf',width = 30,height=10)
plotMDS(dgellist)
dev.off()
 

write.table(lNormexpr, "Full_Normalized.csv",row.names=TRUE, na="",col.names=TRUE, sep=",")
write.table(lNormexpr, "Full_Normalized.txt", sep="\t",row.names=T, col.names=T, quote=F)  


############################################
#
#     DEG Analysis
#
############################################ 
dge <- DGEList(counts=lNormexpr,group=factor(Group$Type))
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)

#######################################
# 1-DEG ( exact test )
#######################################

 et <- exactTest(dge, pair=c(1,2))
 DEG=topTags(et,n=nrow(et))
 DEG=DEG$table
 write.table(DEG, "DEG_FULL.csv",row.names=TRUE, na="",col.names=TRUE, sep=",")
 write.table(DEG, "DEG_FULL.txt", sep="\t",row.names=T, col.names=T, quote=F)    
 
#######################################
# 2-DEG ( quasi-likelihood (QL) F-test or likelihood ratio test )
#######################################
  group <- factor(Group$Type)
  design <- model.matrix(~group)
  fit <- glmQLFit(dge, design)
  qlf.2vs1 <- glmQLFTest(fit, coef=2)
  DEG1=topTags(qlf.2vs1,sort.by="PValue",n=nrow(qlf.2vs1))
  DEG1=DEG1$table
  write.table(DEG1, "DEG_FULL.csv",row.names=TRUE, na="",col.names=TRUE, sep=",")
  write.table(DEG1, "DEG_FULL.txt", sep="\t",row.names=T, col.names=T, quote=F)    
  
#######################################
# 3-DEG  ( Generalized linear models (GLMs) )
#######################################
  group <- factor(Group$Type)
  design <- model.matrix(~group)
  fit <- glmFit(dge, design)
  lrt.2vs1 <- glmLRT(fit, coef=2)
  DEG2=topTags(lrt.2vs1,sort.by="PValue",n=nrow(lrt.2vs1))
  DEG2=DEG2$table
  write.table(DEG2, "DEG_FULL.csv",row.names=TRUE, na="",col.names=TRUE, sep=",")
  write.table(DEG2, "DEG_FULL.txt", sep="\t",row.names=T, col.names=T, quote=F)    
  









