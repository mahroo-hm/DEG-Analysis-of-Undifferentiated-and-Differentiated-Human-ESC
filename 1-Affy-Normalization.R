
cat("\014")

library(affy)
library(limma)
#library(simpleaffy)
library(MASS)
library(readxl)
library(WGCNA)


#############################
#Read Data
#############################
#BiocManager::install("hgu133plus2cdf")

setwd("/Users/macbook/Desktop/LBB/R codes/1-Affymetrix/GSE55945/Data")

dat<-ReadAffy()


#############################
#Normalization
#############################
dat2<-rma(dat)
#dat2<-justRMA()
dat2
dat.m<-exprs(dat2)

dat1<-exprs(dat)

#############################
#Computational QC
#############################
setwd("F:\\Bio_WorkShop\\LBB\\1-Affymetrix\\GSE55945\\QC")


#####library(simpleaffy)
#aqc<-qc(dat)
#pdf("Plot1.pdf", width = 30, height = 20)
#plot(aqc)
#dev.off()
####################


#plotMA3by2(dat.m,device="pdf")   # M=log2(R/G)=log2(R)-log2(G)     A=1/2log2(RG)=1/2(log2(R)+log2(G))     M=log ratio   A=mean average


pdf("BoxPlot_B_Norm.pdf", width = 30, height = 15)
boxplot(data.frame(dat1))
dev.off()




pdf("BoxPlot.pdf", width = 30, height = 15)
boxplot(data.frame(dat.m))
dev.off()


pdf("Hclust.pdf", width = 30, height = 15)
dat.dist<-dist(t(dat.m))
hc <- hclust(dat.dist)

## function to set label color
labelCol <- function(x) {
  if (is.leaf(x)) {
    ## fetch label
    label <- attr(x, "label") 
    ## set label color to red for A and B, to blue otherwise
    attr(x, "nodePar") <- list(lab.col=ifelse(grepl("D_",label), "red", "blue"))
  }
  return(x)
}

## apply labelCol on all nodes of the dendrogram
d <- dendrapply(as.dendrogram(hc), labelCol)
op <- par(mfrow =  c(1,1), mar = c(10,4,4,4))
plot(d)

dev.off()



# Nonmetric Multidimensional Scaling (NMDS):  The samples that belong to the same group can be easily distinguished
#from the other groups using these two axes, so the NMDS plot supports our
#view that the quality of the replicates is good.
pdf("NMDS.pdf", width = 30, height = 15)

mds<-isoMDS(dat.dist)
plot(mds$points[,1], mds$points[,2], main="NMDS",xlab="Dimension 1", ylab="Dimension 2",type="n")
text(mds$points[,1], mds$points[,2],rownames(mds$points), cex=1.2)  
dev.off()


#setwd("C:\\Users\\123\\Desktop\\Bio_WorkShop\\Affymetrix\\GSE55945")
#write.csv(dat.m,file = "GSE55945_Normalized.csv")

'*******************************************************************
Gene Map
*******************************************************************'
setwd("F:\\Bio_WorkShop\\LBB\\1-Affymetrix\\GSE55945")
#  Read Dis-regulate genes
y <- as.data.frame(dat.m)
ProbeID=row.names(y)


#  Read Annotation File
x <- read_excel("GPL570.xlsx")
dim(x)

#match probeId to Genes
oo=match(ProbeID,x$ID)
x1=x[oo,]
final=cbind(ProbeID,x1[,2])
final=cbind(final,y[,1:19])
dim(final)

final1=na.omit(final)
dim(final1)

#################################################
# Use the collapseRows for  replicate Gene Names
a1=final1[,-c(1:2)]
a1=apply(a1,2,as.numeric )

a2=final1[,2]
a3=unique(final1[,1])
row.names(a1)=a3
aa=collapseRows(a1,a2,a3,method = "Average")
bb=aa$datETcollapsed
final2=bb

write.table(final2, file = "Full_Normalized.csv",row.names=TRUE, na="",col.names=TRUE, sep=",")
write.table(final2,"Full_Normalized.txt", sep="\t",row.names=T, col.names=T, quote=F)

