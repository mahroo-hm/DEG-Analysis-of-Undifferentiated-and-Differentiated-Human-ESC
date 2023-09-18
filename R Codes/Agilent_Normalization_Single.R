'*******************************************************************
                           Clear screen
*******************************************************************'

cat("\014")

library(limma)
library(pracma)
library(genefilter)
library(WGCNA)
library(amap)
library(MASS)
library(readxl)

'*******************************************************************
                          Set Path
*******************************************************************'

setwd("F:\\Bio_WorkShop\\LBB\\2-Agilent\\GSE108214")

'*******************************************************************
                          Read Data
*******************************************************************'

files_case <- dir(pattern="*\\.txt$")
dat <- read.maimages(files_case,source="agilent",green.only=TRUE)

setwd("F:\\Bio_WorkShop\\LBB\\2-Agilent\\QC")
names(dat)

pdf("Plot_Before.pdf", width = 30, height = 20)
#png("plot.png", width = 480, height = 480)
plotDensities(dat)
dev.off()

'*******************************************************************
                           Normalization
*******************************************************************'

dat2<-backgroundCorrect(dat,method="normexp")
dat2 <- normalizeBetweenArrays(dat2, method="quantile")
dim(dat2)

#Calculate Average Duplicate ProbeId 
dat2 <- avereps(dat2, ID=dat2$genes$ProbeName)
dim(dat2)

#Expresisson Data
dat.m<-dat2$E
rownames(dat.m)<-dat2$genes$ProbeName
dat.m=na.omit(dat.m)


'*******************************************************************
                          Quality Control
*******************************************************************'
setwd("C:\\Users\\123\\Desktop\\Bio_WorkShop\\Agilent\\QC")

pdf("plotDensity_After.pdf", width = 30, height = 20)
plotDensities(dat2)
dev.off()

pdf("boxplot.pdf", width = 30, height = 20)
boxplot(data.frame(dat.m))
dev.off()



pdf("histogram.pdf", width = 30, height = 20)
hist(dat2$E,col="red")  
dev.off()



#plotMA3by2(dat.m,device="pdf")   # M=log2(R/G)=log2(R)-log2(G)     A=1/2log2(RG)=1/2(log2(R)+log2(G))     M=log ratio   A=mean average

pdf("hclust.pdf", width = 8, height = 6)
dat.dist<-dist(t(dat.m))
hc <- hclust(dat.dist)

## function to set label color
labelCol <- function(x) {
  if (is.leaf(x)) {
    ## fetch label
    label <- attr(x, "label") 
    ## set label color to red for A and B, to blue otherwise
    attr(x, "nodePar") <- list(lab.col=ifelse(grepl("D",label), "red", "blue"))
  }
  return(x)
}

## apply labelCol on all nodes of the dendrogram
d <- dendrapply(as.dendrogram(hc), labelCol)
op <- par(mfrow =  c(1,1), mar = c(10,4,4,4))
plot(d)

dev.off()



'*******************************************************************
Gene Map
*******************************************************************'

setwd("C:\\Users\\123\\Desktop\\Bio_WorkShop\\Agilent")


#  Read Disregulate_genes
y <- as.data.frame(dat.m)
ProbeID=row.names(y)


#  Read Annotation File
x <- read_excel("GPL6480-9577.xls")
dim(x)

#match probeId to Genes
oo=match(ProbeID,x$ID)
x1=x[oo,]
final=cbind(ProbeID,x1[,2])
final=cbind(final,y[,1:22])
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


