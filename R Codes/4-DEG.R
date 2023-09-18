

cat("\014")
library(limma)
library(readxl)

setwd("F:\\Bio_WorkShop\\LBB\\3-DEG")
dat.m=read.table('Full_Normalized.txt',sep='\t')

Trait=read_excel('Trait.xlsx')

groups<- Trait$Type
groups<-as.factor(groups)
design<-model.matrix(~groups)

fit<-lmFit(dat.m, design)
fit<-eBayes(fit)

topTable(fit, coef=2)
tt<-topTable(fit, coef=2, n=nrow(dat.m))

rn11<-rownames(tt)[tt$adj.P.Val  <=0.05  ]
length(rn11)


rn1<-rownames(tt)[tt$adj.P.Val  <=0.05 & tt$logFC>=1 ]
length(rn1)
rn2<-rownames(tt)[tt$adj.P.Val  <=0.05 & tt$logFC<= -1 ]
length(rn2)
rn=c(rn1,rn2)


dat.s<-dat.m[rn,]
dim(dat.s)


UP=tt[rn1,]
DOWN=tt[rn2,]


'*******************************************************************
                          Save Data dat.c OR dat.f
*******************************************************************'
write.table(tt, "DEG_FullFC.csv",row.names=TRUE, na="",col.names=TRUE, sep=",")
write.table(tt, "DEG_FullFC.txt", sep="\t",row.names=T, col.names=T, quote=F)    

write.table(dat.s, "Disregulate_Genes.csv",row.names=TRUE, na="",col.names=TRUE, sep=",")
write.table(dat.s, "Disregulate_Genes.txt", sep="\t",row.names=T, col.names=T, quote=F)  


write.table(UP, "UP.csv",row.names=TRUE, na="",col.names=TRUE, sep=",")
write.table(UP, "UP.txt", sep="\t",row.names=T, col.names=T, quote=F)  


write.table(DOWN, "DOWN.csv",row.names=TRUE, na="",col.names=TRUE, sep=",")
write.table(DOWN, "DOWN.txt", sep="\t",row.names=T, col.names=T, quote=F)  

