
cat('\14')
library("readxl")
library(limma);
library("edgeR")
library("rtracklayer")
library("plyr")

##########################
setwd("F:\\3-Dr.Imani\\2023\\GSE71956\\CD4CD8")
dat=read_excel("Full.xlsx")

col1=dat$ID_REF
dat1=dat[,-1]
dat1=as.data.frame(dat1)
row.names(dat1)=col1



###MAP probeID to GENES
row1=row.names(dat1)
GPL=read_excel("GPL10558.xlsx")

indx=which(row1  %in% GPL$ID)
dat2=cbind(dat1,GPL$ILMN_Gene[indx])
dat3=na.omit(dat2)
col1=colnames(dat3)
col1[50]="symbol"
colnames(dat3)=col1

#################################################
# Use the collapseRows for  replicate Gene Names
library(WGCNA)
a1=dat3[,-c(50)]
a2=row.names(dat3)
a3=dat3[,50]
aa=collapseRows(a1,a3,a2,method = "Average")
bb=aa$datETcollapsed
dat.m=bb

########################

gencode_file = 'gencode.v43.annotation.gtf.gz'
gtf = import.gff(gencode_file, format = 'gtf', genome = 'GRCh38.p13', feature.type = 'exon')

grl = reduce(split(gtf, elementMetadata(gtf)$gene_id))
gene_lengths = ldply(grl, function(x) {
  #sum up the length of individual exons
  return(c('gene_length' = sum(width(x))))
}, .id = 'ensembl_gene_id')


genetype = unique(elementMetadata(gtf)[, c('gene_id', 'gene_type','gene_name')])
colnames(genetype)[1] = 'ensembl_gene_id'
gene_lengths = merge(genetype, gene_lengths)

gene_lengths$ensembl_gene_id = gsub('\\.[0-9]*', '', gene_lengths$ensembl_gene_id)
gene_lengths



###################################################
#Extract Expression and add Gene Type and Gene_Symbol
dat.m=as.data.frame(dat.m)
row1=row.names(dat.m)
df=cbind(gene_lengths@listData$gene_name,gene_lengths@listData$gene_type)
df=as.data.frame(df)
df=df[!duplicated(df$V1),]
df=df[order(df$V1),]
  
indx=which(df$V1 %in% row1 )

GENENAME=df$V1[indx]
indx3=which(df$V1 %in% GENENAME )
indx2=which(row1 %in% GENENAME )

Gene_Symbol  <- df$V1[indx3]
Gene_Symbol <-as.data.frame(Gene_Symbol)
Gene_Type  <- df$V2[indx3]
Gene_Type <-as.data.frame(Gene_Type)

full=cbind(Gene_Type,Gene_Symbol,dat.m[indx2,])
indx4=which(full$Gene_Type %in% "protein_coding")
full_protein=full[indx4,-c(1,2)]

setwd('F:\\3-Dr.Imani\\2023\\GSE71956\\CD4CD8\\IOBR\\protein_coding')
write.table(full, "Final_Expression.csv",row.names=TRUE, na="",col.names=TRUE, sep=",")
write.table(full, "Final_Expression.txt", sep="\t",row.names=T, col.names=T, quote=F) 

write.table(full_protein, "Protein_Codig_Expression.csv",row.names=TRUE, na="",col.names=TRUE, sep=",")
write.table(full_protein, "Protein_Codig_Expression.txt", sep="\t",row.names=T, col.names=T, quote=F) 



