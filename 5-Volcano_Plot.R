
library(calibrate)
library(readxl)

cat('\014')
setwd("F:\\Bio_WorkShop\\LBB\\3-DEG")
res = read.csv("DEG_FullFC.csv");

# Make a basic volcano plot
jpeg("Volcano_Plot.jpg", width = 9, height = 5, units = 'in',  res = 500)

with(res, plot(logFC, -log10(adj.P.Val), pch=20, main="Volcano plot", xlim=c(-5,7), cex=.8))
with(subset(res, adj.P.Val<.05), points(logFC, -log10(adj.P.Val), pch=20, col="green", cex=0.8))
with(subset(res, adj.P.Val<.05 & logFC >1), points(logFC, -log10(adj.P.Val), pch=20, col="red", cex=0.8))
with(subset(res, adj.P.Val<.05 & logFC< -1), points(logFC, -log10(adj.P.Val), pch=20, col="blue", cex=0.8))

legend("bottomright", legend=c("logFC > 1", "logFC < -1","no-significant"),
       col=c("red", "blue","black"),  pch=20, cex=0.8,text.width = 2,text.font = 2, x.intersp=0.8, y.intersp=0.8)

dev.off()

