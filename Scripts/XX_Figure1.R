

setwd("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/RData")

library("pheatmap")
library("ComplexHeatmap")

load("Heatmap1A.RData")

frec<-NULL
for(i in 1:nrow(Mscore.matrix)){
  frec<-c(frec,sum(ifelse(Mscore.matrix[i,]>=1.65,1,0)))
}
frec<-frec/ncol(Mscore.matrix)

#FrecPlot<-cbind(length(frec):1,frec/ncol(Mscore.matrix))
#FrecPlot<-as.data.frame(FrecPlot)
#FrecPlot<-cbind(FrecPlot,Modules.ann$color)
#colnames(FrecPlot)<-c("x","frecuency","color")


tiff(filename=paste0(opt$resultsPath,sep="","/Clustering_Figure1A__1.tiff"),res = 300,width = 10,height = 3.5,units="in")
p1<-pheatmap(Mscore.matrix[,rownames(ann.pats)],
         show_rownames = F,show_colnames = F, 
         cluster_rows = F,cluster_cols = F,
         breaks=seq(-2.5,2.5,length.out = 100),
         annotation_row = annrow,
         annotation_col = anncol,
         annotation_color=Module.colors,
         gaps_col = colGaps,
         gaps_row = rowGaps,
         border_color = NA,
         color = colorRampPalette(c("#003366","white","#CC3300"))(100))

p2 = rowAnnotation(frecuency = anno_barplot(frec,ylim=c(0,0.7),which = "row",
                                            width=unit(2.5,"cm"),axis=T,
                                            border=FALSE,bar_width=1.2,
                                            gp = gpar(border = NA,
                                                      lty = 0,
                                                      color =Modules.ann$color,
                                                      fill = Modules.ann$color)
                                            ))
p1+p2

dev.off()

###################################

load("switchPlot.RData")










