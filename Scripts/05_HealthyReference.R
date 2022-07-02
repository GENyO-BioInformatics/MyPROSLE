##############################
## MyPROSLE 
## R version R version 4.0.4
## 
##############################
## Calculate Mscore without Healthy sample (by patient-paint similarity)

load("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/RData/Clustering.RData")
set.seed(123456788)
source(paste0(opt$scriptPath,"/Resources.R",sep=""))

load.libraries(set="05")

## ................................................................STEP 1
## Preparing reference (Datasets with good number of genes, remove some micro-arrays)


sel.reference<-c(1,4,5,6,8,9) #3: test, #1,4,5,6,8,9: reference
sel.test<-3
genes<-unique(as.character(unlist(Modules.list)))
common.genes<-Reduce(intersect, list(genes,
                                     rownames(DATA$dataset1$SLE),
                                     rownames(DATA$dataset3$SLE),
                                     rownames(DATA$dataset4$SLE),
                                     rownames(DATA$dataset5$SLE),
                                     rownames(DATA$dataset6$SLE),
                                     rownames(DATA$dataset8$SLE),
                                     rownames(DATA$dataset9$SLE)))
rm(genes)

##-------
## Reference gene-expression
Reference<-data.frame(matrix(data=0,ncol=0,nrow=length(common.genes)))
rownames(Reference)<-common.genes

for(i in 1:length(sel.reference)){
  tmp<-DATA[[sel.reference[i]]]$SLE
  tmp<-tmp[common.genes,]
  Reference<-cbind(Reference,tmp)
}
test<-DATA[[sel.test]]$SLE[common.genes,]

Reference.normalized<-apply(Reference,2,NormalizeSamples)
rm(Reference)

##-------
## Reference Mscore
Reference.mscore<-data.frame(matrix(data=0,ncol=0,nrow=nrow(DATA.Mscore$dataset1$SLE)))
rownames(Reference.mscore)<-rownames(DATA.Mscore$dataset1$SLE)
for(i in 1:length(sel.reference)){
  tmp<-DATA.Mscore[[sel.reference[i]]]$SLE
  Reference.mscore<-cbind(Reference.mscore,tmp)
}


##-------
## Get imputed Mscores
Mscores<-data.frame(matrix(data=0,ncol=0,nrow=nrow(DATA.Mscore$dataset1$SLE)))
rownames(Mscores)<-rownames(DATA.Mscore$dataset1$SLE)

for(r in 1:ncol(test)){ ## each patient
  pat<-test[,r]; names(pat)<-rownames(test)
  x<-Get.nearSample(patient=pat,
                    Reference.normalized=Reference.normalized,
                    Reference.mscore=Reference.mscore,
                    k=5)
  Mscores<-cbind(Mscores,x)
}
colnames(Mscores)<-colnames(test)

##-------
## Comparing with real Mscores
test<-DATA.Mscore[[3]]$SLE

real<-cbind(Modules.ann[rownames(test),"color"],test)
colnames(real)[1]<-"color"

imputed<-cbind(Modules.ann[rownames(Mscores),"color"],Mscores)
colnames(imputed)[1]<-"color"

real=melt(real)
imputed=melt(imputed)

plotTest<-cbind(real,imputed$value)
colnames(plotTest)<-c("color","sample","real","imputed")


dev.off()
tiff(filename=paste0(opt$resultsPath,sep="","/HealthyImputation_.tiff"),res = 300,width = 3.7,height = 2.6,units="in")
p1<-ggplot(plotTest,aes(x=real,y=imputed, color=color)) + geom_abline(alpha=0.7) + 
  geom_point(size=0.7) + 
  geom_text(x=-9.5, y=5, label=paste0("pearson = ",sep="",as.character(round(cor(plotTest$real,plotTest$imputed),digits = 2))),color="black") +
  scale_color_manual(values = c(unique(plotTest$color)))+theme_classic()

plot(p1)
dev.off()

## ................................................................STEP 2
## Preparing reference and save  reference necessary objects

sel.reference<-c(1,3,4,5,6,8,9) #3: test, #1,4,5,6,8,9: reference
genes<-unique(as.character(unlist(Modules.list)))
common.genes<-Reduce(intersect, list(genes,
                                     rownames(DATA$dataset1$SLE),
                                     rownames(DATA$dataset3$SLE),
                                     rownames(DATA$dataset4$SLE),
                                     rownames(DATA$dataset5$SLE),
                                     rownames(DATA$dataset6$SLE),
                                     rownames(DATA$dataset8$SLE),
                                     rownames(DATA$dataset9$SLE)))
rm(genes)

##-------
## Reference gene-expression
Reference<-data.frame(matrix(data=0,ncol=0,nrow=length(common.genes)))
rownames(Reference)<-common.genes

for(i in 1:length(sel.reference)){
  tmp<-DATA[[sel.reference[i]]]$SLE
  tmp<-tmp[common.genes,]
  Reference<-cbind(Reference,tmp)
}
cat(paste0("\nSamples in SLE-reference: ",sep="",ncol(Reference)))

Reference.normalized<-apply(Reference,2,NormalizeSamples) ## Matrix with normalized samples by mean/sd (Reference)
rm(Reference)

##-------
## Reference Mscore
Reference.mscore<-data.frame(matrix(data=0,ncol=0,nrow=nrow(DATA.Mscore$dataset1$SLE)))
rownames(Reference.mscore)<-rownames(DATA.Mscore$dataset1$SLE)
for(i in 1:length(sel.reference)){
  tmp<-DATA.Mscore[[sel.reference[i]]]$SLE
  Reference.mscore<-cbind(Reference.mscore,tmp)
}

setwd(opt$rdataPath)

rm(list=setdiff(ls(),c("Reference.mscore","Reference.normalized")))

##-------
## Save references
save.image("References.RData")



## ................................................................ RUN EXAMPLES
## With and without healthy samples, dataset and individual patient

RESULTS1<-Get.MyPROSLE(modules.list = Modules.list,
                       Patient = DATA$dataset3$SLE,
                       Healthy = DATA$dataset3$H,method = "mean",
                       pathReference = paste0(opt$rdataPath,sep="","/References.RData"),
                       nk = 5)

RESULTS2<-Get.MyPROSLE(modules.list = Modules.list,
                       Patient = DATA$dataset3$SLE,
                       Healthy = NULL,method = "mean",
                       pathReference = paste0(opt$rdataPath,sep="","/References.RData"),
                       nk = 5)

Pat = DATA$dataset3$SLE[,1]; names(Pat)<-rownames(DATA$dataset3$SLE)

RESULTS1<-Get.MyPROSLE(modules.list = Modules.list,
                       Patient = Pat,
                       Healthy = DATA$dataset3$H,method = "mean",
                       pathReference = paste0(opt$rdataPath,sep="","/References.RData"),
                       nk = 5)

RESULTS2<-Get.MyPROSLE(modules.list = Modules.list,
                       Patient = Pat,
                       Healthy = NULL,method = "mean",
                       pathReference = paste0(opt$rdataPath,sep="","/References.RData"),
                       nk = 5)


















