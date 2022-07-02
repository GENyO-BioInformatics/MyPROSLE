
## 

setwd("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/RData")

load("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/RData/Nephritis.RData")

samples<-c(rownames(training),rownames(testing))

rm(list=setdiff(ls(),"samples"))


load("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/RData/Clustering.RData")

Modules.list<-Modules.list[rownames(Modules.ann)]

## Generate random subsets 

SLE<-DATA[[1]]$SLE
SLE<-SLE[,samples]

H<-DATA[[1]]$Healthy
subsets<-list()
for(i in 1:20){
  subsets[[i]]<-sample(1:ncol(H),size = round(ncol(H)*0.10,digits = 0))
}


SUBSETS<-list()
## Calculate M-scores for each subset()
for(sub in 1:length(subsets)){
  SLE.sub <- Get.MyPROSLE(modules.list = Modules.list, Patient=SLE,
                          Healthy = H[,subsets[[sub]]], method="mean")
 SUBSETS[[sub]]<-SLE.sub
}


save.image("pruebaConsistency.RData")


MODELS<-readRDS("modelsMYPROSLE.rds")



library("caret")

x<-t(SUBSETS[[1]])
x<-predict(object = MODELS$Nephritis,newdata = x,type="prob")
x<-as.numeric(x$YES)

for(i in 2:length(SUBSETS)){
  
  tmp<-t(SUBSETS[[i]])
  tmp<-predict(object = MODELS$Nephritis,newdata = tmp,type="prob")
  tmp<-as.numeric(tmp$YES)
  x<-cbind(x,tmp)
  
}


comb<-expand.grid(1:20,1:20)

cors<-NULL
for(i in 1:nrow(comb)){
  res<-cor(x[,comb$Var1[i]],x[,comb$Var2[i]])
  cors<-c(cors,res)
  
}










