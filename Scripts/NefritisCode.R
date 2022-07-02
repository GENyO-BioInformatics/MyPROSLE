
setwd("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/RData")

library("caret")
library("MLeval")

clin<-readRDS(file = "Nephritis.rds")

load("Clustering.RData")

data<-DATA.Mscore$dataset1$SLE

data<-data[,clin$SampleID]
data<-data.frame("group"=clin$group,t(data))


source("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/scripts/Resources.R")

n.perm.bias<-10

RESULTS<-list()

## Get n train/test sets (to avoid sample selection bias)

colnames(data)<-FixSymbol(vect=colnames(data),symbol=c("/"," "),toSymbol=c(".",""))
  
list.perm.bias<-list()
for(pb in 1:n.perm.bias){
  list.perm.bias[[pb]]<-createDataPartition(y = data$group, p = 0.8, list = FALSE)
}
  
results.perm<-list()
time1<-Sys.time()
for(mdl in 1:length(list.perm.bias)){
    
    # inTrain<-createDataPartition(y = data$group, p = 0.8, list = FALSE)
    ## training<-data[inTrain,]
    ## testing<-data[-inTrain,]
    
    training <- data[ list.perm.bias[[mdl]],]
    testing <- data[-list.perm.bias[[mdl]],]
      
    res.i<-getML.cat(training = training,testing = testing,
                      kfold = 10,repeats = 30,sel.feature = NULL)
      
    results.perm[[mdl]]<-res.i
 
  }
  time2<-Sys.time()
  print(time2-time1)
  

  RESULTS<-results.perm

  save.image("Nephritis.RData")
  
  #########################################
  
  prioML(lista = RESULTS,nameplot = "NeprMods.tiff")
  
  modelNefr<-RESULTS[[7]]$models$rf
  
  train<-RESULTS[[7]]$models$rf$trainingData
  colnames(train)[1]<-"group"
  sel<-ifelse(rownames(data) %in% rownames(train),F,T)
  test<-data[sel,]
  
  modelStats(model.i = modelNefr,test = test)
  
  imp<-varImp(tabModel)
  imp<-imp$importance
  sel<-order(imp$Overall,decreasing=T)
  
  imp<-data.frame(rownames(imp)[sel],
                  imp[sel,"Overall"])
  colnames(imp)<-c("genes","importance")
  
  
  
  
  
  