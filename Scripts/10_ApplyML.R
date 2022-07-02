##############################
## MyPROSLE 
## R version R version 4.0.4
## 
##############################
## PredictiveModels - Get predictive models


load("D:/DATA/WORK/Toro.et.al.MyPROSLE/RData/ClinicalVariables.RData")

source("D:/DATA/WORK/Toro.et.al.MyPROSLE/Scripts/Resource_models.R")

##-------------------------------------------------- (STEP 1)
## Get best models for each variable

n.perm.bias<-10

RESULTS<-list()
for(var.i in 1:length(variable.list)){
  
  print(names(variable.list)[var.i])
  
  ## Get n train/test sets (to avoid sample selection bias)
  data<-variable.list[[var.i]]$mscore
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
    
    if(variable.list[[var.i]]$variable.type=="categoric"){
      
      res.i<-getML.cat(training = training,testing = testing,
                       kfold = 10,repeats = 10,sel.feature = NULL)
      
      results.perm[[mdl]]<-res.i
      
    }else{
      
      res.i<-getML.num(training = training,testing = testing,
                       kfold = 10,repeats = 10,sel.feature = NULL)
      results.perm[[mdl]]<-res.i
      
    }

  }
  time2<-Sys.time()
  print(time2-time1)
  
  RESULTS[[var.i]]<-results.perm
  names(RESULTS)[var.i]<-names(variable.list)[var.i]
  
  saveRDS(RESULTS,"MscoreModels.rds")
  
}
setwd("D:/DATA/WORK/Toro.et.al.MyPROSLE")

save.image("Models2022.RData")


##-------------------------------------------------- (STEP 2)
## Select models

setwd("D:/DATA/WORK/Toro.et.al.MyPROSLE")

#readRDS("MscoreModels.rds")

memory.limit(size = 10000000000)


load("Models2022.RData")

source("D:/DATA/WORK/Toro.et.al.MyPROSLE/Scripts/Resource_models.R")



