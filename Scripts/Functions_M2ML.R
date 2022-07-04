
require("caret")
require("mlbench")
require("pROC")
require("rpart")
require("randomForest")
require("nnet")
require("caretEnsemble")
require("MLeval")
require("pROC")
require("ROCR")
require("stringr")

##-------
## Fix gramatical symbols 
FixSymbol<-function(vect,
                    symbol=c("/"," "), 
                    toSymbol=c(".","")){
  
  for(i in 1:length(symbol)){
    vect<-gsub(pattern = symbol[i],replacement = toSymbol[i],vect)
  }
  return(vect)
}

#-------
## Function to retrieve Mscore for one patient and one path against healthy control distribution
Get.path.mscore<-function(path,            ## Genes from one path (module)
                          Patient,         ## Numeric vector of gene expression values of a patient (with Gene_Symbol as names)
                          Healthy,         ## Gene expression matrix of healthy samples
                          method = "mean") ## proportion (N? of high dysregulated genes/ N? genes in the paths), mean (mean os zscores for genes of the path)
{
  genes.path<-Reduce(intersect, list(rownames(Healthy),names(Patient),as.character(unlist(path))))
  path.mscore<-0
  
  if(length(genes.path)>2){ ## 3 genes at least in each path
    tmpRef<-Healthy[genes.path,]
    tmpPat<-Patient[genes.path]
    
    # Zscore calculation by module () 
    sel<-ifelse(apply(tmpRef,1,sd)==0,F,T)
    
    if("TRUE" %in% sel){
      if(table(sel)[["TRUE"]]>2){ ## only in modules with more than two correctly meassured genes
        tmpRef<-tmpRef[sel,]
        tmpPat<-tmpPat[sel]
        
        Zscore.genes<-(tmpPat-apply(tmpRef,1,mean))/apply(tmpRef,1,sd)
        
        if(method=="proportion"){
          genes.up<-(length(Zscore.genes[Zscore.genes>=1.65]))/length(Zscore.genes)
          genes.down<-length(Zscore.genes[Zscore.genes<=(-1.65)])/length(Zscore.genes)
          
          if(genes.up>=genes.down){
            path.mscore<-genes.up
          }else{
            path.mscore<-(-genes.down)
          }
        }else if(method=="mean"){
          path.mscore<-mean(Zscore.genes,na.rm=T)
        }
      }
    }
  }
  names(path.mscore)<-names(path)
  return(path.mscore)
}

#-------
## Get Personalyzed dysregulated profiles for SLE patients with respect to healthy distribution
Get.MyPROSLE<-function(modules.list,       ## List of paths (modules) with their genes
                       Patient,            ## Numeric vector of gene expression values of a patient (with Gene_Symbol as names) or dataframe with multiple patients
                       Healthy,            ## Gene expression matrix of healthy samples
                       method="mean",      ## proportion (N? of high dysregulated genes/ N? genes in the paths), mean (mean os zscores for genes of the path)
                       pathReference=NULL, ## Only for Healthy=NULL, path for Reference.RData    
                       nk=5)               ## Only for Healthy=NULL, Number of most similar samples to impute mscore
{
  
  require("parallel")
  res<-NULL
  cl<-makeCluster(detectCores())
  
  if(is.vector(Patient)){ ## Only one patient
    cat("\nMeasuring Mscores for one sample...")
    
    if(is.null(Healthy)){
      cat("\nNo healthy samples detected, calculating Mscores by similarity between SLE-samples")
      
      load(pathReference)
      
      genes<-intersect(names(Patient),rownames(Reference.normalized))
      Patient<-Patient[genes]
      
      res<-Get.nearSample(patient=Patient,
                          Reference.normalized=Reference.normalized[genes,],
                          Reference.mscore=Reference.mscore,
                          k=nk)
      res<-as.data.frame(res);colnames(res)<-"Mscores"  
      
    }else{
      cat("\nHealthy samples detected, calculating Mscores by zscores")
      res<-parLapply(cl, modules.list,Healthy=Healthy,Patient=Patient,method=method, Get.path.mscore)
      res<-as.data.frame(do.call("rbind",res)); colnames(res)<-"Mscores"      
    }
    
  }else{ ## Several patients
    if(is.null(Healthy)){
      cat("\nNo healthy samples detected,calculating Mscores by similarity between SLE-samples for ",sep="",ncol(Patient)," patients...")
      
      load(pathReference)
      
      genes<-intersect(rownames(Patient),rownames(Reference.normalized))
      Patient<-Patient[genes,]
      
      Mscores<-data.frame(matrix(data=0,ncol=0,nrow=nrow(Reference.mscore)))
      rownames(Mscores)<-rownames(Reference.mscore)
      
      pb = txtProgressBar(min = 0, max = ncol(Patient), initial = 0)
      for(r in 1:ncol(Patient)){ ## each patient
        setTxtProgressBar(pb,r)
        pat<-Patient[,r]; names(pat)<-rownames(Patient)
        
        x<-Get.nearSample(patient=pat,
                          Reference.normalized=Reference.normalized[genes,],
                          Reference.mscore=Reference.mscore,
                          k=nk)
        Mscores<-cbind(Mscores,x)
      }
      colnames(Mscores)<-colnames(Patient)
      res<-Mscores
      
    }else{
      
      cat("\nHealthy samples detected, measuring Mscores by zscore for ",sep="",ncol(Patient)," patients...")
      
      res<-data.frame(matrix(nrow=length(modules.list),ncol=ncol(Patient)))
      colnames(res) = colnames(Patient); rownames(res)<-names(modules.list)
      cat("\n")
      pb = txtProgressBar(min = 0, max = ncol(Patient), initial = 0) 
      for(i in 1:ncol(Patient)){
        setTxtProgressBar(pb,i)
        pat = Patient[,i]
        names(pat) = rownames(Patient)
        res.i<-parLapply(cl, modules.list,Healthy=Healthy,Patient=pat,method=method, Get.path.mscore)
        res.i<-as.data.frame(do.call("rbind",res.i)); colnames(res.i)<-colnames(Patient)[i]
        
        res[,i]<-res.i
      }
    }
  }
  stopCluster(cl)
  return(res)
}

##-------
## Normalize by zscore
NormalizeSamples<-function(x){
  xm<-mean(x,na.rm=T)
  xsd<-sd(x,na.rm=T)
  for(i in 1:length(x)){
    if(!is.na(x[i])){
      x[i]<-(x[i]-xm)/xsd
    }
  }
  return(x)
}

##-------
## Get near samples 
Get.nearSample<-function(patient, ## 
                         Reference.normalized,
                         Reference.mscore,
                         k=5){
  require("class")
  
  patient<-patient[intersect(names(patient),rownames(Reference.normalized))]
  patient<-NormalizeSamples(patient)
  
  res<-NULL
  for(i in 1:ncol(Reference.normalized)){
    res<- c(res,dist(rbind(patient, Reference.normalized[,i])))
  }
  names(res)<-colnames(Reference.normalized)
  res<-res[order(res,decreasing = F)]
  res<-names(res)[1:k]
  
  ##-------
  tmp.mscore<-Reference.mscore[,res]
  tmp.mscore<-apply(tmp.mscore,1,mean)
  
  return(tmp.mscore)
}


##----------
## M2ML Function
M2ML<-function(exp.data,   
               geneset,
               metadata,
               var2predict,
               var.type="cat",
               outerfolds=10,
               prob=0.8,
               innerfold=10,
               repeats=10,
               ref.panel=NULL,
               seed=123456788,
               nk=5
){
  
  require("caret")
  require("mlbench")
  require("pROC")
  require("rpart")
  require("randomForest")
  require("nnet")
  require("caretEnsemble")
  require("MLeval")
  require("pROC")
  require("ROCR")
  
  
  ##--------- (Step1)
  ## Calculate M-scores
  sel<-as.numeric(as.character(exp.data$class))
  if(length(sel[sel==0])>=2 & is.null(ref.panel)){ # with controls
    healthy<-exp.data$exp[,ifelse(sel==0,T,F)]
  }else{
    healthy<-NULL
    if(is.null(ref.panel)){
      error<-"Not enough healthy samples or not reference pannel selected"
      return(error)
    }
  }
  disease<-exp.data$exp[,ifelse(sel==1,T,F)]
  
  Mscores<-Get.MyPROSLE(modules.list = geneset, 
                        Patient=disease,
                        Healthy = healthy, 
                        method="mean",
                        pathReference=ref.panel,
                        nk=nk)
  
  
  ##--------- (Step2)
  ## Get predictive model
  
  if(var2predict %in% colnames(metadata)){
    
    samples<-intersect(rownames(metadata),colnames(Mscores))
    Mscores<-Mscores[,samples]
    metadata<-metadata[samples,]
    
    M<-data.frame("group"=metadata[,var2predict],as.data.frame(t(Mscores)))
    rownames(M)<-FixSymbol(vect = rownames(M),symbol = c("/"," ","-"),toSymbol = c(".",".","."))
    
    sel<-ifelse(is.na(M$group),F,T)
    M<-M[sel,]
    
    result.ML<-neastkfoldML(data = M,
                            var.type=var.type,
                            outerfolds = outerfolds,
                            kfold = innerfold,
                            prob=prob,
                            repeats = repeats)
    
    return(result.ML)
    
  }else{
    error<-"Variable to predict not found in metadata"
    return(error)
  }
  
}



##--------
## Get performance results for different algorithms
neastkfoldML<-function(data,
                 outerfolds,
                 kfold,
                 repeats,
                 prob=0.8,
                 var.type){
  
  list.perm.bias<-list()
  for(pb in 1:outerfolds){
    list.perm.bias[[pb]]<-createDataPartition(y = data$group, p = prob, list = FALSE)
  }
  
  ##-------------- (Step1)
  ## Outer fold
  outerML<-list()
  for(i in 1:length(list.perm.bias)){ 
    
    ##-------
    ## Train/ test
    training <- data[list.perm.bias[[i]],]
    testing <- data[-list.perm.bias[[i]],]
    
    my_control <- trainControl(
      method="repeatedcv",
      number=kfold, ## inner fold
      savePredictions="final",
      classProbs=TRUE,
      index=createResample(training$group, kfold),
      repeats=repeats)
    
    
    ##-------
    ## Get models (for tunning parameters)
    
    ## num or cat variable
    if(var.type=="cat"){
      
      model_list <- caretList(
        group~., data=training,
        trControl=my_control,
        methodList=c("glm","lda"),
        tuneList=list(
          xgbTree=caretModelSpec(method="xgbTree", tuneGrid=expand.grid(max_depth = c(2, 3, 4, 5, 6, 8, 10),nrounds = 50,eta = c(0.01,0.05,0.1),gamma = c(0,1),colsample_bytree=c(0.1, 0.4), min_child_weight=c(1, 10), subsample=c(0.5, 1))),
          rf=caretModelSpec(method="rf", tuneGrid=data.frame(.mtry=seq(2,round(sqrt(ncol(training)-1)),1))),
          knn=caretModelSpec(method="knn", tuneGrid=data.frame(.k=seq(2,50,2))),
          svmLinear=caretModelSpec(method="svmLinear", tuneLength = 15),
          svmRadial=caretModelSpec(method="svmRadial", tuneLength = 15),
          nnet=caretModelSpec(method="nnet", tuneGrid=expand.grid(.size = seq(from = 1, to = 4, by = 1),.decay = seq(from = 0.1, to = 0.5, by = 0.1))),
          nb=caretModelSpec(method="nb", tuneGrid=expand.grid(fL=c(0,0.5,1.0), usekernel = TRUE, adjust=c(0.5,1.0)))))
      
      ##--------
      ## Get performance
      model_pred<-list()
      cm<-list()
      for(mp in 1:length(model_list)){ ## Get prediction in outer test (for algorithm selection)
        
        p <- predict(model_list[[mp]], newdata=testing,type="prob")
        if(class(p)=="numeric"){
          cm.model<-confusionMatrix(as.factor(ifelse(p<=0.5,"YES","NO")),as.factor(testing$group))
          pred<-data.frame("NO"=p,"YES"=1-p,"obs"=testing$group)
        }else{ ##data.frame
          cm.model<-confusionMatrix(as.factor(ifelse(p$NO<=0.5,"YES","NO")),as.factor(testing$group))
          pred<-data.frame("NO"=p$NO,"YES"=1-p$NO,"obs"=testing$group)
        }
        pred<-na.omit(pred)
        model_pred[[mp]]<-pred
        cm[[mp]]<-cm.model
      }
      names(model_pred)<-names(model_list)
      
      res<-list(model_list,model_pred,cm)
      names(res)<-c("models","preds","cm")
      
      outerML[[i]]<-res
      
    }else { ## num
      
      model_list <- caretList(
        group~., data=training,
        trControl=my_control,
        methodList=c("glm","lm"),
        tuneList=list(
          xgbTree=caretModelSpec(method="xgbTree", tuneGrid=expand.grid(max_depth = c(2, 3, 4, 5, 6, 8, 10),nrounds = 50,eta = c(0.01,0.05,0.1),gamma = c(0,1),colsample_bytree=c(0.1, 0.4), min_child_weight=c(1, 10), subsample=c(0.5, 1))),
          rf=caretModelSpec(method="rf", tuneGrid=data.frame(.mtry=seq(2,round(sqrt(ncol(training)-1)),1))),
          knn=caretModelSpec(method="knn", tuneGrid=data.frame(.k=seq(2,50,2))),
          svmLinear=caretModelSpec(method="svmLinear", tuneLength = 15),
          svmRadial=caretModelSpec(method="svmRadial", tuneLength = 15),
          nnet=caretModelSpec(method="nnet", tuneGrid=expand.grid(.size = seq(from = 1, to = 4, by = 1),.decay = seq(from = 0.1, to = 0.5, by = 0.1))),
          lars=caretModelSpec(method="lars", tuneGrid=expand.grid(.fraction=seq(.01,.99,length=40))),
          rpart=caretModelSpec(method="rpart", tuneGrid=expand.grid(cp = seq(0, .02, .0001)))))
      
      ##-------
      ## Get performance
      model_pred<-list()
      for(mp in 1:length(model_list)){ ## Get prediction in outer test (for algorithm selection)
        p <- predict(model_list[[mp]], newdata=testing)
        
        pred<-data.frame("pred"=as.numeric(p),"obs"=as.numeric(testing$group))
        model_pred[[mp]]<-pred
        
      }
      names(model_pred)<-names(model_list)
      
      res<-list(model_list,model_pred)
      names(res)<-c("models","preds")
      outerML[[i]]<-res

    } 
    
  } # outer
  
  ##----------- (Step2)
  ## Select best algorithm
  if(var.type=="cat"){
    STATS<-matrix(data=0,ncol=length(outerML[[1]]$models),nrow=13)
    colnames(STATS)<-names(outerML[[1]]$models)
    
    for(mod in 1:ncol(STATS)){
      
      tmp<-outerML[[1]]$preds[[mod]]
      cm.model<-outerML[[1]]$cm[[mod]]
      stats<-c(cm.model$overall["Accuracy"],
               cm.model$byClass[c("Balanced Accuracy","Precision","Recall","F1",
                                  "Prevalence","Detection Rate","Detection Prevalence",
                                  "Sensitivity","Specificity","Pos Pred Value","Neg Pred Value")])
      
      count<-ifelse(is.na(stats)==T,0,1)
      
      for(i in 2:length(outerML)){
        tmp<-rbind(tmp,outerML[[i]]$preds[[mod]])
        
        cm.model<-outerML[[i]]$cm[[mod]]
        stats2<-c(cm.model$overall["Accuracy"],
                  cm.model$byClass[c("Balanced Accuracy","Precision","Recall","F1",
                                     "Prevalence","Detection Rate","Detection Prevalence",
                                     "Sensitivity","Specificity","Pos Pred Value","Neg Pred Value")])
        
        count<-count+ifelse(is.na(stats2)==T,0,1)
        stats2[is.na(stats2)]<-0
        stats<-stats+stats2
      }
      stats<-stats/count
      
      
      x<-NA
      tryCatch({
        x <- evalm(tmp)
        x<-x[[7]][[1]]["AUC-ROC","Score"]
      }, error=function(e){})
      stats<-c(x,stats)
      names(stats)[1]<-"AUC"
      
      STATS[,mod]<-stats
    }
    
    rownames(STATS)<-names(stats)
    
    STATS<-STATS[,order(STATS["AUC",],decreasing=T)]
    
    ##-------
    ## Get best Tunning
    
    bestTune<-outerML[[1]]$models[[colnames(STATS)[1]]]$bestTune
    for(i in 2:length(outerML)){
      bestTune<-rbind(bestTune,outerML[[i]]$models[[colnames(STATS)[1]]]$bestTune)
    }
    
    bestTune<-apply(bestTune,2,mean)
    bestTune<-as.data.frame(t(bestTune),nrow=1)
    colnames(bestTune)<-paste0(".",colnames(bestTune))
    
    ##-------
    ## Get final model (with all samples): data
    
    fit.model<-train(group ~ .,data=data,method=colnames(STATS)[1],
                     tuneGrid=bestTune)
    
    stats<-STATS[,1]
    
    preds<-NULL
    for(pr in 1:length(outerML)){
      preds<-rbind(preds,outerML[[pr]]$preds[[colnames(STATS)[1]]])
    }
    
    auc <- evalm(preds)
    
    ##-------
    ## Return model
    
    model<-list(fit.model,stats,auc)
    names(model)<-c("model","stats","auc")
    
    return(model)
    
  }else{ #num
    
    STATS<-rep(0,length(outerML[[1]]$models))
    names(STATS)<-names(outerML[[1]]$models)
    
    for(mod in 1:length(outerML[[1]]$models)){
      
      tmp<-outerML[[1]]$preds[[mod]]
      for(pr in 2:length(outerML)){
        tmp<-rbind(tmp,outerML[[pr]]$preds[[mod]])
      }
      STATS[mod]<-cor(tmp[,1],tmp[,2])
    }
    
    STATS<-STATS[order(STATS,decreasing = T)]
    
    
    ##-------
    ## Get best Tunning
    
    bestTune<-outerML[[1]]$models[[names(STATS)[1]]]$bestTune
    for(i in 2:length(outerML)){
      bestTune<-rbind(bestTune,outerML[[i]]$models[[names(STATS)[1]]]$bestTune)
    }
    
    bestTune<-apply(bestTune,2,mean)
    bestTune<-as.data.frame(t(bestTune),nrow=1)
    colnames(bestTune)<-paste0(".",colnames(bestTune))
    
    
    fit.model<-train(group ~ .,data=data,method=names(STATS)[1],
                     tuneGrid=bestTune)
    
    stats<-STATS[1]
    
    ##-------
    ## Return model
    
    model<-list(fit.model,stats)
    names(model)<-c("model","stats")
    
    return(model)
    
    #
  }
  
}

