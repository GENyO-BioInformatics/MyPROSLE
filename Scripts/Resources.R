##############################
## MyPROSLE 
## R version R version 4.0.4
## Toro-Dominguez, Daniel (daniel.toro@genyo.es)
##############################
## Functions

## ................................................................

#-------
## Chech installed packages
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)>0){
    cat(paste0("\nInstalling ", paste(new.pkg,collapse = ", ")))
    if (!requireNamespace("BiocManager", quietly = TRUE)){
      install.packages("BiocManager")
    }
    bioconductor_packages <- BiocManager::available()
    bioconductor_packages <- bioconductor_packages[bioconductor_packages %in% new.pkg]
    if (length(bioconductor_packages) > 0){
      BiocManager::install(bioconductor_packages, dependencies = TRUE,ask=FALSE)
    }
    new.pkg = new.pkg[!(new.pkg %in% bioconductor_packages)]
    if (length(new.pkg)>0){
      install.packages(new.pkg, dependencies = TRUE)
    }
  }
  res <- lapply(pkg,load.packages)
} ## Character vector contains package names

#-------
## Load packages
load.packages <- function(pkg){
  cat(paste0("\nLoading ",pkg))
  suppressMessages(require(pkg,character.only = T))
} ## Character vector contains package names

#-------
## Log Normalization
norm.log<-function(data){
  
  qx <- as.numeric(quantile(data, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) { data[which(data <= 0)] <- NaN
  data <- log2(data) }
  return(data)
} ## Gene expression matrix

#-------
## Summary of variance of genes
summ.var<-function(data){ ## Gene expression matrix
  
  require(caret)
  data<-t(data)
  nzv <- nearZeroVar(data, saveMetrics= TRUE)
  return(nzv)
}

#-------
## Retrieve datasets size info
data.info<-function(data.list){ 
  
  require("crayon")
  
  totalSLE<-0
  totalH<-0
  for(i in 1:length(data.list)){
   sle<-dim(data.list[[i]]$SLE)[2]
   h<-dim(data.list[[i]]$H)[2]
   cat(paste0("\n",names(DATA[i]),": nºSLE (",sep="",sle,"), nºH (",h,")")) 
   totalSLE<-totalSLE+sle
   totalH<-totalH+h
  }
  cat(crayon::bold(paste0("\nTotal: nºSLE (",sep="",totalSLE,"), nºH (",totalH,")")) )
}

#-------
## Function to add variables (options) to a list 
AddOpt<-function(option.list,name.option,variable.value){
  option.list[[name.option]]<-variable.value
  return(option.list)
}

#-------
## Calculate median of expression by gene
calculate_medians = function(gene, ## Gene to annotate
                             genome, ## Table with gene-probe set conections
                             expressionMatrix){ ## Gene-expression data
  
  probe_ids = genome[genome$toGenes==gene,"fromGenes"] ## Select probes for each gene
  if (length(probe_ids)>1){
    res =  colMedians(as.matrix(expressionMatrix[probe_ids,])) ## Median of all probes 
  } else{
    res <- as.numeric(expressionMatrix[probe_ids,])  
  }
  return(res)
}

#-------
## Function to annotate genes from a gene-expression matrix
annotateGenes<-function(data,
                        toGenes='external_gene_name',
                        fromGenes='ensembl_gene_id'){
  
  require("biomaRt")
  require("parallel")
  
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #,host="uswest.ensembl.org",ensemblRedirect = FALSE)
  genome = getBM(attributes=c(toGenes,fromGenes),mart = ensembl)
  colnames(genome)<-c("toGenes","fromGenes")
  genome <- genome[genome$fromGenes != "",]
  genome = genome[genome$fromGenes %in% rownames(data),]
  data = data[genome$fromGenes,]
  finalGenes = unique(genome$toGenes)
  
  if(as.character(Sys.info()["sysname"])=="Windows"){
    nCores = 1
  }else{
    nCores = detectCores()
  }
  
  temp = mclapply(finalGenes,calculate_medians,genome = genome,expressionMatrix = data,mc.cores = nCores)
  temp = as.data.frame(do.call("rbind", temp))
  rownames(temp) = finalGenes
  colnames(temp) = colnames(data)
  
  return(temp)
}

#-------
## Function to retrieve Mscore for one patient and one path against healthy control distribution
Get.path.mscore<-function(path,            ## Genes from one path (module)
                          Patient,         ## Numeric vector of gene expression values of a patient (with Gene_Symbol as names)
                          Healthy,         ## Gene expression matrix of healthy samples
                          method = "mean") ## proportion (Nº of high dysregulated genes/ Nº genes in the paths), mean (mean os zscores for genes of the path)
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
Get.MyPROSLE<-function(modules.list,      ## List of paths (modules) with their genes
                      Patient,            ## Numeric vector of gene expression values of a patient (with Gene_Symbol as names) or dataframe with multiple patients
                      Healthy,            ## Gene expression matrix of healthy samples
                      method,             ## proportion (Nº of high dysregulated genes/ Nº genes in the paths), mean (mean os zscores for genes of the path)
                      pathReference=NULL, ## Only for Healthy=NULL, path for Reference.RData    
                      nk=5)                ## Only for Healthy=NULL, Number of most similar samples to impute mscore
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


#-------
##Get Distance matrix from Mscore Matrix
#GetDist<-function(list.data, ## list with Mscore matrix
#                  int){ ## position of Mscore matrix in the list
#  require("SNFtool")
#  data<-list.data[[1]][[int]]
# res<-(dist2(as.matrix(data),as.matrix(data)))^(1/2)
#  return(res)
#}


#-------
## Get cluster stability
ClusterStats<-function(list.i, ## permutation list
                       Dist.list, ## List of distance matrices
                       combinations){ ## Matrix with parameter combinations
  require("SNFtool")
  require("NbClust")
  
  nDatasets<-as.numeric(unlist(list.i))
  Times = 30; 	# Number of Iterations, usually (10~20)
  
  ListaNB<-list()
  #pb = txtProgressBar(min = 0, max = nrow(combinations), initial = 0)
  for(i in 1:nrow(combinations)){ # combinations
    #setTxtProgressBar(pb,i)
    ## Parameter selection
    K<<-combinations[i,1]
    Alpha<<-combinations[i,2]
    
    ## Dataset selection
    sel.samples<-sample(1:length(Dist.list),nDatasets,replace=F)
    
    tmp.list<-list()
    for(l.i in 1:length(Dist.list)){
      data<-Dist.list[[l.i]]
      w<-affinityMatrix(data,K,Alpha)
      tmp.list[[l.i]]<-w
    }
    tmp.list<-tmp.list[sel.samples]
    
    #tmp.list<-Dist.list[sel.samples]
    #tmp.list<-lapply(tmp.list,function(x){
    #res<-affinityMatrix(as.matrix(x), K, Alpha)
    #})
    
    W = SNF(tmp.list, K, Times)
    Nb<-NbClust(data = W, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 15,
                method = "complete", index = "ch", alphaBeale = 0.1)
    ListaNB[[i]]<-Nb$All.index
    
  } # combinations
  
  ## Get summary
  nbclusts<-lapply(ListaNB,function(x){
    res<-as.numeric(names(x[order(x,decreasing=T)][1]))
  })
  nbclusts<-unlist(nbclusts)
  res<-table(nbclusts)
  
  res<-c(as.numeric(res["2"]),as.numeric(res["3"]),as.numeric(res["4"]),as.numeric(res["5"]),
         as.numeric(res["6"]),as.numeric(res["7"]),as.numeric(res["8"]),as.numeric(res["9"]),
         as.numeric(res["10"]),as.numeric(res["11"]),as.numeric(res["12"]),as.numeric(res["13"]),
         as.numeric(res["14"]),as.numeric(res["15"]))
  res[is.na(res)]<-0
  
  res<-data.frame(cbind(rep(nDatasets,14),2:15,res))
  colnames(res)<-c("datasets","clusters","nbindex")
  
  return(res)
  
}

#-------
## Get Outliers
Get.Outliers<-function(m,thresh){
  require("raster")
  m1<-mean(as.numeric(m[,1]),na.rm=T)
  m2<-mean(as.numeric(m[,2]),na.rm=T)
  
  res<-NULL
  for(i in 1:nrow(m)){
    distance<-pointDistance(p1=m[i,],p2=c(m1,m2),lonlat = F)
    
    if(distance>=thresh){
      res<-c(res,F)
    }else{
      res<-c(res,T)
    }
  }
  return(res)
}

#-------
## Get gaps for heatmap
GetGap<-function(vect){ ## vect is a order vector
  gaps<-NULL
  id<-vect[1]
  for(i in 1:length(vect)){
    if(vect[i]!=id){ ## new id
      gaps<-c(gaps,i)
      id<-vect[i] ##
    }  
  }
  gaps<-gaps-1
  return(gaps)
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
##-------
## Function to get significances for time events (survival analysis)
SurvAnalysis<-function(rem.matrix,     ## REM matrix (REM1, REM2)
                       mscores.matrix,   ## Mscores matrix (Mscores, Mscores.m)
                       thresh,           ## Threshold for Mscore (to define classes)
                       pval.th=0.05,     ## Significance threshold
                       ann=NULL,         ## Annotation matrix for individual analysis
                       dir = "Positive", ## "Positive", "Negative" or "Absolute" for considering high dysregulation
                       Type="Last"){     ## "Last" or "Next"
  
  require("survival")
  require("dplyr")
  require("survminer")
  
  ## Select mscores for suitable patients (contained in rem.matrix)
  M<-as.data.frame(t(mscores.matrix[,rem.matrix$Sample]))
  
  ## Select time-variable
  switch(Type, 
         Last={
           cat("\nPerforming analysis for time from LAST active state...")
           M<-cbind(as.numeric(rem.matrix$TfLS),M); colnames(M)[1]<-"Time"
         },
         Next={
           cat("\nPerforming analysis for time from NEXT active state...")
           M<-cbind(as.numeric(rem.matrix$TfNS),M); colnames(M)[1]<-"Time"
         })
  
  ## Only for individual analysis (for the scatter plot)
  if(!is.null(ann)){
    ann[colnames(M),]
    Pvalues<-NULL
  }else{
    PlotList<-list()
  }
  
  nsig<-1
  for(i in 2:ncol(M)){
    
    name<-colnames(M)[i]
    data<-M[,c(1,i)]; colnames(data)[2]<-c("Event")
    
    ## Select dysregulation direction 
    switch(dir,
           Positive={
             data[,2]<-ifelse(data[,2]>=thresh,"High","Low" )
           },
           Negative={
             data[,2]<-ifelse(data[,2]<=(thresh),"High","Low" )
           },
           Absolute={
             data[,2]<-ifelse(abs(data[,2])>=thresh,"High","Low" )
           })
    
    ## Get statistical significance
    Censored<-c(rep(1,nrow(data))); data<-cbind(data,Censored)
    survObject <<- Surv(data$Time, data$Censored)
    
    if(length(unique(data$Event))>1){
      
      fit <<- survfit(survObject ~ Event, data = data)
      pval<-surv_pvalue(fit,data); pval<-pval$pval
      
      if(is.null(ann)){
        
        if(pval<=pval.th){
          p<<-ggsurvplot(fit, data = data, pval = TRUE,pval.size=3,title=name,ylab=NULL,size=0.7,
                         font.title=c(10,"bold","black"),pval.coord=c(250,0.9),
                         font.tickslab = c(8, "plain", "black"))
          #print(p)
          PlotList[[nsig]]<-p
          nsig<-nsig+1
        }
        
      }else{
        Pvalues<-c(Pvalues,pval)
      }
      
      
      
    }else{
      if(!is.null(ann)){
        Pvalues<-c(Pvalues,NA)
      }
    }
    
   
  }
  
  
  if(!is.null(ann)){
    cat("\nReturning table with pvalues")
    ann<-cbind(Pvalues,ann)
    return(ann)
  }else{
    cat("\nReturning list with information for significant Survival plots")
    return(PlotList) 
  }
  
}

##-------
## Function to summarize module mscores in the SLE-signatures
GetSignatures<-function(listM,   ## List of datasets with Mscores or matrix
                        ann){     ## Annotation matrix
  
  Signatures<-unique(ann$category)
  
  if(class(listM)=="list"){
    res<-list()
    for(i in 1:length(listM)){
      
      tmpSLE<-listM[[i]]$SLE
      SignMatrix<-data.frame(matrix(ncol=ncol(tmpSLE),nrow=0))
      
      for(sign in 1:length(Signatures)){
        
        features<-rownames(ann[ifelse(ann$category==Signatures[sign],T,F),])
        tmp<-tmpSLE[features,]
        tmp<-apply(tmp,2,mean)
        
        SignMatrix<-rbind(SignMatrix,tmp)
      }
      rownames(SignMatrix)<-Signatures
      colnames(SignMatrix)<-colnames(listM[[i]]$SLE)
      res[[i]]<-SignMatrix
      names(res)[i]<-paste0("dataset",sep="",i)
    }
    
  }else{
      tmpSLE<-listM
      SignMatrix<-data.frame(matrix(ncol=ncol(tmpSLE),nrow=0))
      
      for(sign in 1:length(Signatures)){
        
        features<-rownames(ann[ifelse(ann$category==Signatures[sign],T,F),])
        tmp<-tmpSLE[features,]
        tmp<-apply(tmp,2,mean)
        
        SignMatrix<-rbind(SignMatrix,tmp)
      }
      rownames(SignMatrix)<-Signatures
      colnames(SignMatrix)<-colnames(tmpSLE)
      res<-SignMatrix
  }

  return(res)
}

##-------
## Function to compare values between two groups
GroupComparison<-function(mdata,             ## data.frame with features in rows and samples in columns
                          group,             ## vector indicating sample classes
                          annot=NULL,        ## (Optional) vector of annotation. names must coincide with rownames from mdata
                          namePlot,          ## Name for the plot
                          Width,             ## Width in inches
                          Height,            ## Height in inches
                          pval.th=0.01,      ## Pvalue threshold to consider significance
                          method="wilcox",   ## Method for the statistical comparison
                          Ylim=c(-5,5),      ## Limits for Y value in plot
                          textSize=6,        ## font size for pvalue
                          colori=c("cadetblue2","darkseagreen4"), ## Colors for each group
                          groupName=c("V","L"), ## Name for each group
                          pointPval=3){      ## Height to put pval label
  
  ## Get significant comparisons
  sign<-NULL
  nms<-NULL
  val<-NULL
  for(i in 1:nrow(mdata)){ ## Test for each feature
    
    Values<-as.numeric(mdata[i,])
    switch(method,
           wilcox={
             
             res<-min(wilcox.test(Values~group,alternative="greater")$p.value,
                      wilcox.test(Values~group,alternative="less")$p.value)
             
             if(res<=pval.th){
               sign<-c(sign,i)
               nms<-c(nms,rownames(mdata)[i])
               val<-c(val,res)
             }
           },
           other={
             ## Implement other test (if neccesary)
           })
  }
  
  ## Get values for boxplots
  meltData<-data.frame(matrix(ncol=3,nrow=0))
  colnames(meltData)<-c("Group","Mscore","Colors")
  
  for(i in 1:length(sign)){
    if(is.null(annot)){
      nme<-paste0("Cl_",sep="",nms[i])
    }else{
      x<-annot[nms[i]]
      nme<-paste0("Cl_",x,"_",sep="",nms[i])
    }
    
    Values<-mdata[nms[i],]
    
    tmp<-as.data.frame(cbind(Group,as.numeric(Values)))
    tmp[,1]<-ifelse(tmp[,1]==1,groupName[1],groupName[2]) 
    
    colors<-ifelse(tmp[,1]==groupName[1],colori[1],colori[2])
    tmp<-cbind(tmp,colors)
    
    tmp$Group<-paste0(nme,sep="","_",tmp$Group)
    colnames(tmp)<-c("Group","Mscore","Colors")
    meltData<-rbind(meltData,tmp)
    
  }
  
  ## Boxplots
  meltData$Group<-factor(meltData$Group,levels = unique(meltData$Group))
  
  ## points
  points<-seq(from=1.5,to=(2*length(sign))+1,by=2)
  
  p1<-ggplot(meltData,aes(x=Group,y=Mscore,fill=Colors))+
    geom_hline(yintercept = 0, color="grey50")+
    theme_bw() +theme(axis.text.x = element_text(angle = 90))+
    geom_jitter(size=2.5,color=meltData$Colors,alpha=0.8,shape=16)+
    geom_boxplot(alpha=0.5,outlier.alpha=0) + ylim(Ylim[1],Ylim[2]) + 
    scale_fill_manual(values=c("cadetblue2","darkseagreen4"))
  
  for(i in 1:length(points)){
    p1<-p1 + geom_text(x=points[i], 
                       y=pointPval, 
                       label=as.character(round(val[i],digits = 4)),
                       size=textSize,
                       angle = 0,color="black", check_overlap = TRUE)
    #p1<-p1 + geom_segment(aes(x=points[i]-0.5,y=pointPval-0.3,
    #                          xend=points[i]+0.5,yend=pointPval-0.3))
  }
  
  tiff(filename=paste0(opt$resultsPath,sep="","/",namePlot,".tiff"),res = 300,width = Width, height = Height,units="in")
  plot(p1)
  dev.off()
  
}

##-------
## Split data in a balanced way between multiple classes
splitData<-function(class.vector,  ## Vector with class labels
                    p){            ## %, size of the subsets (in 0-1 scale)
  require("caret")
  res.vect<-as.numeric(createDataPartition(class.vector, p = p, list=F,times=1))
  return(res.vect)
}

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

##-------
## Reduce colinear features for regression models
vif_func<-function(in_frame,thresh=10,trace=T,...){
  
  library(fmsb)
  
  if(any(!'data.frame' %in% class(in_frame))) in_frame<-data.frame(in_frame)
  
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  var_names <- names(in_frame)
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in <- formula(paste(val, '~', form))
    vif_init<-rbind(vif_init, c(val, VIF(lm(form_in, data = in_frame, ...))))
  }
  vif_max<-max(as.numeric(vif_init[,2]), na.rm = TRUE)
  
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(var_names)
  }
  else{
    
    in_dat<-in_frame
    
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      
      vif_vals<-NULL
      var_names <- names(in_dat)
      
      for(val in var_names){
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = '+')
        form_in <- formula(paste(val, '~', form))
        vif_add<-VIF(lm(form_in, data = in_dat, ...))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2]), na.rm = TRUE))[1]
      
      vif_max<-as.numeric(vif_vals[max_row,2])
      
      if(vif_max<thresh) break
      
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        flush.console()
      }
      
      in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
      
    }
    
    return(names(in_dat))
    
  }
  
}


##

## Linear model Predictor
Getlm<-function(data,           ## Data with class to predict in the first columns (named as Group), features in columns and samples in rows
                kfold=10,       ## Either the number of folds or number of resampling iterations
                repeatedCv=30,  ## For repeated k-fold cross-validation only: the number of complete sets of folds to compute
                bootstrap=0.8,  ## Percent to divide training and test (0-1 scale)
                perm.bias=10,   ## Number of permutations in train and test selection
                max.features=15, ## Maximun number of features to test (one by one). NULL for all
                seed=12345){ 
  
  require("caret")
  require("stringr")
  require("parallel")
  require("car")
  
  set.seed(seed)
  
  fitControl <- trainControl(method = "repeatedcv",
                             number = kfold,
                             repeats = repeatedCv,
                             allowParallel = T,
                             verboseIter = T,
                             savePredictions = TRUE,
                             classProbs = F)

  ## Allow parallel
  clusters = makeCluster(detectCores()-1)
  registerDoParallel(clusters)
  
  ## Fix colnames from data
  colnames(data)<-FixSymbol(vect=colnames(data),symbol=c("/"," "),toSymbol=c(".",""))
  
  ## Filter variables by vif
  cat("\nFiltering features by VIF...")
  rem.vars<-vif_func(in_frame = data[,-1])
  
  if(length(rem.vars)==ncol(data)-1){
    
    cat("\nAny available feature")
    results<-NULL
  }else{
    data<-data[,c("values",rem.vars)]
    
    ##------
    ## Get Feature importance
    model <- train(values~., data = data,method="lm",trControl = fitControl)
    imp<-varImp(model)
    features.imp<-as.data.frame(cbind(rownames(imp$importance),imp$importance))
    colnames(features.imp)[1]<-"Features"
    features.imp$Overall<-as.numeric(features.imp$Overall)
    features.imp<-features.imp[order(features.imp$Overall,decreasing=T),]
    features.imp$Features<-factor(x = features.imp$Features,levels = unique(features.imp$Features))
    
    importance<-as.character(features.imp$Features[features.imp$Overall!=0])
    
    featureSets<-1:length(importance)
    if(!is.null(max.features)){
      if(max.features<length(featureSets)){
        featureSets<-featureSets[1:max.features]
      }
    }
    
    
    ##------
    ## Getting model with minimal features
    cat("\nGetting model with minimal features...\n")
    
    perm.bias.vals<-list()
    for(b in 1:perm.bias){
      perm.bias.vals[[b]]<-sample(1:nrow(data),size = round(nrow(data)*bootstrap,digits = 0),replace = F)
      
    }
    
    resTable.lm<-data.frame(matrix(data=NA,nrow=length(perm.bias.vals),ncol=length(featureSets)))
    colnames(resTable.lm)<-featureSets
    rownames(resTable.lm)<-1:length(perm.bias.vals)
    
    for(f in 1:length(featureSets)){ ## feature number
      tmp<-data[,c("values",as.character(importance[featureSets[1:f]]))]
      for(b in 1:length(perm.bias.vals)){ ## selection bias
        TRAIN<-tmp[as.numeric(perm.bias.vals[b][[1]]),]
        TEST<-tmp[-as.numeric(perm.bias.vals[b][[1]]),]
        
        model.lm <- train(values~., data = TRAIN,method="lm",trControl = fitControl)
        res.lm<-predict(model.lm, newdata = TEST)
        resTable.lm[b,f]<-cor(TEST$values,res.lm)
      }
    }
    
    ##------
    ## Getting final model
    cat("\nGetting final model\n")
    
    nFeatures<-apply(resTable.lm,2,mean)
    nFeatures<-as.numeric(names(nFeatures[order(nFeatures,decreasing = T)])[1])
    
    subset<-as.numeric(rownames(resTable.lm)[ order(resTable.lm[,as.character(nFeatures)],decreasing=T)[1]])
    #subset<-sample(1:nrow(data),size = round(nrow(data)*bootstrap,digits = 0),replace = F)
    
    tmp<-data[,c("values",as.character(importance[1:nFeatures]))]
    features.selected<-as.character(importance[1:nFeatures])
    
    TRAIN<-tmp[as.numeric(perm.bias.vals[subset][[1]]),]
    TEST<-tmp[-as.numeric(perm.bias.vals[subset][[1]]),]
    
    model.lm<- train(values ~., data = TRAIN,trControl = fitControl,method="lm")
    model.glm<- train(values ~., data = TRAIN,trControl = fitControl,method="glm")
    
    if(model.lm$results$Rsquared>=model.glm$results$Rsquared){
      res.lm<-predict(model.lm, newdata = TEST)
      cm.model<-cor(TEST$values,res.lm)
      fit.model<-model.lm
      #cm.model<-model.lm$results$Rsquared
    }else{
      res.glm<-predict(model.glm, newdata = TEST)
      cm.model<-c(cor(TEST$values,res.glm),as.numeric(apply(resTable.lm,2,mean)[nFeatures]))
      names(cm.model)<-c("max.cor.Test","mean.Rs")
      fit.model<-model.glm
      #cm.model<-model.glm$results$Rsquared
    }
    #cm.model<-cor(TEST$values,res)
    
    results<-list(fit.model,features.selected,cm.model)
    names(results)<-c("model","features","stats")
    
    stopCluster(clusters)
    registerDoSEQ()
    
    return(results)
  }

}


##-------
## Classification model for Worsening
GetModel<-function(data,           ## Data with class to predict in the first columns (named as Group), features in columns and samples in rows
                   method="rf",    ## Method to be applied in classification
                   bootstrap=0.8,  ## Percent to divide training and test (0-1 scale)
                   kfold=10,       ## Either the number of folds or number of resampling iterations
                   repeatedCv=20,  ## For repeated k-fold cross-validation only: the number of complete sets of folds to compute
                   perm.bias=10,   ## Number of permutations in train and test selection
                   seed=123,
                   accuracy.type="Balanced Accuracy", ## "Balanced Accuracy" or "Accuracy to generate the better model
                   featureSets=c(1,2,3,4,5,6,7,8,9,10,15,20,25)){ ## Subsets of features to calculate minimun number of features

  require("caret")
  require("stringr")
  
  ## Set parameters to fit the model
  set.seed(seed)
  
  fitControl <- trainControl(method = "repeatedcv",
                             number = kfold,
                             repeats = repeatedCv,
                             allowParallel = T,
                             verboseIter = T,
                             savePredictions = TRUE,
                             classProbs = T)
  
  ## Allow parallel
  clusters = makeCluster(detectCores()-1)
  registerDoParallel(clusters)
  
  ## Fix colnames from data
  colnames(data)<-FixSymbol(vect=colnames(data),symbol=c("/"," "),toSymbol=c(".",""))
  
  ##------  
  ## Step 1: Get best parameters
  
  ## split TRAIN and TEST
  sel<-splitData(data$Group,p = bootstrap)
  TRAIN<-data[sel,]; TEST<-data[-sel,]
  
  ## Tunning parameters
  cat("\nTuning parameters for Random Forest algorithm...\n")
  # Number of features to consider at every split (mtry)
  n<-10
  if(n>round(sqrt(ncol(data))+2)){
    n<-round(sqrt(ncol(data))+2)
  }
  
  mtry<-unique(seq(from=2,to=round(sqrt(ncol(data))+2,digits=0),by=round((round(sqrt(ncol(data))+2))/n,digits=0)))
  # Number of trees in random forest (ntree)
  ntree<-seq(from=200,to=2000,by=100)
  #ntree<-as.list(ntree)
  tunegrid<-expand.grid(mtry=mtry)
           
  res<-list()
  pb = txtProgressBar(min = 0, max = length(ntree), initial = 0)
  for(i in 1:length(ntree)){
     model<- invisible(train(Group ~., data = TRAIN,method="rf",
                   trControl = fitControl,tuneGrid=tunegrid,ntree=ntree[i]))
     res[[i]]<-model
     setTxtProgressBar(pb,i)
  }
           
  ## Search best parameters (high accuracy)
  accuracy<-0
  trees<-NULL
  importance<-NULL
  for(i in 1:length(res)){
             
    tmp<-res[[i]]$results
    tmp<-tmp[order(tmp$Accuracy,decreasing=T),]
    if(tmp$Accuracy[1]>accuracy){
       mtry<-tmp$mtry[1]
       trees<-ntree[i]
       accuracy<-tmp$Accuracy[1]
               
       imp<-varImp(res[[i]])
       features.imp<-as.data.frame(cbind(rownames(imp$importance),imp$importance))
       colnames(features.imp)[1]<-"Features"
       features.imp$Overall<-as.numeric(features.imp$Overall)
       features.imp<-features.imp[order(features.imp$Overall,decreasing=T),]
       features.imp$Features<-factor(x = features.imp$Features,levels = unique(features.imp$Features))
               
       importance<-as.character(features.imp$Features)
     }
  }
           
  ## Best mtry and ntree selected and feature importance assigned

  cat("\nOptimal parameters computed\n")
  
  ##------  
  ## Step 2: Test number of features and TRAN-TEST selection bias
  
  cat("\nGetting model with minimal features...\n")
  
  splitN<-function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
  sel<-as.numeric(createDataPartition(data$Group, p = bootstrap, list=F,times=perm.bias))
  perm.bias.vals<-splitN(sel,perm.bias)
           
  resTable<-data.frame(matrix(data=NA,nrow=length(perm.bias.vals),ncol=length(featureSets)))
  colnames(resTable)<-featureSets
  rownames(resTable)<-1:length(perm.bias.vals)
           
           
  for(i in 1:length(featureSets)){ ## feature number
             
     tmp<-data[,c("Group",as.character(importance[featureSets[1:i]]))]
             
     for(b in 1:length(perm.bias.vals)){ ## selection bias
               
         TRAIN<-tmp[as.numeric(perm.bias.vals[b][[1]]),]
         TEST<-tmp[-as.numeric(perm.bias.vals[b][[1]]),]
               
         model<- invisible(train(Group ~., data = TRAIN,method="rf",
                       trControl = fitControl,tuneGrid=expand.grid(mtry=mtry),ntree=trees))
               
         res<-predict(model, newdata = TEST,type = "prob")
         
         pred<-as.factor(ifelse(res[,1]>=0.5, colnames(res)[1], colnames(res)[2]))
         real<-as.factor(TEST$Group)
         cm.model<-confusionMatrix(pred,real)              
         if(accuracy.type=="Accuracy"){
           acc<-as.numeric(cm.model$overall[accuracy.type])
         }else if(accuracy.type=="Balanced Accuracy"){
           acc<-as.numeric(cm.model$byClass[accuracy.type])
         }
         resTable[b,i]<-acc
         }
     }
           
  ##------  
  ## Step 3: Get final model
  cat("Extracting Final Model...\n")
           
  nFeatures<-apply(resTable,2,mean)
  nFeatures<-as.numeric(names(nFeatures[order(nFeatures,decreasing = T)])[1])
           
  subset<-as.numeric(rownames(resTable)[ order(resTable[,as.character(nFeatures)],decreasing=T)[1]])
           
  tmp<-data[,c("Group",as.character(importance[1:nFeatures]))]
  features.selected<-as.character(importance[1:nFeatures])
           
  TRAIN<-tmp[as.numeric(perm.bias.vals[subset][[1]]),]
  TEST<-tmp[-as.numeric(perm.bias.vals[subset][[1]]),]
           
  fit.model<- train(Group ~., data = TRAIN,method="rf",
                    trControl = fitControl,tuneGrid=expand.grid(mtry=mtry),ntree=trees)
           
  res<-predict(fit.model, newdata = TEST,type = "prob")
  cm.model<-confusionMatrix(as.factor(ifelse(res[,1]>=0.5, colnames(res)[1], colnames(res)[2])),
                                     as.factor(TEST$Group))   
  results<-list(fit.model,features.selected,cm.model)
  names(results)<-c("model","features","stats")
  
  stopCluster(clusters)
  registerDoSEQ()
  
  return(results)
}


##-------
## Anova test for categorical (yes/no) variables with respect to Mscores
ClinicalAnalyser<-function(variables, ## Name of variables (synonymous)
                           M,         ## M values for all patients
                           data,      ## List of Clinical matrices
                           adjust="fdr"){     
  
  data.toPredict<-NULL
  
  variable.values<-NULL
  for(i in 1:length(data)){
    for(v in 1:length(variables)){
      if(variables[v] %in% colnames(data[i][[1]])){
        x<-data[i][[1]][,variables[i]]
        names(x)<-rownames(data[i][[1]])
        
        variable.values<-c(variable.values,x)
        
      }}}
  variable.values<-variable.values[!duplicated(names(variable.values))]
  
  variable.values<-variable.values[ifelse(is.na(variable.values)==T,F,T)]
  modules<-rep(1,nrow(M))
  vectFC<-rep(1,nrow(M))
  
  if(length(table(variable.values))>1){ ## Two categories
    if(table(variable.values)[1]>=5 & table(variable.values)[2]>=5){ ## At least 5 points for category
      
      Mvalues<-M[,names(variable.values)]
      variable.values<-as.character(variable.values)
      
      data.toPredict<-as.data.frame(cbind(variable.values,as.data.frame(t(Mvalues))))
      colnames(data.toPredict)[1]<-"Group"
      
      modules<-rep(1,nrow(Mvalues))
      for(mod in 1:nrow(Mvalues)){ ## Run each module
        
        ## Anova test 
        mtemp<-as.numeric(Mvalues[mod,])
        mtemp<-as.data.frame(cbind(mtemp,variable.values))
        colnames(mtemp)<-c("value","Group")
        
        mtemp$value<-as.numeric(as.character(mtemp$value))
        mtemp$Group<-paste0("Group_",sep="",as.character(mtemp$Group))
        
        aov <- anova(lm(value ~ Group, data = mtemp))
        
        #if(as.numeric(aov$Pr[1])<=0.05){
        modules[mod]<-as.numeric(aov$Pr[1])
        
        sel<-ifelse(mtemp$Group=="Group_Yes",T,F)
        
        
        if(mean(mtemp$value[sel])<mean(mtemp$value[!sel])){
          vectFC[mod]<- -1
        }
        
        #}
      }## Para cada modulo
      
      names(modules)<-rownames(Mvalues)
      modules<-p.adjust(modules,method = adjust)
      
      
      #gplot<-ggplot(mtemp,aes(x=group,y=value))+geom_boxplot()+theme_bw()
      #print(gplot)
      
    }}
  modules<-(-log10(modules))
  modules<-modules*vectFC
  
  res<-list(modules,data.toPredict); names(res)<-c("modules","melt")
  
  return(res)
}


##-------
## Add elements to model list
AddModel<-function(newModel,models,name){
  
  x<-(length(models)+1)
  models[[x]]<-newModel; names(models)[x]<-name
  return(models)
}

##-------
## Create breaks and custom color ranges
Custom.colorRange<-function(hits,colorseq,size=1000){
  
  breaks=seq(min(hits),max(hits),length.out = size)
  cols<-NULL
  i=1
  while(i<length(hits)){
    sel<-ifelse(breaks>=hits[i] & breaks<hits[i+1],T,F)
    n<-table(sel)["TRUE"]
    if((i+1)==(length(hits))){
      n<-n+1
    }
    tmp<-colorRampPalette(colors=c(colorseq[i],colorseq[i+1]))(n = n)
    cols<-c(cols,tmp)
    i<-i+1
  }
  vals<-list(breaks,cols); names(vals)<-c("breaks","colors")
  return(vals)
}






