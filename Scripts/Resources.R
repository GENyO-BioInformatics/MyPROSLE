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
## Load libraries
load.libraries <- function(set="all"){
  
 switch(set,
        "all"={
          pkgs<-c("tmod","caret","matrixStats","biomaRt","NOISeq",
                  "NbClust","SNFtool","tidyr","parallel","ggplot2",
                  "ConsensusClusterPlus","raster","pheatmap","ggbreak",
                  "reshape","survival","dplyr","survminer","ggrepel",
                  "doParallel","stringr","factoextra")
          check.packages(pkgs)
        },
        "00"={
          pkgs<-c("tmod","caret","matrixStats","biomaRt","NOISeq",
                  "NbClust","SNFtool","tidyr","parallel","ggplot2",
                  "ConsensusClusterPlus","raster","pheatmap","ggbreak",
                  "reshape","survival","dplyr","survminer","ggrepel",
                  "doParallel","stringr")
          check.packages(pkgs)
        },
        "01"={
          pkgs<-c("caret","matrixStats","NbClust","SNFtool","factoextra",
                  "tidyr","parallel","ggplot2","ConsensusClusterPlus",
                  "raster","pheatmap")
          check.packages(pkgs)
        },
        "02"={
          pkgs<-c("ggplot2")
          check.packages(pkgs)
        },
        "03"={
          pkgs<-c("ggplot2")
          check.packages(pkgs)
        },
        "04"={
          pkgs<-c("ggplot2","ggbreak","pheatmap")
          check.packages(pkgs)
        },
        "05"={
          pkgs<-c("matrixStats","raster","pheatmap","reshape",
                  "parallel","ggplot2")
          check.packages(pkgs)
        },
        "06"={
          pkgs<-c("ggplot2","reshape","survival","dplyr",
                  "survminer","ggrepel","doParallel","stringr")
          check.packages(pkgs)
        },
        "07"={
          pkgs<-c("ggplot2","reshape","survival","dplyr",
                  "survminer","ggrepel","doParallel","stringr")
          check.packages(pkgs)
        }
        ) 
  
  
}
  
#-------
## Function to add variables (options) to a list 
AddOpt<-function(option.list,name.option,variable.value){
  option.list[[name.option]]<-variable.value
  return(option.list)
}

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
## Retrieve datasets size info
data.info<-function(data.list){ 
  
  require("crayon")
  
  totalSLE<-0
  totalH<-0
  for(i in 1:length(data.list)){
    sle<-dim(data.list[[i]]$SLE)[2]
    h<-dim(data.list[[i]]$H)[2]
    cat(paste0("\n",names(DATA[i]),": n?SLE (",sep="",sle,"), n?H (",h,")")) 
    totalSLE<-totalSLE+sle
    totalH<-totalH+h
  }
  cat(crayon::bold(paste0("\nTotal: n?SLE (",sep="",totalSLE,"), n?H (",totalH,")")) )
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
Get.MyPROSLE<-function(modules.list,      ## List of paths (modules) with their genes
                       Patient,            ## Numeric vector of gene expression values of a patient (with Gene_Symbol as names) or dataframe with multiple patients
                       Healthy,            ## Gene expression matrix of healthy samples
                       method,             ## proportion (N? of high dysregulated genes/ N? genes in the paths), mean (mean os zscores for genes of the path)
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
## Extract suitable clinical data
MergeClinical<-function(variableIDs,
                        clin,
                        minSize=10,
                        Mscore.list,
                        Msig.list,
                        variable.type="categoric",
                        val2val = c("0"="NO","8"="YES")){
  
  ##------
  ## Get clinical values for all clinical datasets
  variable.values<-NULL
  for(clin.i in 1:length(clin)){
    for(var.i in 1:length(variableIDs)){
      if(variableIDs[var.i] %in% colnames(clin[[clin.i]])){
        x<-clin[[clin.i]][,variableIDs[var.i]]
        names(x)<-rownames(clin[[clin.i]])
        variable.values<-c(variable.values,x)
      }
    }
  }
  variable.values<-variable.values[!duplicated(names(variable.values))]
  variable.values<-variable.values[ifelse(is.na(variable.values)==T,F,T)]
  
  ##------
  ## Transform data values
  if(variable.type!="categoric"){
    nme<-names(variable.values)
    variable.values<-as.numeric(as.character(variable.values))
    names(variable.values)<-nme
  }else{
    for(var.i in 1:length(variable.values)){
      for(v2 in 1:length(val2val)){
        if(variable.values[var.i]==names(val2val)[v2]){
          variable.values[var.i]<-val2val[v2]
        }
      }
    }
  }
  variable.values<-variable.values[ifelse(variable.values!="X",T,F)]
  
  ##------
  ## Select variables with a minimun 
  
  if(variable.type=="categoric"){
    if(length(table(variable.values))>1){ ## Two categories
      if(table(variable.values)[1]>=minSize & table(variable.values)[2]>=minSize){ ## At least 'minSize' points for category
        
        Mscore<-as.data.frame(matrix(nrow=nrow(Mscore.list[[1]])))
        Msig<-as.data.frame(matrix(nrow=nrow(Msig.list[[1]])))
        
        for(dat in 1:length(Mscore.list)){
          if(length(intersect(colnames(Mscore.list[[dat]]),names(variable.values)))>0){
            Mscore<-cbind(Mscore,Mscore.list[[dat]][,intersect(colnames(Mscore.list[[dat]]),names(variable.values))]) 
          }
          if(length(intersect(colnames(Msig.list[[dat]]),names(variable.values)))>0){
            Msig<-cbind(Msig,Msig.list[[dat]][,intersect(colnames(Msig.list[[dat]]),names(variable.values))]) 
          }
        }
        Mscore<-Mscore[,names(variable.values)]
        Msig<-Msig[,names(variable.values)]
        
        Mscore<-as.data.frame(t(Mscore))
        Msig<-as.data.frame(t(Msig))
        
        Mscore<-cbind(variable.values,Mscore); colnames(Mscore)[1]<-"group"
        Msig<-cbind(variable.values,Msig); colnames(Msig)[1]<-"group"
        
        res<-list(variable.type,Mscore,Msig)
        names(res)<-c("variable.type","mscore","msig")
        
        return(res)
        
      }else{
        return(NULL)
      }
    }else{
      return(NULL)
    }
    
  }else{ ## Numeric values
    
    if(length(variable.values)>=minSize){ ## At least 'minSize' points
      
      Mscore<-as.data.frame(matrix(nrow=nrow(Mscore.list[[1]])))
      Msig<-as.data.frame(matrix(nrow=nrow(Msig.list[[1]])))
      
      for(dat in 1:length(Mscore.list)){
        if(length(intersect(colnames(Mscore.list[[dat]]),names(variable.values)))>0){
          Mscore<-cbind(Mscore,Mscore.list[[dat]][,intersect(colnames(Mscore.list[[dat]]),names(variable.values))]) 
        }
        if(length(intersect(colnames(Msig.list[[dat]]),names(variable.values)))>0){
          Msig<-cbind(Msig,Msig.list[[dat]][,intersect(colnames(Msig.list[[dat]]),names(variable.values))]) 
        }
      }
      Mscore<-Mscore[,names(variable.values)]
      Msig<-Msig[,names(variable.values)]
      
      Mscore<-as.data.frame(t(Mscore))
      Msig<-as.data.frame(t(Msig))
      
      Mscore<-cbind(variable.values,Mscore); colnames(Mscore)[1]<-"group"
      Msig<-cbind(variable.values,Msig); colnames(Msig)[1]<-"group"
      
      res<-list(variable.type,Mscore,Msig)
      names(res)<-c("variable.type","mscore","msig")
      
      return(res)
      
    }else{
      return(NULL)
    }
    
    
  }
  
  
  
  
}


##-------
## Get Candidates to Ensemble models
toEnsemble<-function(mcor,minCor=0.25){
  candidates<-colnames(mcor)[1]
  forbiden<-rownames(mcor)[ifelse(mcor[,1]<=minCor & mcor[,1]>=0,F,T)]
  
  for(nc in 1:ncol(mcor)){
    if(!colnames(mcor)[nc] %in% forbiden){ 
      candidates<-c(candidates,colnames(mcor)[nc])
      forbiden<-c(forbiden,rownames(mcor)[ifelse(mcor[,nc]<=minCor & mcor[,1]>=0,F,T)])
      forbiden<-unique(forbiden)
    }
  }
  return(unique(candidates))
}

##-------
## Get Stats for the models
modelStats<-function(model.i,test){
  
  p <- predict(model.i, newdata=test,type="prob")
  
  if(class(p)=="numeric"){
    cm.model<-confusionMatrix(as.factor(ifelse(p<=0.5,"YES","NO")),as.factor(test$group))
    pred<-data.frame("NO"=p,"YES"=1-p,"obs"=test$group)
  }else{ ##data.frame
    cm.model<-confusionMatrix(as.factor(ifelse(p$NO<=0.5,"YES","NO")),as.factor(test$group))
    pred<-data.frame("NO"=p$NO,"YES"=1-p$NO,"obs"=test$group)
    
  }
  pred<-na.omit(pred)
  
  stats<-as.numeric(cm.model$overall["Accuracy"])
  
  x<-NA
  tryCatch({
    x <- evalm(pred)
    x<-x[[7]][[1]]["AUC-ROC","Score"]
  }, error=function(e){})
  stats<-c(stats,x)
  names(stats)<-c("Accuracy","AUC")
  
  stats<-c(stats,cm.model$byClass[c("Balanced Accuracy","Precision","Recall","F1",
                                    "Prevalence","Detection Rate","Detection Prevalence",
                                    "Sensitivity","Specificity","Pos Pred Value","Neg Pred Value")])
  return(stats)
}

# https://topepo.github.io/caret/available-models.html
##-------
## Get Prediction models for categorical variables
getML.cat<-function(training,
                    testing,
                    kfold=10,
                    repeats=10,
                    sel.feature=NULL,
                    set.seed=123456788){
  
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
  
  ##------
  ## Feature selection on training set
  if(!is.null(sel.feature)){
    
    ## DESARROLLAR UN PAR DE METODOS, NO PARA ESTE CASO  
    
  }
  
  ##-------
  ## TrainControl
  ## Se divide el training en 10 partes aleatorias (9 train, 1 test)
  ## Se repite 10 veces para calcular hiperparametros optimos
  
  my_control <- trainControl(
    method="repeatedcv",
    number=kfold,
    savePredictions="final",
    classProbs=TRUE,
    index=createResample(training$group, kfold),
    repeats=repeats)
  
  ##-------
  ## GetModels ## MEJORA A HACER, INCLUIR MAS MODELOS Y SELECCIONAR DESDE LA FUNCION CUALES SE QUIEREN
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
      nb=caretModelSpec(method="nb", tuneGrid=expand.grid(fL=c(0,0.5,1.0), usekernel = TRUE, adjust=c(0.5,1.0))))
    #
  )
  
  
  ##-------
  ## Extract stats
  res<-as.data.frame(matrix(data=0,ncol=length(model_list),nrow=13))
  colnames(res)<-names(model_list)
  for(cm in 1:length(model_list)){
    res[,cm]<-modelStats(model.i = model_list[[cm]],test = testing)
  }
  rownames(res)<-c("Accuracy","AUC","Balanced Accuracy","Precision","Recall","F1",
                   "Prevalence","Detection Rate","Detection Prevalence",
                   "Sensitivity","Specificity","Pos Pred Value","Neg Pred Value")
  
  res<-res[,order(-res["AUC",], -res["Balanced Accuracy",])]
  
  
  ##-------
  ## Aggregate best two models
  cor.bestModels<-modelCor(resamples(model_list))
  cor.bestModels<-cor.bestModels[colnames(res),colnames(res)]
  tags<-toEnsemble(cor.bestModels)
  
  if(length(tags)>1){
    
    listM<-model_list[ifelse(names(model_list) %in% tags,T,F)]
    aggregatedModel <- caretEnsemble(
      listM,
      trControl=trainControl(
        number=10,
        verboseIter = TRUE
      ))
    
    ensemble<-modelStats(model.i = aggregatedModel,test = testing)
    res<-cbind(ensemble,res)
    res<-res[,order(-res["AUC",], -res["Balanced Accuracy",])]
    
    model_list[["ensemble"]]<-aggregatedModel
  }
  
  model_list<-model_list[colnames(res)]
  
  models.results<-list(res,model_list)
  names(models.results)<-c("stats","models")
  
  return(models.results)
  
}


##--------
## Get Prediction models for numerical variables
getML.num<-function(training,
                    testing,
                    kfold=10,
                    repeats=10,
                    sel.feature=NULL,
                    set.seed=123456788){
  
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
  
  ##------
  ## Feature selection on training set
  if(!is.null(sel.feature)){
    
    ## DESARROLLAR UN PAR DE METODOS, NO PARA ESTE CASO  
    
  }
  
  my_control <- trainControl(
    method="repeatedcv",
    number=kfold,
    savePredictions="final",
    classProbs=TRUE,
    index=createResample(training$group, kfold),
    repeats=repeats)
  
  ##-------
  
  
  
  ## GetModels ## MEJORA A HACER, INCLUIR MAS MODELOS Y SELECCIONAR DESDE LA FUNCION CUALES SE QUIEREN
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
      rpart=caretModelSpec(method="rpart", tuneGrid=expand.grid(cp = seq(0, .02, .0001))))
    
  )
  
  ##-------
  ## Extract stats
  res<-as.data.frame(matrix(data=0,ncol=length(model_list),nrow=4))
  colnames(res)<-names(model_list)
  for(cm in 1:length(model_list)){
    x<-model_list[[cm]]$results
    x<-x[order(-x$Rsquared,x$RMSE),c("Rsquared","RMSE","MAE")]
    x<-as.numeric(x[1,])
    x<-c(cor(testing$group,predict(model_list[[cm]],newdata = testing),method = "pearson"),x)
    res[,cm]<-x
  }
  rownames(res)<-c("Correlation","Rsquared","RMSE","MAE")
  
  res<-res[,order(-res["Correlation",], -res["Rsquared",])]
  
  
  ##-------
  ## Aggregate best two models
  cor.bestModels<-modelCor(resamples(model_list))
  cor.bestModels<-cor.bestModels[colnames(res),colnames(res)]
  tags<-toEnsemble(cor.bestModels)
  
  if(length(tags)>1){
    
    listM<-model_list[ifelse(names(model_list) %in% tags,T,F)]
    aggregatedModel <- caretEnsemble(
      listM,
      trControl=trainControl(
        number=10,
        verboseIter = TRUE
      ))
    
    x<-as.numeric(aggregatedModel$error[,c("Rsquared","RMSE","MAE")])
    x<-c(cor(testing$group,predict(aggregatedModel,newdata = testing),method = "pearson"),x)
    ensemble<-x; names(ensemble)<-c("Correlation","Rsquared","RMSE","MAE")
    
    res<-cbind(ensemble,res)
    res<-res[,order(-res["Correlation",], -res["Rsquared",])]
    
    model_list[["ensemble"]]<-aggregatedModel
  }
  
  model_list<-model_list[colnames(res)]
  
  models.results<-list(res,model_list)
  names(models.results)<-c("stats","models")
  
  return(models.results)
  
}


prioML<-function(lista,nameplot=NULL){
  
  require("ggplot2")
  require("reshape")
  
  # Get best models..................
  RES<-matrix(data=0,ncol=9,nrow=13)
  colnames(RES)<-c("rf","nb","xgbTree","svmLinear","svmRadial","nnet","glm","knn","lda")
  rownames(RES)<-rownames(lista[[1]]$stats)
  RES<-as.data.frame(RES)
  
  for(perm in 1:length(lista)){
    for(i in 1:nrow(RES)){
      for(j in 1:ncol(RES)){
        if(is.null(lista[[perm]]$stats[i,colnames(RES)[j]])==FALSE){
          RES[i,j]<-RES[i,j]+lista[[perm]]$stats[i,colnames(RES)[j]]
        }
      }
    }
  }
  
  RES<-RES/length(lista)
  
  
  RES$stats<-rownames(RES)
  #print(RES)
  RES<-melt(RES)
  
  if(is.null(nameplot)==FALSE){
    tiff(filename = nameplot,width = 1500,height = 1200,res = 300, units = "in")
    p1<-ggplot(RES,aes(x=stats,y=value,group=variable,color=variable))+theme_bw()+
      geom_line()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      scale_color_manual(values=c("rf"="#FF6633","nb"="#990066",
                                  "xgbTree"="#3399CC","svmLinear"="#33FF99",
                                  "svmRadial"="#009900","nnet"="firebrick3",
                                  "glm"="#CC9900","knn"="#FFFF33","lda"="#330000"))+
      geom_hline(yintercept = 0.7,linetype="dashed",color="blue")+ylim(0,1)
    
    plot(p1)
    dev.off()
  }
  
  tmp<-RES[RES$stats=="AUC",]
  tmp<-tmp[order(tmp$value,decreasing=T),]
  
  #if(tmp$value[1]>=0.6){ ## Filtro de 0.6 AUC medio
  tmp<-as.character(tmp$variable[1])
  
  ## Get most hight model
  m<-data.frame(lista[[1]]$stats[,tmp])
  rownames(m)<-rownames(lista[[1]]$stats)
  
  for(i in 2:length(lista)){
    m<-cbind(m,data.frame(lista[[i]]$stats[,tmp]))
  }
  colnames(m)<-1:length(lista)
  
  m<-m[,order(-m[2,],-m[3,])]
  best<-as.numeric(as.character(colnames(m)[1]))
  
  model<-lista[[best]]$models[[tmp]]
  stats<-as.data.frame(m[,1])
  rownames(stats)<-rownames(m)
  colnames(stats)<-"stats"
  
  res<-list(model,stats)
  names(res)<-c("model","stats")
  
  return(res)
  
  #}else{
  
  #  return("NotSelected")
  #}
  
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






