###########################################
## Load Tabalumab dataset
## 
###########################################
##

##------------------------------------------------------------------------- STEP 0
## Set functions and environment 

## Define main work folder
pathMain<-"/home/daniel/Desktop/WORK/tabalumab/TmabData"
pathData<-"/home/daniel/Desktop/WORK/tabalumab/TmabData/Data"

set.seed(1234567)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)

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

packages<-c("affy","oligo","limma","GEOquery","stringr","biomaRt","matrixStats","parallel","car","caret","lars","glmnet")
check.packages(packages)
rm(packages)


##------------------------------------------------------------------------- STEP 1
## Load data 

setwd(pathMain)

## Get metadata from NCBI GEO
gset = getGEO("GSE88887", GSEMatrix =TRUE, destdir=getwd())
clin<-phenoData(gset[[1]])
clin<-pData(clin)

## Filter metadata columns
clin<-clin[,c("subject_id:ch1","group:ch1","time:ch1","treatment (q2w=once in 2 weeks; q4w=once in 4 weeks):ch1")]
colnames(clin)<-c("SubjectID","State","Time","Treatment")

## Correct terms in metadata
x<-clin$State
x[x=="Systemic lu pus erythematosus (SLE)"]<-"SLE"
clin$State<-x; rm(x)
x<-clin$Treatment
x[is.na(x)]<-"None"
clin$Treatment<-x; rm(x)

setwd(pathData)

selected<-rownames(clin)

## Get name of the files of the selected samples
celfiles <- list.files(pathData) # full = TRUE
nameFiles<-NULL
for(i in 1:length(selected)){
  x<-str_detect(celfiles,selected[i])
  x<-celfiles[x]
  if(!is.null(x) & !is.na(x)){
    nameFiles<-c(nameFiles,x)
  }
}

x<-setdiff(celfiles,nameFiles)
file.remove(x)
rm(gset,celfiles,i,selected,x)

#-------
## Read raw data from folder in subsetds (too huge dataset)

subset.list<-list()
i=0
subset<-1
while(i<length(nameFiles)){
  out<-i+200
  if(out>length(nameFiles)){
    out<-length(nameFiles)
  }
  sub<-nameFiles[i:out]
  subset.list[[subset]]<-sub
  i<-out+1
  subset<-subset+1
}
rm(i,out,sub,subset)

RAW<-list()
for(i in 1:length(subset.list)){
  
  files<-as.character(subset.list[[i]])
  raw<- read.celfiles(files)
  raw<-rma(raw)
  raw<-exprs(raw)
  RAW[[i]]<-raw
}
RAW = as.data.frame(do.call("cbind", RAW))
colnames(RAW)<-rownames(clin)

clin.tab<-clin
Raw.tab<-RAW

setwd(pathMain)
rm(RAW,clin,check.packages,load.packages,pathMain,subset.list,nameFiles,pathData,i,files,raw)

#-------
## Save RData
save.image("Tabalumab.RData")

