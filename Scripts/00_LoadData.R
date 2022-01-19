##############################
## MyPROSLE 
## R version R version 4.0.4
## 
##############################
## Load and process data

## For Github, Step 1 (preprocessing datasets) is avoid, go direct to Step2


## @@@@@@@@@@@@@@@@@@@@@@@@@@@ Input options
## Input general options
## database: "tmod","BloodGen3Module" or "all"

opt<-list(mainPath="D:/DATA/WORK/Toro.et.al.MyPROSLE",
          dataPath="D:/DATA/WORK/Toro.et.al.MyPROSLE/Datasets",
          scriptPath="D:/DATA/WORK/Toro.et.al.MyPROSLE/Scripts",
          database="tmod") 

set.seed(123456788)
source(paste0(opt$scriptPath,"/Resources.R",sep=""))

pkgs<-c("tmod","caret","matrixStats","biomaRt","NOISeq")
check.packages(pkgs)

rm(pkgs)

if (!file.exists(paste0(opt$mainPath,"/RData"))){
  dir.create(paste0(opt$mainPath,"/RData"))
} 
opt<-AddOpt(option.list = opt,name.option = "rdataPath",variable.value = paste0(opt$mainPath,"/RData"))

if (!file.exists(paste0(opt$mainPath,"/Results"))){
  dir.create(paste0(opt$mainPath,"/Results"))
} 
opt<-AddOpt(option.list = opt,name.option = "resultsPath",variable.value = paste0(opt$mainPath,"/Results"))

if (!file.exists(paste0(opt$mainPath,"/Models"))){
  dir.create(paste0(opt$mainPath,"/Models"))
} 
opt<-AddOpt(option.list = opt,name.option = "modelPath",variable.value = paste0(opt$mainPath,"/Models"))

## @@@@@@@@@@@@@@@@@@@@@@@@@@@ 

## ................................................................STEP 1
## Loading/ pre-processing data

#-------
## Gene expression datasets (downloaded from NCBI GEO and/or ADEX)
## ADEX preprocessing pipelines were previously applied to optain gene-expression matrices

DATA<-list()
cat("\nLoading and processing datasets...")

## Dataset 1: Pascual et. al (GSE65391)
data<-read.table(file=paste0(opt$dataPath,sep="","/GSE65391_expressionData.txt"),sep="\t",header=T,row.names = 1)
clin<-read.csv(file=paste0(opt$dataPath,sep="","/Pascual_allClin.csv"),sep="\t")
rownames(clin)<-clin[,1];  clin<-clin[,-1]
clin<-as.data.frame(clin); clin<-t(clin); clin<-as.data.frame(clin)
data<-data[,rownames(clin)];

## Preprocessing (Log normalization and filtering non variable genes)
data<-norm.log(data)
nonVar.genes<-summ.var(data = data)
data<-data[!nonVar.genes$nzv,] ## remove genes with near-zero variance

## Separate Healthy and SLE samples
state<-ifelse(grepl("BAY", clin$Patient),"Healthy","SLE")
H<-data[,state=="Healthy"]
SLE<-data[,state=="SLE"]

dataset1<-list(SLE,H); names(dataset1)<-c("SLE","Healthy")
DATA[["dataset1"]]<-dataset1
rm(clin,data,dataset1,H,SLE,state,nonVar.genes)


#-------
## DATASETS from ADEX
clin<-read.csv(file=paste0(opt$dataPath,sep="","/metadata.tsv"),sep="\t")
rownames(clin)<-clin$Sample

#-------
## Dataset 2: GSE45291_SLE
data<-read.csv(file=paste0(opt$dataPath,sep="","/GSE45291_SLE.tsv"),sep="\t")
rownames(data)<-data$gene; data<-data[,-1]

## Preprocessing (Log normalization and filtering non variable genes)
data<-norm.log(data)
nonVar.genes<-summ.var(data = data)
data<-data[!nonVar.genes$nzv,] ## remove genes with near-zero variance

## Separate Healthy and SLE samples
tmp<-clin[colnames(data),c("GSE","Condition")]
H<-data[,rownames(tmp)[tmp$Condition=="Healthy"]]
SLE<-data[,rownames(tmp)[tmp$Condition=="SLE"]]

dataset2<-list(SLE,H); names(dataset2)<-c("SLE","Healthy")
DATA[["dataset2"]]<-dataset2
rm(data,dataset2,H,SLE,tmp,nonVar.genes)

#-------
## Dataset 3: GSE61635
data<-read.csv(file=paste0(opt$dataPath,sep="","/GSE61635.tsv"),sep="\t")
rownames(data)<-data$gene; data<-data[,-1]

## Preprocessing (Log normalization and filtering non variable genes)
data<-norm.log(data)
nonVar.genes<-summ.var(data = data)
data<-data[!nonVar.genes$nzv,] ## remove genes with near-zero variance

## Separate Healthy and SLE samples
tmp<-clin[colnames(data),c("GSE","Condition")]
H<-data[,rownames(tmp)[tmp$Condition=="Healthy"]]
SLE<-data[,rownames(tmp)[tmp$Condition=="SLE"]]

dataset3<-list(SLE,H); names(dataset3)<-c("SLE","Healthy")
DATA[["dataset3"]]<-dataset3
rm(data,dataset3,H,SLE,tmp,nonVar.genes)

#-------
## Dataset 4: GSE72509
data<-read.csv(file=paste0(opt$dataPath,sep="","/GSE72509.tsv"),sep="\t")
rownames(data)<-data$gene; data<-data[,-1]

## Preprocessing (Log normalization and filtering non variable genes)
data<-norm.log(data)
nonVar.genes<-summ.var(data = data)
data<-data[!nonVar.genes$nzv,] ## remove genes with near-zero variance

## Separate Healthy and SLE samples
tmp<-clin[colnames(data),c("GSE","Condition")]
H<-data[,rownames(tmp)[tmp$Condition=="Healthy"]]
SLE<-data[,rownames(tmp)[tmp$Condition=="SLE"]]

dataset4<-list(SLE,H); names(dataset4)<-c("SLE","Healthy")
DATA[["dataset4"]]<-dataset4
rm(data,dataset4,H,SLE,tmp,nonVar.genes)

#-------
## Dataset 5: GSE108497
data<-read.csv(file=paste0(opt$dataPath,sep="","/GSE108497.tsv"),sep="\t")
rownames(data)<-data$gene; data<-data[,-1]

## Preprocessing (Log normalization and filtering non variable genes)
data<-norm.log(data)
nonVar.genes<-summ.var(data = data)
data<-data[!nonVar.genes$nzv,] ## remove genes with near-zero variance

## Separate Healthy and SLE samples
tmp<-clin[colnames(data),c("GSE","Condition")]
H<-data[,rownames(tmp)[tmp$Condition=="Healthy"]]
SLE<-data[,rownames(tmp)[tmp$Condition=="SLE"]]

dataset5<-list(SLE,H); names(dataset5)<-c("SLE","Healthy")
DATA[["dataset5"]]<-dataset5
rm(data,dataset5,H,SLE,tmp,nonVar.genes)

#-------
## Dataset 6: GSE110169_SLE
data<-read.csv(file=paste0(opt$dataPath,sep="","/GSE110169_SLE.tsv"),sep="\t")
rownames(data)<-data$gene; data<-data[,-1]

## Preprocessing (Log normalization and filtering non variable genes)
data<-norm.log(data)
nonVar.genes<-summ.var(data = data)
data<-data[!nonVar.genes$nzv,] ## remove genes with near-zero variance

## Separate Healthy and SLE samples
tmp<-clin[colnames(data),c("GSE","Condition")]
H<-data[,rownames(tmp)[tmp$Condition=="Healthy"]]
SLE<-data[,rownames(tmp)[tmp$Condition=="SLE"]]

dataset6<-list(SLE,H); names(dataset6)<-c("SLE","Healthy")
DATA[["dataset6"]]<-dataset6
rm(data,dataset6,H,SLE,tmp,nonVar.genes)

#-------
## Dataset 7: GSE110174
data<-read.csv(file=paste0(opt$dataPath,sep="","/GSE110174.tsv"),sep="\t")
rownames(data)<-data$gene; data<-data[,-1]

## Preprocessing (Log normalization and filtering non variable genes)
data<-norm.log(data)
nonVar.genes<-summ.var(data = data)
data<-data[!nonVar.genes$nzv,] ## remove genes with near-zero variance

## Separate Healthy and SLE samples
tmp<-clin[colnames(data),c("GSE","Condition")]
H<-data[,rownames(tmp)[tmp$Condition=="Healthy"]]
SLE<-data[,rownames(tmp)[tmp$Condition=="SLE"]]

dataset7<-list(SLE,H); names(dataset7)<-c("SLE","Healthy")
DATA[["dataset7"]]<-dataset7
rm(clin,data,dataset7,H,SLE,tmp,nonVar.genes)

#-------
## Dataset 8: PRECISESADS (2 sets of patients and 1 of healthy controls, without batch effect)

## First SLE set (and controls)
Count.precisesads<-read.table(file=paste0(opt$dataPath,sep="","/Count.PRECISESADS.csv"),
                              header=T,sep=";",row.names=1,dec=",",stringsAsFactors = F)
Count.precisesads<-type.convert(x = Count.precisesads,as.is=F)

## Filter non-expressed genes
Count.precisesads<-Count.precisesads[rownames(Count.precisesads)[rowCounts(Count.precisesads>=10)>=10],]

## Annotate to Gene symbol
Count.precisesads<-annotateGenes(data=Count.precisesads,
                        toGenes='external_gene_name',
                        fromGenes='ensembl_gene_id')

## Tmm and Log2 normalization and filtering non variable gene
data<-tmm(Count.precisesads)
data<-log2(data+1)
nonVar.genes<-summ.var(data = data)
data<-data[!nonVar.genes$nzv,] ## remove genes with near-zero variance

colnames(data)<-gsub(pattern = "X",replacement = "",x = colnames(data))

precisesads1<-data
rm(Count.precisesads,data,nonVar.genes)

## Second SLE set
Count.precisesads2<-read.table(file=paste0(opt$dataPath,sep="","/Count.PRECISESADS2.csv"),
                              header=T,sep=";",row.names=1,dec=",",stringsAsFactors = F)
Count.precisesads2<-type.convert(x = Count.precisesads2,as.is=F)

## Filter non-expressed genes
Count.precisesads2<-Count.precisesads2[rownames(Count.precisesads2)[rowCounts(Count.precisesads2>=10)>=10],]

## Annotate to Gene symbol
Count.precisesads2<-annotateGenes(data=Count.precisesads2,
                                 toGenes='external_gene_name',
                                 fromGenes='ensembl_gene_id')

## Tmm and Log2 normalization and filtering non variable gene
data<-tmm(Count.precisesads2)
data<-log2(data+1)
nonVar.genes<-summ.var(data = data)
data<-data[!nonVar.genes$nzv,] ## remove genes with near-zero variance

colnames(data)<-gsub(pattern = "X",replacement = "",x = colnames(data))

precisesads2<-data
rm(Count.precisesads2,data,nonVar.genes)

## Select common genes
genes<-as.character(intersect(rownames(precisesads1),rownames(precisesads2)))
precisesads1<-as.data.frame(precisesads1); precisesads1<-precisesads1[genes,]
precisesads2<-as.data.frame(precisesads2); precisesads2<-precisesads2[genes,]

## Select SLE and Healthy samples
clin1<-read.csv(file=paste0(opt$dataPath,sep="","/Metadata.PRECISESADS.csv"),header=T,sep=";",row.names=1,dec=",")
clin1<-clin1[colnames(precisesads1),c(2,8)]

SLE<-precisesads1[,rownames(clin1[ifelse(clin1$Diagnosis=="SLE",T,F),])]
H<-precisesads1[,rownames(clin1[ifelse(clin1$Diagnosis=="CTRL",T,F),])]  
  
clin2<-read.csv(file=paste0(opt$dataPath,sep="","/Metadata.PRECISESADS2.csv"),header=T,sep=";",row.names=1,dec=",")
clin2<-clin2[colnames(precisesads2),c(1,5)]

SLE2<-precisesads2[,rownames(clin2[ifelse(clin2$Diagnosis=="SLE",T,F),])]
SLE<-cbind(SLE,SLE2)

dataset8<-list(SLE,H); names(dataset8)<-c("SLE","Healthy")
DATA[["dataset8"]]<-dataset8
rm(clin1,clin2,dataset8,H,SLE,genes,SLE2,precisesads1,precisesads2)

#-------
## Dataset 9: Petri
data<-read.csv(file=paste0(opt$dataPath,sep="","/PetriALL.txt"),sep="\t")
rownames(data)<-data$GeneSymbol; data<-data[,-1]

## Preprocessing (Log normalization and filtering non variable genes)
data<-norm.log(data)
nonVar.genes<-summ.var(data = data)
data<-data[!nonVar.genes$nzv,] ## remove genes with near-zero variance

## Separate Healthy and SLE samples
colnames(data)<-gsub(pattern = "X",replacement = "",x = colnames(data))

clin<-read.csv(file=paste0(opt$dataPath,sep="","/Metadata.petri.csv"),sep=";")
rownames(clin)<-clin$GZ_Filenames

H<-data[,rownames(clin)[clin$Diagnosis=="Healthy"]]
SLE<-data[,rownames(clin)[clin$Diagnosis=="SLE"]]


dataset9<-list(SLE,H); names(dataset9)<-c("SLE","Healthy")
DATA[["dataset9"]]<-dataset9

rm(clin,data,dataset9,H,SLE,nonVar.genes)

data.info(data.list=DATA)

# Total: nºSLE (3048), nºH (697)

#-------
## Save RData object

save.image(paste0(opt$rdataPath,sep="","/Datasets.RData"))


## ................................................................STEP 2
## Get Mscores

#-------
## Select Gene Co-expression modules database
## Any database can be used (Must be transformed to pathway-named list of vector of genes)

#For Github: load Datasets.RData

switch (opt$database,
    tmod= {
      cat("\nUssing tmod databases...")
      
      data(tmod)
      Modules.list <- tmod$MODULES2GENES
      Modules.ann <- tmod$MODULES
      Modules.ann<-Modules.ann[,c(1,2)]
      colnames(Modules.ann)<-c("ID","Function")
      Modules.list<-Modules.list[rownames(Modules.ann)]
      rm(tmod)
      },
    BloodGen3Module = {
      cat("\nUssing BloodGen3Module database...")
      
      ## File extract from: https://github.com/Drinchai/BloodGen3Module/blob/master/R/sysdata.rda
      load(paste0(opt$dataPath,sep="","/sysdata.rda"))
      
      Modules.list<-list()
      modules<-unique(Gen3_ann$Module)
      for(i in 1:length(modules)){
        tmp<-Module_listGen3[ifelse(Module_listGen3$Module==modules[i],T,F),]
        Modules.list[[modules[i]]]<-as.character(tmp$Gene)
      }
      Modules.ann<-Gen3_ann[,c(1,2)]
      colnames(Modules.ann)<-c("ID","Function")
      Modules.list<-Modules.list[rownames(Modules.ann)]
      rm(tmp,modules,Gen3_ann,Module_listGen3,i,color)
    },
    all = {
    cat("\nUssing tmod and BloodGen3Module databases...")
    
    ## tmod
    data(tmod)
    Modules.list <- tmod$MODULES2GENES
    Modules.ann <- tmod$MODULES
    Modules.ann<-Modules.ann[,c(1,2)]
    colnames(Modules.ann)<-c("ID","Function")
    Modules.list<-Modules.list[rownames(Modules.ann)]
    rm(tmod)
    
    ## BloodGen3Module
    load(paste0(opt$dataPath,sep="","/sysdata.rda"))
    Modules.list2<-list()
    modules<-unique(Gen3_ann$Module)
    for(i in 1:length(modules)){
      tmp<-Module_listGen3[ifelse(Module_listGen3$Module==modules[i],T,F),]
      Modules.list2[[modules[i]]]<-as.character(tmp$Gene)
    }
    Modules.ann2<-Gen3_ann[,c(1,2)]
    colnames(Modules.ann2)<-c("ID","Function")
    Modules.list2<-Modules.list2[rownames(Modules.ann2)]
    rm(tmp,modules,Gen3_ann,Module_listGen3,i,color)
    
    Modules.list<-c(Modules.list,Modules.list2)
    Modules.ann<-rbind(Modules.ann,Modules.ann2)
    rm(Modules.ann2,Modules.list2)
  }
)

cat(paste0("\nNº of gene modules: ",sep="",nrow(Modules.ann)))

#-------
DATA.Mscore<-list()

for(i in 1:length(DATA)){ ## Get Mscores for all SLE samples and Healthy controls
  SLE<-Get.MyPROSLE(modules.list = Modules.list, Patient=DATA[[i]]$SLE,
                    Healthy = DATA[[i]]$Healthy, method="mean")
  H<-Get.MyPROSLE(modules.list = Modules.list, Patient=DATA[[i]]$Healthy,
                    Healthy = DATA[[i]]$Healthy, method="mean")
  
  x<-list(SLE,H); names(x)<-c("SLE","H")
  DATA.Mscore[[paste0("dataset",i,sep="")]]<-x
}

rm(SLE,H,i,x)

#-------
## Save RData object

save.image(paste0(opt$rdataPath,sep="","/Mscores.RData"))






