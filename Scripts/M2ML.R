##############################
## M2ML
## R version R version 4.2.0
## 
##############################
## 

## ................................................................
## Example
## Download example data from:
## https://drive.google.com/drive/folders/1C8qcAEN4wYtqpaWr8-3aL0RiUDXrkkUK?usp=sharing


setwd("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/M2ML")

##----------
## Input 

## 1) Functions
source("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/M2ML/Functions_M2ML.R")

## 2) List with Gene expression dataset (named as: exp) and vector of classes (named as: class)
datasets<-readRDS("datasets.rds")

## 3) Functional genesets (list with genes for each biological function)
load("ModulesInfo.RData")
geneset<-Modules.list
rm(Modules.ann,Modules.list)

## 4) Dataframe with clinical variables to build prediction models (rownames must be the same that colnames of the gene expression dataset (samples ids))
metadata<-readRDS("metadata.rds")


##----------
## Get predictive models
## M2ML variables
# exp.data: List with Gene expression dataset (named as: exp) and vector of classes (named as: class). Vector of classes: 0 for healthy, 1 for case samples. Only case samples could be selected (but a ref.panel must be provided)
# geneset: Functional genesets (list with genes for each biological function). Gene IDs from genesets and from datasets must be the same
# metadata: data.frame, samples in rows and clinical outcomes to predict in columns (rownames must be the same that colnames of the gene expression dataset (samples ids))
# var2predict: # var2predict: name of the variable from metadata to predict (colname)
# var.type: "cat" / "num" for categorical or numerical variables
# outerfolds: Number of outer train/test splits of the exp.data.
# prob: proportion train/test split (from 0 to 1). I.e. 0.8: 80% and 20% of the samples are randomly selected for the train set and test set, respectively 
# innerfold: number of folds for parameter tuning.
# repeats: number of iterations for parameter tuning
# ref.panel: optional (only if healthy samples are not included), path of RData with normalized samples to calculate M-scores by patient-patient similarity
# nk: optional (only if healthy samples are not included), number of most similar samples to impute mscore based on patient reference.


## Example 1: Categorical variable (YES/NO)
model.example1<-M2ML(exp.data=datasets, 
                     geneset=geneset,
                     metadata=metadata,
                     var2predict="var1.example",
                     var.type = "cat",
                     outerfolds=3,
                     prob=0.8,
                     innerfold=10,
                     repeats=3,
                     ref.panel=NULL,
                     nk=5)

## Example 2: Continuous variable
model.example2<-M2ML(exp.data=datasets,
                     geneset=geneset,
                     metadata=metadata,
                     var2predict="var2.example",
                     var.type = "num",
                     outerfolds=3,
                     prob=0.8,
                     innerfold=10,
                     repeats=3,
                     ref.panel=NULL,
                     nk=5)



## Example 3: Without healthy controls
sel<-ifelse(datasets$class==1,T,F)
dataset.noH<-list("exp"=datasets$exp[,sel],
                  "class"=datasets$class[sel])

model.example3<-M2ML(exp.data=dataset.noH,
                     geneset=geneset,
                     metadata=metadata,
                     var2predict="var2.example",
                     var.type = "num",
                     outerfolds=3,
                     prob=0.8,
                     innerfold=10,
                     repeats=3,
                     ref.panel="References.Rdata",
                     nk=3)




