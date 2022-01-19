##############################
## MyPROSLE 
## R version R version 4.0.4
## 
##############################
## Clinical associations

load("D:/DATA/WORK/Toro.et.al.MyPROSLE/RData/Clustering.RData")
set.seed(123456788)

pkgs<-c("tmod","caret","matrixStats","biomaRt","NOISeq",
        "NbClust","SNFtool","tidyr","parallel","ggplot2",
        "ConsensusClusterPlus","raster","pheatmap","reshape",
        "survival","dplyr","survminer","ggrepel","doParallel",
        "stringr")
check.packages(pkgs)
rm(pkgs)

rm(list=setdiff(ls(),c("Modules.ann","DATA.Mscore","DATA.Mscore","opt","Module.colors")))

source(paste0(opt$scriptPath,"/Resources.R",sep=""))

#opt<-AddOpt(option.list = opt,name.option = "modelPath",variable.value = paste0(opt$mainPath,"/Models"))

load(paste0(opt$modelPath,sep="","/Models.RData"))

GetModels<-FALSE

## ................................................................STEP 1
## Symptomathology  ## 24

## Load data
load(paste0(opt$dataPath,sep="","/Clin_pascualSLE.RData"))
load(paste0(opt$dataPath,sep="","/Clin_precisesadsSLE.RData"))
load(paste0(opt$dataPath,sep="","/Clin_petriSLE.RData"))
rm(M.pascual,M.petri,M.precisesads)

DATA.Msig<-GetSignatures(listM = DATA.Mscore,ann = Modules.ann)

Ms<-cbind(DATA.Mscore$dataset1$SLE,DATA.Mscore$dataset8$SLE,DATA.Mscore$dataset9$SLE)
Ms<-Ms[rownames(Modules.ann),]
Msig<-cbind(DATA.Msig$dataset1,DATA.Msig$dataset8,DATA.Msig$dataset9)

##-------
## Fix clinical info

Clin.pascual.tmp<-Clin.pascual[,c("Seizure","Psychosis","OrganicBrainSyndrome",
                                  "Visual_Disturbance","CranialNerve_Disorder","SLE_Headache",
                                  "Cva","Vasculitis","Arthritis","Myositis","UrinaryCast",
                                  "Hematuria","Proteinuria","Pyuria","NewRash","Alopecia",
                                  "Mucosal_ulcers","Pleurisy","Pericarditis","LowComplement",
                                  "Increased_DNAbinding","Fever","Thrombocytopenia","Leukopenia")]

for(i in 1:ncol(Clin.pascual.tmp)){
  Clin.pascual.tmp[,i]<-ifelse(is.na(Clin.pascual.tmp[,i]),NA,
                               ifelse(as.character(Clin.pascual.tmp[,i])=="0","NO","Yes"))
}

Clin.petri.visit.tmp<-Clin.petri.visit[,c("CNS.1","CNS.2","CNS.3","CNS.4","CNS.5",
                                          "CNS.6","CNS.7","Vasc.1","MS.1","MS.2",
                                          "Renal.1","Renal.2","Renal.3","Renal.4",
                                          "Skin.1","Skin.2","Skin.3","Sero.1","Sero.2",
                                          "Imm.1","Imm2","Cons.1","Heme.1","Heme.2")]

for(i in 1:ncol(Clin.petri.visit.tmp)){
  Clin.petri.visit.tmp[,i]<-ifelse(is.na(Clin.petri.visit.tmp[,i]),NA,
                                   ifelse(as.character(Clin.petri.visit.tmp[,i])=="0","NO","Yes"))
}

dataTMP<-list(Clin.pascual.tmp,Clin.petri.visit.tmp) #,Clin.precisesads.tmp

adj="fdr"

##-------
## Get significance and models for each variable

variable.list<-list()
variable.list[[1]]<-c("Seizure","CNS.1"); names(variable.list)[1]<-"CNS_Seizure"
variable.list[[2]]<-c("Psychosis","CNS.2"); names(variable.list)[2]<-"CNS_Psychosis"
variable.list[[3]]<-c("OrganicBrainSyndrome","CNS.3"); names(variable.list)[3]<-"CNS_OrganicBrainSyndrome"
variable.list[[4]]<-c("Visual_Disturbance","CNS.4"); names(variable.list)[4]<-"CNS_VisualDisturbance"
variable.list[[5]]<-c("CranialNerve_Disorder","CNS.5"); names(variable.list)[5]<-"CNS_CranialNerveDisorder"
variable.list[[6]]<-c("SLE_Headache","CNS.6"); names(variable.list)[6]<-"CNS_SLEHeadache"
variable.list[[7]]<-c("Cva","CNS.7"); names(variable.list)[7]<-"CNS_CerebrovascularAccident"
variable.list[[8]]<-c("Seizure","CNS.1","Psychosis","CNS.2","OrganicBrainSyndrome","CNS.3",
                      "Visual_Disturbance","CNS.4","CranialNerve_Disorder","CNS.5",
                      "SLE_Headache","CNS.6","Cva","CNS.7"); names(variable.list)[8]<-"CNS_Summary"
variable.list[[9]]<-c("Vasculitis","Vasc.1"); names(variable.list)[9]<-"Vascular_Vasculitis"
variable.list[[10]]<-c("Arthritis","MS.1"); names(variable.list)[10]<-"Musculoskeletal_Arthritis"
variable.list[[11]]<-c("Myositis","MS.2"); names(variable.list)[11]<-"Musculoskeletal_Myositis"
variable.list[[12]]<-c("Arthritis","MS.1","Myositis","MS.2"); names(variable.list)[12]<-"Musculoskeletal_Summary"
variable.list[[13]]<-c("UrinaryCast","Renal.1"); names(variable.list)[13]<-"Renal_UrinaryCast"
variable.list[[14]]<-c("Hematuria","Renal.2"); names(variable.list)[14]<-"Renal_Haematuria"
variable.list[[15]]<-c("Proteinuria","Renal.3"); names(variable.list)[15]<-"Renal_Proteinuria"
variable.list[[16]]<-c("Pyuria","Renal.4"); names(variable.list)[16]<-"Renal_Pyuria"
variable.list[[17]]<-c("UrinaryCast","Renal.1","Hematuria","Renal.2",
                         "Proteinuria","Renal.3","Pyuria","Renal.4"); names(variable.list)[17]<-"Renal_Summary"
variable.list[[18]]<-c("NewRash","Skin.1"); names(variable.list)[18]<-"Dermal_Rash"
variable.list[[19]]<-c("Alopecia","Skin.2"); names(variable.list)[19]<-"Dermal_Alopecia"
variable.list[[20]]<-c("Mucosal_ulcers","Skin.3"); names(variable.list)[20]<-"Dermal_MucosalUlcers"
variable.list[[21]]<-c("NewRash","Skin.1","Alopecia","Skin.2","Mucosal_ulcers","Skin.3"); names(variable.list)[21]<-"Dermal_Summary"
variable.list[[22]]<-c("Pleurisy","Sero.1"); names(variable.list)[22]<-"Serosal_Pleurisy"
variable.list[[23]]<-c("Pericarditis","Sero.2"); names(variable.list)[23]<-"Serosal_Pericarditis"
variable.list[[24]]<-c("Pleurisy","Sero.1","Pericarditis","Sero.2"); names(variable.list)[24]<-"Serosal_Summary"
variable.list[[25]]<-c("LowComplement","Imm.1"); names(variable.list)[25]<-"Immunological_LowComplement"
variable.list[[26]]<-c("Increased_DNAbinding","Imm2"); names(variable.list)[26]<-"Immunological_IncreasedDNAbinding"
variable.list[[27]]<-c("LowComplement","Imm.1","Increased_DNAbinding","Imm2"); names(variable.list)[27]<-"Immunological_Summary"
variable.list[[28]]<-c("Fever","Cons.1"); names(variable.list)[28]<-"Constitutional_Fewer"
variable.list[[29]]<-c("Thrombocytopenia","Heme.1"); names(variable.list)[29]<-"Haematological_Trombocytopenia"
variable.list[[30]]<-c("Leukopenia","Heme.2"); names(variable.list)[30]<-"Haematological_Leukopenia"
variable.list[[31]]<-c("Thrombocytopenia","Heme.1","Leukopenia","Heme.2"); names(variable.list)[31]<-"Haematological_Summary"

## 24 Clinical manifestations

RESULTS<-data.frame(matrix(ncol=0,nrow=nrow(Modules.ann)))

##-------
## Get significance by anova and classifier for each clinical variable
for(i in 1:length(variable.list)){

  cat(paste0("\nGetting results from ",sep="",names(variable.list[i])))
  vars<-as.character(variable.list[[i]])
  res<-ClinicalAnalyser(variables=vars,M=Ms,data=dataTMP,adjust = adj)
  
  RESULTS<-cbind(RESULTS,res$modules)
  colnames(RESULTS)[ncol(RESULTS)]<-names(variable.list[i])
  
  if(!is.null(res$melt)){

    if(GetModels){
      if(length(table(res$melt$Group))>1){ ## Two categories
        if(table(res$melt$Group)[1]>=10 & table(res$melt$Group)[2]>=10){ ## At least 10 points for category
         
          tmp.model<-GetModel(data=res$melt,method="rf",bootstrap=0.8,kfold=10,repeatedCv=20,
                      perm.bias=10,featureSets=c(1,2,3,4,5,6,7,8,9,10,15,20,25,30),seed=12345,
                      accuracy.type ="Balanced Accuracy")
        
          MODELS<-AddModel(newModel = tmp.model,models=MODELS,name = names(variable.list)[i])
          print(MODELS[[length(MODELS)]])
        }
      }
    }
  }
  
}

rownames(RESULTS)<-rownames(Modules.ann)

##-------
## Parameters for visualization
hits<-c(-100,-50,-30,-20,-10,-3,-1.2,0,1.2,3,10,20,40,50,100)
colorseq<-c("#000033","#000033","#000099","#0000CC","#3366CC","#3399FF",
            "white","white","white","#FFCC33","#FF6633","#FF3333","#993333",
            "#990033","#FF3399")

bks<-Custom.colorRange(hits=hits,colorseq=colorseq,size=1000)



## Filtering non significant results
TMP<-abs(RESULTS)
RESULTS.filter.sym<-RESULTS[ifelse(apply(TMP,1,max)>1.3,T,F),
                        ifelse(apply(TMP,2,max)>1.3,T,F)]

ann<-data.frame(Modules.ann$path10)
rownames(ann)<-rownames(Modules.ann)
colnames(ann)<-"Signature"


tiff(filename=paste0(opt$resultsPath,sep="","/Sympthoms_.tiff"),res = 300,width = 9,height = 3,units="in")
p<-pheatmap(t(RESULTS.filter.sym),show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = T, fontsize = 7,
         breaks=bks$breaks,color=bks$colors,border_color = "black",
         annotation_col = ann,annotation_colors = Module.colors)
print(p)
dev.off()
#
rm(Clin.pascual,Clin.pascual.tmp,Clin.petri.patient,Clin.petri.visit,
   Clin.petri.visit.tmp,Clin.precisesads,dataTMP,TMP,
   bks,i,res,hits,vars,colorseq,variable.list,p)

## ................................................................STEP 2
## Auto-antibodies (PRECISESADS data)  ## 23
##

## Load data
clin.aa<-read.csv(paste0(opt$dataPath,sep="","//auto_antibodies.txt"),header=T,sep="\t",dec=",")
rownames(clin.aa)<-clin.aa$OMICID

sel<-ifelse(str_detect(string = colnames(clin.aa),pattern = "CALL",negate = F),T,F)

clin.aa<-clin.aa[intersect(colnames(Ms),rownames(clin.aa)),sel]

for(i in 1:ncol(clin.aa)){
  clin.aa[,i]<-ifelse(is.na(clin.aa[,i]),NA,
                      ifelse(as.character(clin.aa[,i])=="negative","NO","Yes"))
}

dataTMP<-list(clin.aa)

RESULTS.aa<-data.frame(matrix(ncol=0,nrow=nrow(Modules.ann)))

##-------
## Get significance by anova and classifier for each clinical variable
for(i in 1:ncol(clin.aa)){
  cat(paste0("\nGetting results from ",sep="",colnames(clin.aa)[i]))
  
  vars<-c(colnames(clin.aa)[i]) 
  res<-ClinicalAnalyser(variables=vars,M=Ms,data=dataTMP,adjust = adj)
  RESULTS.aa<-cbind(RESULTS.aa,res$modules)
  colnames(RESULTS.aa)[ncol(RESULTS.aa)]<-colnames(clin.aa)[i]
 
  if(GetModels){
    if(!is.null(res$melt)){
    
      if(length(table(res$melt$Group))>1){ ## Two categories
        if(table(res$melt$Group)[1]>=10 & table(res$melt$Group)[2]>=10){ ## At least 10 points for category
        
          tmp.model<-GetModel(data=res$melt,method="rf",bootstrap=0.8,kfold=10,repeatedCv=20,
                            perm.bias=10,featureSets=c(1,2,3,4,5,6,7,8,9,10,15,20,25,30),seed=12345)
        
          MODELS<-AddModel(newModel = tmp.model,models=MODELS,name = colnames(clin.aa)[i])
        }
      }
    }
  }

}
rownames(RESULTS.aa)<-rownames(Modules.ann)


##-------
## Parameters for visualization
hits<-c(-100,-50,-30,-20,-10,-3,-1.2,0,1.2,3,10,20,40,50,100)
colorseq<-c("#000033","#000033","#000099","#0000CC","#3366CC","#3399FF",
            "white","white","white","#FFCC33","#FF6633","#FF3333","#993333",
            "#990033","#FF3399")

bks<-Custom.colorRange(hits=hits,colorseq=colorseq,size=1000)

## Filtering non significant results
TMP<-abs(RESULTS.aa)
RESULTS.filter.aa<-RESULTS.aa[ifelse(apply(TMP,1,max)>1.3,T,F),
                        ifelse(apply(TMP,2,max)>1.3,T,F)]



ann<-data.frame(Modules.ann$path10)
rownames(ann)<-rownames(Modules.ann)
colnames(ann)<-"Signature"

tiff(filename=paste0(opt$resultsPath,sep="","/Autoantibodies_.tiff"),res = 300,width = 9,height = 2.5,units="in")
p<-pheatmap(t(RESULTS.filter.aa),show_rownames = T,show_colnames = F, 
         cluster_rows = T,cluster_cols = T, fontsize = 7,
         breaks=bks$breaks,color=bks$colors,border_color = "black",
         annotation_col = ann,annotation_colors = Module.colors)
print(p)
dev.off()
#

rm(bks,clin.aa,dataTMP,p,res,TMP,adj,colorseq,hits,i,sel,vars)



## ................................................................STEP 3
## Cytokines (PRECISESADS data) ## 18 

##-------
## Load data
clin.cyt<-read.csv(paste0(opt$dataPath,sep="","/cytokines.txt"),header=T,sep="\t",dec=".")
rownames(clin.cyt)<-clin.cyt$OMICID

clin.cyt<-clin.cyt[intersect(colnames(Ms),rownames(clin.cyt)),c(106:ncol(clin.cyt))]


##-------
## Filtering cytokines
datos<-NULL
for(i in 1:ncol(clin.cyt)){
  datos<-c(datos,table(is.na(clin.cyt[,i]))["FALSE"])
}
names(datos)<-colnames(clin.cyt)

opcion=2

if(opcion==1){ ## Remove only duplicates
  rmCol<-c("BLC_ELISA","CRP","FAS_LIGAND_ELISA","GDF_15_ELISA",
           "IL_1_RA_ELISA","IL_1_RII_ELISA","IL_6_CALL","IP_10_ELISA",
           "MCP_2_ELISA","MCP_4_ELISA","MIP_1_BETA_ELISA","MMP_2",
           "MMP_8_ELISA","TARC_ELISA","TNF_ALPHA","TNF_RI_ELISA")
  sel<-!(colnames(clin.cyt) %in% rmCol)
  
}else{ ## Select the most measured cytokines
  varSel<-c("BAFF_ELISA","BLC","CRP_ELISA","FAS_LIGAND","GDF_15",
            "IL_1_RA","IL_1_RII","IL_6_ELISA","IP_10","MCP_2","MCP_4",
            "MIP_1_BETA","MMP_2_ELISA","MMP_8","TARC","TGF_BETA_ELISA",
            "TNF_ALPHA_ELISA","TNF_RI")
  sel<-(colnames(clin.cyt) %in% varSel)
  
}

clin.cyt<-clin.cyt[,sel]
rm(datos,sel)


#GetModels<-TRUE
##-------
## Get Cytokine correlation
RESULTS_cyt<-matrix(data=NA,ncol=ncol(clin.cyt),nrow=nrow(Ms))
rownames(RESULTS_cyt)<-rownames(Ms)  
colnames(RESULTS_cyt)<-colnames(clin.cyt)

for(i in 1:ncol(clin.cyt)){
  
  values<-clin.cyt[,i]
  pats<-rownames(clin.cyt)[ifelse(is.na(values),F,T)]
  values<-values[ifelse(is.na(values),F,T)]
  
  M<-Ms[,pats]
  
  modulos<-NULL
  pvals<-NULL
  
  ## Get model
  if(GetModels){
    data<-cbind(values,t(Ms[,pats])) ## use Msig to avoid high correlated gene-modules
    data<-as.data.frame(data)

    tmp.model<-Getlm(data=data,bootstrap=0.8,kfold=10,repeatedCv=30,
                     perm.bias=10,max.features=15,seed=12345)
    
    
      MODELS<-AddModel(newModel = tmp.model,models=MODELS,name = colnames(clin.cyt)[i])
  }
  
  for(mod in 1:nrow(M)){
    yval<-as.numeric(M[mod,])
    
    if(sd(values)!=0 & sd(yval)!=0){
      cortest<-cor.test(values,yval,method = "pearson")
      modulos<-c(modulos,as.numeric(cortest$estimate))
      pvals<-c(pvals,cortest$p.value)
      
    }else{
      modulos<-c(modulos,0)
      pvals<-c(pvals,1)
    }
  }
  res<-cbind(modulos,pvals)
  colnames(res)<-c("Cor","Pval")  
  rownames(res)<-rownames(M)  
  res<-as.data.frame(res)
  
  res$Pval<-p.adjust(res$Pval,method="fdr")
  
  
  RESULTS_cyt[,i]<-(-log10(as.numeric(res$Pval)))*ifelse(as.numeric(res$Cor)>=0,1,-1)
  
}

##------
## Select models manually with an accuracy threshold
MODELS<-MODELS[c(1,2,3, 10, 11, 18, 19, 4, 6, 8, 9, 14)]
## 5 models: rsquared >0.5, accuracies > 0.7
## 5 models: accuracies >0.6


##-------
## Filter and plotting
TMP<-abs(RESULTS_cyt)
RESULTS.filter.cyt<-RESULTS_cyt[ifelse(apply(TMP,1,max)>1.3,T,F),
                           ifelse(apply(TMP,2,max)>1.3,T,F)]

## Parameters for visualization
hits<-c(-100,-50,-30,-20,-10,-3,-1.2,0,1.2,3,10,20,40,50,100)
colorseq<-c("#000033","#000033","#000099","#0000CC","#3366CC","#3399FF",
            "white","white","white","#FFCC33","#FF6633","#FF3333","#993333",
            "#990033","#FF3399")

bks<-Custom.colorRange(hits=hits,colorseq=colorseq,size=1000)

ann<-data.frame(Modules.ann$path10)
rownames(ann)<-rownames(Modules.ann)
colnames(ann)<-"Signature"

tiff(filename=paste0(opt$resultsPath,sep="","/Cytokines_.tiff"),res = 300,width = 9,height = 2.5,units="in")
p<-pheatmap(t(RESULTS.filter.cyt),show_rownames = T,show_colnames = F, 
            cluster_rows = T,cluster_cols = T, fontsize = 7,
            breaks=bks$breaks,color=bks$colors,border_color = "black",
            annotation_col = ann,annotation_colors = Module.colors)
print(p)
dev.off()
#


## ................................................................STEP 4
## Cell lines (PRECISESADS data)  ## 24

clin.cell<-read.csv(paste0(opt$dataPath,sep="","/FlowCytometry.csv"),header=T,sep=";",dec=".")
rownames(clin.cell)<-clin.cell$OMICID
clin.cell<-clin.cell[,-1]

clin.cell<-clin.cell[intersect(colnames(Ms),rownames(clin.cell)),]

sel<-ifelse(is.na(clin.cell$PBMC)==FALSE &is.na(clin.cell$PMN)==FALSE,T,F)

clin.cell<-clin.cell[sel,]

total<-clin.cell$PBMC + clin.cell$PMN

for(i in 1:ncol(clin.cell)){
  clin.cell[,i]<-clin.cell[,i]/total
}
clin.cell<-clin.cell[,-c(18,19,20,21,22)]

RESULTS_cell<-matrix(data=NA,ncol=ncol(clin.cell),nrow=nrow(Ms))
rownames(RESULTS_cell)<-rownames(Ms)  
colnames(RESULTS_cell)<-colnames(clin.cell)

for(i in 1:ncol(clin.cell)){
  
  values<-clin.cell[,i]
  pats<-rownames(clin.cell)[ifelse(is.na(values),F,T)]
  values<-values[ifelse(is.na(values),F,T)]
  
  M<-Ms[,pats]
  
  modulos<-NULL
  pvals<-NULL
  
  ## Get model
  if(GetModels){
    data<-cbind(values,t(Ms[,pats])) ## use Msig to avoid high correlated gene-modules
    data<-as.data.frame(data)
    
    tmp.model<-Getlm(data=data,bootstrap=0.8,kfold=10,repeatedCv=30,
                     perm.bias=10,max.features=15,seed=12345)
    
    
      MODELS<-AddModel(newModel = tmp.model,models=MODELS,name = colnames(clin.cell)[i])

  }
  
  for(mod in 1:nrow(M)){
    yval<-as.numeric(M[mod,])
    
    if(sd(values)!=0 & sd(yval)!=0){
      cortest<-cor.test(values,yval,method = "pearson")
      modulos<-c(modulos,as.numeric(cortest$estimate))
      pvals<-c(pvals,cortest$p.value)
      
    }else{
      modulos<-c(modulos,0)
      pvals<-c(pvals,1)
    }
  }
  res<-cbind(modulos,pvals)
  colnames(res)<-c("Cor","Pval")  
  rownames(res)<-rownames(M)  
  res<-as.data.frame(res)
  
  res$Pval<-p.adjust(res$Pval,method="fdr")
  
  
  RESULTS_cell[,i]<-(-log10(as.numeric(res$Pval)))*ifelse(as.numeric(res$Cor)>=0,1,-1)
  
}



##------
## Select models manually with an accuracy threshold
MODELS<-MODELS[c(1:12,13,15,18,20,21,22,24,26,29,14,16,17,23,25)]
# 1:12 (correctos)
# 13, 15, 18, 20, 21, 22, 24, 26, 29 #(9) rsq>=0.5
# 14, 16, 17, 23, 25,  ## (5) cor test >0.6 + rsq >0.2

#save.image(file=paste0(opt$mainPath,sep="","/temp.RData"))

##-------
## Parameters for visualization
hits<-c(-100,-50,-30,-20,-10,-3,-1.2,0,1.2,3,10,20,40,50,100)
colorseq<-c("#000033","#000033","#000099","#0000CC","#3366CC","#3399FF",
            "white","white","white","#FFCC33","#FF6633","#FF3333","#993333",
            "#990033","#FF3399")

bks<-Custom.colorRange(hits=hits,colorseq=colorseq,size=1000)

## Filtering non significant results
TMP<-abs(RESULTS_cell)
RESULTS.filter.cell<-RESULTS_cell[ifelse(apply(TMP,1,max)>1.3,T,F),
                              ifelse(apply(TMP,2,max)>1.3,T,F)]

ann<-data.frame(Modules.ann$path10)
rownames(ann)<-rownames(Modules.ann)
colnames(ann)<-"Signature"

tiff(filename=paste0(opt$resultsPath,sep="","/Cells_.tiff"),res = 300,width = 9,height = 2.5,units="in")
p<-pheatmap(t(RESULTS.filter.cell),show_rownames = T,show_colnames = F, 
            cluster_rows = T,cluster_cols = T, fontsize = 7,
            breaks=bks$breaks,color=bks$colors,border_color = "black",
            annotation_col = ann,annotation_colors = Module.colors)
print(p)
dev.off()
#
GetModels=FALSE

rm(clin.cell,clin.cyt,cortest,data,DATA.Mscore,DATA.Msig,M,Ms,Msig,p,res,RESULTS,
   RESULTS_cell,RESULTS_cyt,RESULTS.aa,TMP,tmp.model,colorseq,hits,i,mod,modulos,
   opcion,pats,pvals,sel,total,values,varSel,yval)

##------
## Save object with plots
save.image(file=paste0(opt$rdataPath,sep="","/ClinicalAssociations.RData"))


##------
## Save models

rm(list=setdiff(ls(),c("MODELS")))

setwd("D:/DATA/WORK/Toro.et.al.MyPROSLE/Models")

save.image("Models.RData")



