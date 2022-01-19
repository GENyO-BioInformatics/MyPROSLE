##############################
## MyPROSLE 
## R version R version 4.0.4
## 
##############################
## Drug response

load("D:/DATA/WORK/Toro.et.al.MyPROSLE/RData/Clustering.RData")
set.seed(123456788)

pkgs<-c("tmod","caret","matrixStats","biomaRt","NOISeq",
        "NbClust","SNFtool","tidyr","parallel","ggplot2",
        "ConsensusClusterPlus","raster","pheatmap","reshape",
        "survival","dplyr","survminer","ggrepel","doParallel",
        "stringr")
check.packages(pkgs)
rm(pkgs)

rm(list=setdiff(ls(),c("MODELS","Module.colors","Modules.ann","Modules.list","opt")))

source(paste0(opt$scriptPath,"/Resources.R",sep=""))

## ................................................................STEP 1
## Loading Tabalumab dataset

load(paste0(opt$dataPath,sep="","/Tabalumab.RData"))

#-------
## Log Normalization
Raw.tab<-norm.log(Raw.tab)

#-------
## Annotation probe sets to gene symbol
rownames(Raw.tab)<-gsub("(.hg).*","\\1",rownames(Raw.tab))

## to avoid problems with biomaRt R package
httr::set_config(httr::config(ssl_verifypeer = FALSE))

## Run gene anotation to gene symbol
data<-annotateGenes(data=Raw.tab,
                    toGenes ='external_gene_name',
                    fromGenes = 'affy_hta_2_0')

rm(Raw.tab)

save.image(paste0(opt$dataPath,sep="","/Tabalumab.annotated.RData"))

#-------
## Separe Healthy and SLE samples
data<-data[,rownames(clin.tab)]

sel<-ifelse(clin.tab$State=="Normal",T,F)

H<-data[,sel]
SLE<-data[,!sel]

cat(paste0("\nNumber of controls: ",sep="",ncol(H)))

clin.tab<-clin.tab[!sel,]

#-------
## Correct terms in metadata
x<-clin.tab$State
x[x=="Systemic lupus erythematosus (SLE)"]<-"SLE"
clin.tab$State<-x; rm(x)
x<-clin.tab$Treatment
x[is.na(x)]<-"None"
clin.tab$Treatment<-x; rm(x)


##------
## Remove non informative samples
clin.tab<-clin.tab[ifelse(clin.tab$Time=="baseline" & clin.tab$Treatment!="None",F,T),]

##------
## Select patients with three visits
pats<-unique(clin.tab$SubjectID)
npats<-0
selected3<-NULL

for(i in 1:length(pats)){
  tmpM<-clin.tab[clin.tab$SubjectID==pats[i],]
  tmp<-as.character(tmpM$Time)
  
  if("baseline" %in% tmp & "week16" %in% tmp & "week52" %in% tmp){
    npats<-npats+1
    selected3<-c(selected3,rownames(tmpM))
  }
}
cat(paste0("\nNumero de pacientes: ",sep="",npats))
rm(npats,i,tmpM,tmp)
## 417

clin.tab<-clin.tab[selected3,]
SLE<-SLE[,selected3]


## ................................................................STEP 2
## Get M-scores for selected patients

Modules.list<-Modules.list[rownames(Modules.ann)]

SLE.mscores<-Get.MyPROSLE(modules.list = Modules.list, Patient=SLE,
                  Healthy = H, method="mean")

rm(sel,H,SLE,data)

SLE.msig<-GetSignatures(listM = SLE.mscores,ann = Modules.ann)


## ................................................................STEP 3
## Transcriptional response based on dysregulated signatures at baseline

base.samples<-clin.tab[ifelse(clin.tab$Time=="baseline",T,F),]
BASELINE<-SLE.msig[,rownames(base.samples)]
colnames(BASELINE)<-base.samples$SubjectID

##------
## Calculate delta between times and baseline
DELTA<-data.frame(matrix(ncol=5,nrow=0))
colnames(DELTA)<-c("GSM","SubjectID","Time","Treatment","Response")

for(i in 1:nrow(base.samples)){
  
  tmpclin<-clin.tab[clin.tab$SubjectID==base.samples[i,1],]
  tmp<-SLE.msig[,c(rownames(tmpclin)[1],rownames(tmpclin)[2],rownames(tmpclin)[3])]
  tmp<-as.data.frame(tmp)
  
  
  if(nrow(tmp[tmp[,1]>=1.65,])>=1){
    
    t1<-(nrow(tmp[tmp[,2]>=1.65,])-nrow(tmp[tmp[,1]>=1.65,]))
    DELTA<-rbind(DELTA,c(colnames(tmp)[2],tmpclin[2,1],tmpclin[2,3],tmpclin[2,4],t1))
    
    t2<-(nrow(tmp[tmp[,3]>=1.65,])-nrow(tmp[tmp[,1]>=1.65,]))
    DELTA<-rbind(DELTA,c(colnames(tmp)[3],tmpclin[3,1],tmpclin[3,3],tmpclin[3,4],t2))
    
  }
  #t1<-M[,rownames(tmp)[2]]-M[,rownames(tmp)[1]]
  #t2<-M[,rownames(tmp)[3]]-M[,rownames(tmp)[1]]
}

colnames(DELTA)<-c("GSM","SubjectID","Time","Treatment","Response")
DELTA$Response<-as.numeric(as.character(DELTA$Response))

DELTA<-cbind(DELTA,t(BASELINE[,DELTA$SubjectID]))

colnames(DELTA)<-FixSymbol(vect=colnames(DELTA),symbol=c("/"," "),toSymbol=c(".",""))

##------
## Response based on number of dysregulated signatures

signatures<-colnames(DELTA)[6:ncol(DELTA)]

Events<-NULL

for(i in 1:length(signatures)){
  
  placebo<-NULL
  q4<-NULL
  q2<-NULL
  tmpDelta<-DELTA[ifelse(DELTA[,signatures[i]]>=1.65,T,F),c("Time","Treatment","Response",signatures[i])]
  
  if(nrow(tmpDelta)>0){
    ## No resp, Resp
    placebo<-c(
      sum(ifelse(tmpDelta$Treatment=="Placebo" & tmpDelta$Response>=0,1,0)),
      sum(ifelse(tmpDelta$Treatment=="Placebo" & tmpDelta$Response<0,1,0)))
    q4<-c(
      sum(ifelse(tmpDelta$Treatment=="LY 120 mg Q4W" & tmpDelta$Response>=0,1,0)),
      sum(ifelse(tmpDelta$Treatment=="LY 120 mg Q4W" & tmpDelta$Response<0,1,0)))
    q2<-c(
      sum(ifelse(tmpDelta$Treatment=="LY 120 mg Q2W" & tmpDelta$Response>=0,1,0)),
      sum(ifelse(tmpDelta$Treatment=="LY 120 mg Q2W" & tmpDelta$Response<0,1,0)))
    
    
    m<-rbind(placebo,q4,q2,signatures[i])
    m<-t(m)
    rownames(m)<-c("Nonresponse","Response")
    colnames(m)<-c("Placebo","Q4","Q2","Signature")
    #print(colnames(DELTA)[i])
    #print(m)
    
    ## We select signatures with at least 5 response and non-response events
    if(sum(as.numeric(m[1,1:3]))>=5 & sum(as.numeric(m[2,1:3]))>=5){
      Events<-rbind(Events,m)
    }
  }
  
}

Events<-cbind(rownames(Events),Events)
rownames(Events)<-NULL
Events<-as.data.frame(Events)
Events$Placebo<-as.numeric(as.character(Events$Placebo))
Events$Q4<-as.numeric(as.character(Events$Q4))
Events$Q2<-as.numeric(as.character(Events$Q2))
colnames(Events)[1]<-"Response"

##------
## Get significance (Fisher test / Hypergeometric test)

signatures<-unique(Events$Signature)
pvals<-NULL

for(i in 1:length(signatures)){
  
  x<-ifelse(Events$Signature==signatures[i],T,F)
  test<-t(as.data.frame(cbind(Events[x,"Placebo"],Events[x,"Q4"]+Events[x,"Q2"])))
  rownames(test)<-c("placebo","treatment")
  colnames(test)<-c("No","Resp")
  test<-t(as.data.frame(test))
  #print(signatures[i])
  print(test)
  
  q<-test["No","placebo"]
  m<-sum(test["No",])
  n<-sum(test)-m
  k<-sum(test[,"placebo"])
  
  rest<-phyper(q,m,n,k,lower.tail = F)
  pvals<-c(pvals,rest)
  
  #rest<-fisher.test(test,alternative = "greater")
  #pvals<-c(pvals,rest$p.value)
}
names(pvals)<-signatures


##------
## Plots

PlotList<-list()


for(i in 1:length(signatures)){
  
  ## Placebo
  x<-ifelse(Events$Signature==signatures[i],T,F)
  m<-Events[x,c("Response","Placebo")]
  colnames(m)<-c("Type","Value")
  m$Value<-as.numeric(as.character(m$Value))
  
  p1<-ggplot(m, aes(x="", y=Value, fill=Type)) + theme(axis.title.x=element_blank(),
                                                   axis.text.x=element_blank(),
                                                   axis.ticks.x=element_blank(),
                                                   panel.background = element_blank(),
                                                   plot.title = element_text(size = 8, face = "bold")) +
    geom_bar(stat="identity", width=1,color="black",alpha=0.7) +
    coord_polar("y", start=0)+ggtitle(signatures[i],)+
    scale_fill_manual(values=c("#455a64", "#E69F00"))+
    geom_text(x=1.66, y=0, 
              label=paste0("Pvalue: ",sep="",round(pvals[i],digits=4)),
              size=3, angle = 0,color="black", check_overlap = TRUE)
  
  ## Treatment
  m<-Events[x,c("Response","Q4","Q2")]
  m$Q4<-as.numeric(as.character(m$Q4))
  m$Q2<-as.numeric(as.character(m$Q2))
  m$Value<-m$Q4+m$Q2; m<-m[,c("Response","Value")]
  colnames(m)<-c("Type","Value")
  
  p2<-ggplot(m, aes(x="", y=Value, fill=Type)) + theme(axis.title.x=element_blank(),
                                                   axis.text.x=element_blank(),
                                                   axis.ticks.x=element_blank(),
                                                   panel.background = element_blank()) +
    geom_bar(stat="identity", width=1,color="black",alpha=0.7) +
    coord_polar("y", start=0)+
    scale_fill_manual(values=c("#455a64", "#E69F00"))+
    geom_text(x=1.66, y=0, 
              label=paste0("Pvalue: ",sep="",round(pvals[i],digits=4)),
              size=3, angle = 0,color="black", check_overlap = TRUE)
  
  pl<-list(p1,p2); names(pl)<-c("Placebo","Treatment")
  
  PlotList[[signatures[i]]]<-pl
  
}

## Piechart represets the % of patients with a High dysregulated signature at baseline, that respond and not respond

dev.off()

## Save plot

tiff(filename=paste0(opt$resultsPath,sep="","/TranscriptionalResponse_.tiff"),res = 300,width = 9.2,height = 3.5,units="in")
p<-ggarrange(PlotList[[1]]$Placebo,PlotList[[2]]$Placebo,PlotList[[3]]$Placebo,
             PlotList[[4]]$Placebo,PlotList[[5]]$Placebo,PlotList[[6]]$Placebo,
             PlotList[[1]]$Treatment,PlotList[[2]]$Treatment,PlotList[[3]]$Treatment,
             PlotList[[4]]$Treatment,PlotList[[5]]$Treatment,PlotList[[6]]$Treatment,
             ncol=6,nrow=2,label.x = c("Placebo","Tabalumab"),common.legend = TRUE,
             heights = c(1,1),widths = c(1,1))

print(p)
dev.off()



## ................................................................STEP 4
## Get classifier for Response / Non-response to tabalumab


DELTA2<-DELTA[ifelse(DELTA$Treatment!="Placebo",T,F),]
pats<-unique(DELTA2$SubjectID)

##------
## Select consistent response or non-response in the two times
patsSel<-NULL
for(i in 1:length(pats)){
  tmp<-DELTA2[DELTA2$SubjectID==pats[i],]
  
  if((tmp$Response[1]>=0 & tmp$Response[2]>=0) |
     (tmp$Response[1]<0 & tmp$Response[2]<0)){
    patsSel<-c(patsSel,tmp$SubjectID[1])
  }
}
DELTA2<-DELTA2[DELTA2$SubjectID %in% patsSel,]

##------
## Assign patient as responders (delta<0) or non-responders (delta(>=0))
x<-1:nrow(DELTA2)
DELTA2<-DELTA2[(x%%2)==0,]

DELTA2$Response<-ifelse(DELTA2$Response>=0,"NoResp","Resp")


#pheatmap(t(DELTA2[,6:14]),show_rownames = T,show_colnames = F, 
#         cluster_rows = T,cluster_cols = T,
#         breaks=seq(-2.5,2.5,length.out = 100),
#         annotation_col = DELTA2[,c(5,5)],
#         border_color = "black",
#         color = colorRampPalette(c("deepskyblue4","white","coral2"))(100))


#pheatmap(t(DELTA2[,c(6,9,10,12)]),show_rownames = T,show_colnames = F, 
#         cluster_rows = T,cluster_cols = T,
#         breaks=seq(-2.5,2.5,length.out = 100),
#         annotation_col = DELTA2[,c(5,5)],
#         border_color = "black",
#         color = colorRampPalette(c("deepskyblue4","white","coral2"))(100))


M<-DELTA2[,5:ncol(DELTA2)]
colnames(M)[1]<-"Group"

tabalumab.model<-GetModel(data=M,method="rf",bootstrap=0.8,kfold=10,repeatedCv=30,
                    perm.bias=10,featureSets=c(1,2,3,4,5,6,7,8,9),seed=12345,
                    accuracy.type ="Balanced Accuracy")

##-------
## Plot

library("caret")

M<-DELTA2[,c(5:14)]

res<-predict(tabalumab.model$model, newdata = M[,2:ncol(M)],type = "prob")
res1<-ifelse(res[,2]>=0.5, "Resp", "NoResp")

M<-cbind(res1,res$Resp,M); colnames(M)[1:3]<-c("Prediction","Probability","Response")


ann.colorResp<-list(Prediction=c(Resp="#eebc4d",NoResp="#7d8c93"),
                       Response=c(Resp="#eebc4d",NoResp="#7d8c93"))

tiff(filename=paste0(opt$resultsPath,sep="","/Tabalumab_Response.tiff"),res = 300,width = 9.2,height = 3,units="in")
p2<-pheatmap(t(M[,-c(1:3)]),show_rownames = T,show_colnames = F, 
             cluster_rows = T,cluster_cols = T,
             breaks=seq(-2.5,2.5,length.out = 100),
             annotation_col = M[,c(2,3)],
             border_color = "black",
             annotation_colors = ann.colorResp,
             color = colorRampPalette(c("deepskyblue4","white","coral2"))(100))
print(p2)
dev.off()

rm(ann.colorResp,base.samples,DATA,m,M,p,p1,p2,pl,PlotList,res,res1,test,tmp,tmpclin,
   tmpDelta,i,k,n,paths,pats,patsSel,placebo,pvals,q,q2,q4,rest,selected3,t1,t2,x)

save.image(paste0(opt$rdataPath,sep="","/DrugRespose.RData"))

##-------
## Save model

#opt<-AddOpt(option.list = opt,name.option = "modelPath",variable.value = paste0(opt$mainPath,"/Models"))

load(paste0(opt$modelPath,sep="","/Models.RData"))

MODELS<-AddModel(newModel = tabalumab.model,models=MODELS,name = "Tabalumab")

# "D:/DATA/WORK/Toro.et.al.MyPROSLE/Models/Models.RData"
rm(list=setdiff(ls(),c("MODELS")))

setwd("D:/DATA/WORK/Toro.et.al.MyPROSLE/Models")

save.image("Models.RData")

rm(list=ls())


## ................................................................STEP 4
## Get specific differences in gene-modules

load("D:/DATA/WORK/Toro.et.al.MyPROSLE/RData/DrugRespose.RData")
set.seed(123456788)

pkgs<-c("tmod","caret","matrixStats","biomaRt","NOISeq",
        "NbClust","SNFtool","tidyr","parallel","ggplot2",
        "ConsensusClusterPlus","raster","pheatmap","reshape",
        "survival","dplyr","survminer","ggrepel","doParallel",
        "stringr")
check.packages(pkgs)
rm(pkgs)

pats<-unique(clin.tab$SubjectID)

Response<-data.frame(matrix(ncol=4,nrow=0))
colnames(Response)<-c("D1","Time","Type","Module")

for(p in 1:length(pats)){ # Patient loop
  
  tmp<-clin.tab[clin.tab$SubjectID==pats[p],]
  
  if("baseline" %in% tmp$Time & "week16" %in% tmp$Time){ ## first time
    tmp.data<-SLE.mscores[,c(rownames(tmp)[tmp$Time %in% "baseline" & tmp$Treatment=="None"],rownames(tmp)[tmp$Time %in% "week16"])]
    
    for(m in 1:nrow(tmp.data)){ ## module loop
      
      if(tmp.data[m,1]>=1.65){
        d1<-tmp.data[m,2]-tmp.data[m,1]
        d1<-c(d1,"week16",as.character(tmp[2,4]),rownames(tmp.data)[m])
        Response<-rbind(Response,d1)
      }
    }
  } 
  
  if("baseline" %in% tmp$Time & "week52" %in% tmp$Time){ ## second time
    tmp.data<-SLE.mscores[,c(rownames(tmp)[tmp$Time %in% "baseline" & tmp$Treatment=="None"],rownames(tmp)[tmp$Time %in% "week52"])]
    
    for(m in 1:nrow(tmp.data)){ ## gene-module loop
      
      if(tmp.data[m,1]>=1.65){
        d1<-tmp.data[m,2]-tmp.data[m,1]
        d1<-c(d1,"week52",as.character(tmp[2,4]),rownames(tmp.data)[m])
        Response<-rbind(Response,d1)
      }
    }
  } 
  
} 

colnames(Response)<-c("D1","Time","Type","Module")
Response$D1<-as.numeric(Response$D1)

##-------
## Get plot and significance for all group comparisons

boxplotList<-list()
int<-1

modules<-unique(Response$Module)

for(m in 1:length(modules)){
  
  M<-Response[ifelse(Response$Module==modules[m],T,F),]
  M[M$Type=="Placebo","Type"]<-"0_Placebo"
  M[M$Type=="LY 120 mg Q2W","Type"]<-"2_120mgQ2W"
  M[M$Type=="LY 120 mg Q4W","Type"]<-"1_120mgQ4W"
  
  M$Type<-paste0(M$Type,"_",M$Time)
  
  combinations<-combn(x=unique(M$Type),m = 2)
  
  filter<-NULL
  comb<-NULL
  for(i in 1:ncol(combinations)){
    
    x<-M[M$Type==combinations[1,i],1]
    y<-M[M$Type==combinations[2,i],1]
    
    if(length(x)>5 & length(y)>5){
      
      test<-t.test(x,y)
      if(test$p.value<=0.01){
        
        title<-paste0(as.character(Modules.ann[modules[m],"category"]),sep=" ",as.character(Modules.ann[modules[m],"ID"]))
        colori<-as.character(Modules.ann[modules[m],"color"])
        
        filter<-c(filter,combinations[1,i],combinations[2,i])
        
        comb<-c(comb,i)
        
        colnames(M)[1]<-"delta"
        
      } # est$p.value<=0.05
    } # length(x)>=3 & length(y)>=3)
  } # combinations
  
  if(!is.null(filter)){
    filter<-unique(filter)
    sel<-M$Type %in% filter
    
    M2<-M[sel,]
    my_comparisons<-list()
    combinations<-combinations[,comb]
    if(length(comb)==1){
      my_comparisons[[1]]<-combinations
    }else{
      for(cb in 1:ncol(combinations)){
        my_comparisons[[cb]]<-c(combinations[1,cb],combinations[2,cb])
      }
    }
    
    gplot<-ggplot(M2,aes(x=Type,y=delta))+theme_classic()+ggtitle(title)+
      geom_jitter(color=colori,alpha=0.6,size=2.5) + 
      geom_boxplot(fill=colori,alpha=0.9,outlier.alpha = 0) +
      theme(plot.title = element_text(size = 10),axis.text.x = element_text(angle = 90)) + 
      stat_compare_means(comparisons = my_comparisons,method="t.test",size = 3.5)
    
    #print(gplot)
    boxplotList[[int]]<-gplot
    int<-int+1
  }
  
} 

##--------
## Make plot

order<-c(1,9,10,
         11,12,13,
         4,6,2,3,
         5,8,
         7,14)

boxplotList<-boxplotList[order]

lengths<-NULL
for(i in 1:length(boxplotList)){
  lengths<-c(lengths,length(unique(boxplotList[[i]]$data$Type)))
}

library("cowplot")

library("gridExtra")

tiff(filename=paste0(opt$resultsPath,sep="","/Boxplots_tabalumab_.tiff"),res = 300,width = 15,height = 15,units="in")

p<-grid.arrange(arrangeGrob(boxplotList[[1]],boxplotList[[2]],boxplotList[[3]],boxplotList[[4]],ncol=4),
                arrangeGrob(boxplotList[[5]],boxplotList[[6]],boxplotList[[7]],boxplotList[[8]],ncol=4),
                arrangeGrob(boxplotList[[9]],boxplotList[[10]],boxplotList[[11]],boxplotList[[12]],boxplotList[[13]],boxplotList[[14]],ncol=6),
                nrow = 3)       
print(p)
dev.off()

##--------
## Save Object
save.image(paste0(opt$rdataPath,sep="","/DrugRespose.RData"))

