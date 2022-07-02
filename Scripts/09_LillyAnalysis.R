##############################
## MyPROSLE 
## R version R version 4.0.4
## 
##############################
## 10_LillyAnalysis

## ................................................................STEP 0
## Loading gene-module connections

set.seed(123456788)

load("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/RData/Clustering.RData")
rm(list=setdiff(ls(),c("Module.colors","Modules.ann","Modules.list","opt")))

setwd("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022")

source("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/scripts/Resources.R")

pkgs<-c("caret","matrixStats","biomaRt","sas7bdat","reshape","ggplot2")
check.packages(pkgs)
rm(pkgs)

## ................................................................STEP 1
## Loading expression dataset

load(paste0(getwd(),sep="","/Datasets/tabalumab/Tabalumab.RData"))

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

#-------
## Separe Healthy and SLE samples
data<-data[,rownames(clin.tab)]

sel<-ifelse(clin.tab$State=="Normal",T,F)

H<-data[,sel]
SLE<-data[,!sel]

cat(paste0("\nNumber of controls: ",sep="",ncol(H)))
# 60

clin.tab<-clin.tab[!sel,]

#-------
## Correct terms in metadata
x<-clin.tab$State
x[x=="Systemic lupus erythematosus (SLE)"]<-"SLE"
clin.tab$State<-x; rm(x)
x<-clin.tab$Treatment
x[is.na(x)]<-"None"
clin.tab$Treatment<-x; rm(x)

cat(paste0("\nNumero de pacientes: ",sep="",length(unique(clin.tab$SubjectID))))
# 1775

## ................................................................STEP 2
## Get M-scores for selected patients

Modules.list<-Modules.list[rownames(Modules.ann)]

SLE.mscores<-Get.MyPROSLE(modules.list = Modules.list, Patient=SLE,
                          Healthy = H, method="mean")

rm(sel,H,SLE,data)

SLE.msig<-GetSignatures(listM = SLE.mscores,ann = Modules.ann)


save.image("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/RData/TabalumabMscores.RData")

rm()

## ................................................................STEP 3
## Get clinical data

load("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/RData/TabalumabMscores.RData")

source("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/scripts/Resources.R")

## Primeras respuestas medidas...............
#clin.resp<-rbind(read.sas7bdat(paste0(getwd(),sep="","/Datasets/tabalumab/1/clinresp.sas7bdat")),
#         read.sas7bdat(paste0(getwd(),sep="","/Datasets/tabalumab/2/clinresp.sas7bdat")))
#clin.resp<-clin.resp[,c("USUBJID","TRT","MAJORCR","MAJRCRBG",
#        "MAJRCRSL","PARTLCR","PRTLCRBG","PRTLCRSL")]


## Medicion SRI......................
clin.resp<-rbind(read.sas7bdat(paste0(getwd(),sep="","/Datasets/tabalumab/1/sri.sas7bdat")),
                 read.sas7bdat(paste0(getwd(),sep="","/Datasets/tabalumab/2/sri.sas7bdat")))
clin.resp<-clin.resp[,c("USUBJID","TRT","VISID","SRI5RESP","SRI4RESP")]
# "SRI5RDUR","SRI5RRAW","SRI4RRAW"

patients<-unique(clin.resp$USUBJID)
Response<-as.data.frame(matrix(ncol=4,nrow=0))
colnames(Response)<-c("USUBJID","TRT","SRI5RESP","SRI4RESP")

for(i in 1:length(patients)){
  
  tmp<-clin.resp[clin.resp$USUBJID==patients[i],]
  if(16 %in% tmp$VISID){
    
    sel<-ifelse(tmp$VISID==16,T,F)
    Response<-rbind(Response,tmp[sel,c("USUBJID","TRT","SRI5RESP","SRI4RESP")])
    
  }
}

ids<-read.table(file=paste0(getwd(),sep="","/Datasets/tabalumab/1/Tabalumab_samples_clinical-usubjid.txt"),sep="\t",header = T)
# 1772 pacientes      


variable.list<-list() # "variable.type" "mscore" "msig" 

#vars<-c("MAJORCR","MAJRCRBG","MAJRCRSL","PARTLCR","PRTLCRBG","PRTLCRSL")
#vars<-c("SRI5RESP","SRI4RESP") #

vars<-list(list("var"="SRI5RESP","trt"="all"),
           list("var"="SRI5RESP","trt"="Q4W"),
           list("var"="SRI5RESP","trt"="Q2W"),
           list("var"="SRI4RESP","trt"="all"),
           list("var"="SRI4RESP","trt"="Q4W"),
           list("var"="SRI4RESP","trt"="Q2W"))


for(i in 1:length(vars)){
  
  # 699 patients treated (50:50) 
  # 327 placebo
  
  tmp<-Response[,c("USUBJID","TRT",vars[[i]]$var)]
  
  switch(vars[[i]]$trt,
         all={
           tmp<-tmp[ifelse(!is.na(tmp[,3]) & tmp$TRT!="Placebo",T,F),]
         },
         Q4W={
           tmp<-tmp[ifelse(!is.na(tmp[,3]) & tmp$TRT=="LY 120 mg Q4W",T,F),]
         },
         Q2W={
           tmp<-tmp[ifelse(!is.na(tmp[,3]) & tmp$TRT=="LY 120 mg Q2W",T,F),]
         })
  
   # 
  tmp[,3]<-ifelse(tmp[,3]==0,"NO","YES") 
  ## NO: 389, YES: 309  SRI5
  ## NO: 279, YES: 419  SRI4
  
  ids.tmp<-ids[ids$Time=="baseline",]
  ids.tmp<-ids.tmp[ids.tmp$Clinical_USUBJID %in% tmp$USUBJID,]
  
  response<-NULL
  for(j in 1:nrow(ids.tmp)){
   response<-c(response,
               tmp[tmp$USUBJID==ids.tmp$Clinical_USUBJID[j],vars[[i]]$var])
    
  }
  ids.tmp$response<-response
  
  mscores<-SLE.mscores[,ids.tmp$GEO_accession]
  msig<-SLE.msig[,ids.tmp$GEO_accession]
  
  
  mscores<-as.data.frame(t(mscores))
  msig<-as.data.frame(t(msig))
  
  mscores<-cbind(ids.tmp$response,mscores)
  msig<-cbind(ids.tmp$response,msig)
  colnames(mscores)[1]<-"group"
  colnames(msig)[1]<-"group"
  
  variable.type<-"categoric"
  
  res<-list(variable.type,mscores,msig)
  names(res)<-c("variable.type","mscore","msig")
  
  variable.list[[i]]<-res
  
}

for(i in 1:length(vars)){
  names(variable.list)[i]<-paste0(vars[[i]]$var,sep="","_",vars[[i]]$trt)
}

save.image(paste0(opt$rdataPath,sep="","/Lillysets.RData"))


## ................................................................STEP 4
## Get Models

#load("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/RData/Lillysets.RData")

source("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/scripts/Resources.R")

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
                       kfold = 10,repeats = 30,sel.feature = NULL)
      
      results.perm[[mdl]]<-res.i
      
    }else{
      
      res.i<-getML.num(training = training,testing = testing,
                       kfold = 10,repeats = 30,sel.feature = NULL)
      results.perm[[mdl]]<-res.i
      
    }
    
  }
  time2<-Sys.time()
  print(time2-time1)
  
  RESULTS[[var.i]]<-results.perm
  names(RESULTS)[var.i]<-names(variable.list)[var.i]
  
  #saveRDS(RESULTS,"MscoreModels.rds")
  
}

setwd("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/RData")

saveRDS(RESULTS,"Models_Lilly_mscore.rds")

set<-c(1,2,4,9,10)
for(i in 1:10){
  print(i)
  print(RESULTS$SRI5RESP_Q4W[[i]]$stats)
}


#for(i in 1:length(RESULTS$MAJRCRBG)){
#  print(RESULTS$MAJRCRBG[[i]]$stats)
#}

## 

MODELS.list<-list()
punt<-1
for(i in 1:length(RESULTS)){
  
  x<-prioML(lista = RESULTS[[i]],nameplot = NULL)
  x$model["AUC",]
  
  if(x$stats["AUC","stats"]>=0.7){
    MODELS.list[[punt]]<-x
    names(MODELS.list)[punt]<-names(RESULTS)[i]
    punt<-punt+1
  }
  
}

myPROSLE.model<-MODELS.list

setwd(opt$rdataPath)

saveRDS(myPROSLE.model,file = "myPROSLEmodels.RData")

save.image("LillyAnaisis.RData")



## 

#load("D:/DATA/WORK/Toro.et.al.MyPROSLE/RData/myPROSLEmodels.RData")

##modelStats(model.i = )

## AUC

x<-RESULTS$SRI5RESP_all[[9]]$models$nnet$trainingData

all<-variable.list$SRI5RESP_all$mscore

sel<-setdiff(rownames(all),rownames(x))
test<-all[sel,]

modelStats(model.i = RESULTS$SRI5RESP_all[[9]]$models$nnet,test = test)

tab.model<-RESULTS$SRI5RESP_all[[9]]$models$nnet
rm(list=setdiff(ls(),"tab.model"))


## ................................................................STEP 5
## Compare probs drug placebo

setwd("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022")

load("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/RData/TabalumabMscores.RData")

source("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/scripts/Resources.R")

## Primeras respuestas medidas...............
#clin.resp<-rbind(read.sas7bdat(paste0(getwd(),sep="","/Datasets/tabalumab/1/clinresp.sas7bdat")),
#         read.sas7bdat(paste0(getwd(),sep="","/Datasets/tabalumab/2/clinresp.sas7bdat")))
#clin.resp<-clin.resp[,c("USUBJID","TRT","MAJORCR","MAJRCRBG",
#        "MAJRCRSL","PARTLCR","PRTLCRBG","PRTLCRSL")]


## Medicion SRI......................
clin.resp<-rbind(read.sas7bdat(paste0(getwd(),sep="","/Datasets/tabalumab/1/sri.sas7bdat")),
                 read.sas7bdat(paste0(getwd(),sep="","/Datasets/tabalumab/2/sri.sas7bdat")))
clin.resp<-clin.resp[,c("USUBJID","TRT","VISID","SRI5RESP","SRI4RESP")]
# "SRI5RDUR","SRI5RRAW","SRI4RRAW"

patients<-unique(clin.resp$USUBJID)
Response<-as.data.frame(matrix(ncol=4,nrow=0))
colnames(Response)<-c("USUBJID","TRT","SRI5RESP","SRI4RESP")

for(i in 1:length(patients)){
  
  tmp<-clin.resp[clin.resp$USUBJID==patients[i],]
  if(16 %in% tmp$VISID){
    
    sel<-ifelse(tmp$VISID==16,T,F)
    Response<-rbind(Response,tmp[sel,c("USUBJID","TRT","SRI5RESP","SRI4RESP")])
    
  }
}

ids<-read.table(file=paste0(getwd(),sep="","/Datasets/tabalumab/1/Tabalumab_samples_clinical-usubjid.txt"),sep="\t",header = T)
# 1772 pacientes      


variable.list<-list() # "variable.type" "mscore" "msig" 

#vars<-c("MAJORCR","MAJRCRBG","MAJRCRSL","PARTLCR","PRTLCRBG","PRTLCRSL")
#vars<-c("SRI5RESP","SRI4RESP") #

vars<-list(list("var"="SRI5RESP","trt"="all"))


for(i in 1:length(vars)){
  
  # 699 patients treated (50:50) 
  # 327 placebo
  
  tmp<-Response[,c("USUBJID","TRT",vars[[i]]$var)]
  
  switch(vars[[i]]$trt,
         all={
           tmp<-tmp[ifelse(!is.na(tmp[,3]),T,F),]
         },
         Q4W={
           tmp<-tmp[ifelse(!is.na(tmp[,3]) & tmp$TRT=="LY 120 mg Q4W",T,F),]
         },
         Q2W={
           tmp<-tmp[ifelse(!is.na(tmp[,3]) & tmp$TRT=="LY 120 mg Q2W",T,F),]
         })
  
  # 
  tmp[,3]<-ifelse(tmp[,3]==0,"NO","YES") 
  ## NO: 389, YES: 309  SRI5
  ## NO: 279, YES: 419  SRI4
  
  ids.tmp<-ids[ids$Time=="baseline",]
  ids.tmp<-ids.tmp[ids.tmp$Clinical_USUBJID %in% tmp$USUBJID,]
  
  response<-NULL
  for(j in 1:nrow(ids.tmp)){
    response<-c(response,
                tmp[tmp$USUBJID==ids.tmp$Clinical_USUBJID[j],vars[[i]]$var])
    
  }
  ids.tmp$response<-response
  
  mscores<-SLE.mscores[,ids.tmp$GEO_accession]
  msig<-SLE.msig[,ids.tmp$GEO_accession]
  
  
  mscores<-as.data.frame(t(mscores))
  msig<-as.data.frame(t(msig))
  
  mscores<-cbind(ids.tmp$response,mscores)
  msig<-cbind(ids.tmp$response,msig)
  colnames(mscores)[1]<-"group"
  colnames(msig)[1]<-"group"
  
  variable.type<-"categoric"
  
  res<-list(variable.type,mscores,msig)
  names(res)<-c("variable.type","mscore","msig")
  
  variable.list[[i]]<-res
  
}

for(i in 1:length(vars)){
  names(variable.list)[i]<-paste0(vars[[i]]$var,sep="","_",vars[[i]]$trt)
}


res<-predict(tab.model,newdata = variable.list$SRI5RESP_all$mscore,type = "prob")

pats<-ids[rownames(res),"Clinical_USUBJID"]
resp=NULL
for(i in 1:length(pats)){
  
  tmp<-clin.resp[ifelse(clin.resp$USUBJID==pats[i],T,F),]
  
  if(tmp$TRT[1]=="Placebo"){
    resp<-c(resp,"Placebo")
  }else{
    resp<-c(resp,"Tabalumab")
  }
  
}


res<-cbind(resp,variable.list$SRI5RESP_all$mscore$group,res)
colnames(res)<-c("Treatment","SRI5","NO","YES")

library("ggplot2")

res<-cbind(paste0(res$Treatment,"_",sep="",res$SRI5),res)
colnames(res)[1]<-"group"

ggplot(res,aes(x=group,y=YES,fill=group))+geom_boxplot(outlier.colour = NA,alpha=0.9)+
  theme_classic()+scale_fill_manual(values=c("#7d8c93","#eebc4d","#7d8c93","#eebc4d"))

pn<-res[ifelse(res$group=="Placebo_NO",T,F),"YES"]
py<-res[ifelse(res$group=="Placebo_YES",T,F),"YES"]

tn<-res[ifelse(res$group=="Tabalumab_NO",T,F),"YES"]
ty<-res[ifelse(res$group=="Tabalumab_YES",T,F),"YES"]

min(wilcox.test(pn,py,alternative="greater")$p.value,
         wilcox.test(pn,py,alternative="less")$p.value)

min(wilcox.test(pn,tn,alternative="greater")$p.value,
    wilcox.test(pn,tn,alternative="less")$p.value)
# 8.95925e-23

min(wilcox.test(py,ty,alternative="greater")$p.value,
    wilcox.test(py,ty,alternative="less")$p.value)
# 3.47516e-21

min(wilcox.test(py,tn,alternative="greater")$p.value,
    wilcox.test(py,tn,alternative="less")$p.value)
# 2.273715e-13

min(wilcox.test(pn,ty,alternative="greater")$p.value,
    wilcox.test(pn,ty,alternative="less")$p.value)
# 9.243757e-22


#####################

library("caret")

imp<-varImp(tab.model)
imp<-data.frame("moduleID"=rownames(imp$importance),
                "importance"=imp$importance$Overall)
imp<-imp[order(imp$importance,decreasing = T),]
rownames(imp)<-imp$moduleID
imp<-cbind(imp,
           Modules.ann[rownames(imp),c("category","color")])

paths<-unique(imp$category)
means<-NULL
for(i in 1:length(paths)){
  
  tmp<-imp[imp$category==paths[i],]
  tmp<-mean(tmp$importance,na.rm = T)
  means<-c(means,tmp)
}

sum.imp<-data.frame("SLE_signature"=paths,
                    "importance"=means)

sum.imp<-sum.imp[order(sum.imp$importance,decreasing=T),]
sum.imp$SLE_signature<-factor(sum.imp$SLE_signature,levels=unique(sum.imp$SLE_signature))

ggplot(sum.imp,aes(x=SLE_signature,y=importance,fill="white"))+
  geom_bar(stat="identity",color="black")+theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

selimp<-rownames(imp)[1:10]

library("pheatmap")

M<-variable.list$SRI5RESP_all$mscore
M<-M[,c("group",selimp)]

ann<-data.frame(M[,1])
rownames(ann)<-rownames(M)

M2<-rbind(M[M$group=="YES",],
         M[M$group=="NO",])

pats<-ifelse(res$Treatment!="Placebo" & res$YES>0.9,T,F)
M2<-M2[pats,]

pheatmap(t(M2[,-1]),annotation_col = ann,
         cluster_cols = T,show_colnames = F,
         breaks = seq(-2,2,length.out = 100))



## ................................................................STEP 6
## Times

pat.list<-list()
pats<-ids.tmp #[ids.tmp$Diagnosis=="SLE",]
patsid<-unique(pats$SubjectID)
punt<-1

for(i in 1:length(patsid)){
  
  tmp<-clin.tab[clin.tab$SubjectID==patsid[i],]
  tmp<-tmp[tmp$Treatment!="None",]
  tmp$Treatment<-ifelse(tmp$Treatment=="Placebo","Placebo","Tabalumab")
  
  if("baseline" %in% tmp$Time & "week52" %in% tmp$Time & "week16" %in% tmp$Time){
    tmp<-cbind("samples"=rownames(tmp),tmp)
    
    tmp<-cbind(tmp,"response"=rep(pats[pats$SubjectID==patsid[i],"response"][1],nrow(tmp)))
    
    mtemp<-SLE.mscores[,tmp$samples]
    msig<-SLE.msig[,tmp$samples]
    
    tmp.pat<-list(tmp,mtemp,msig)
    names(tmp.pat)<-c("pat","mscore","msig")
    pat.list[[punt]]<-tmp.pat
    punt<-punt+1
  }

}



plot.list<-list()
for(var in 1:length(selimp)){
  meltM<-matrix(ncol=3,nrow=0)
  colnames(meltM)<-c("time","response","value")
  
  for(i in 1:length(pat.list)){
    
    if(pat.list[[i]]$pat$Treatment[2]!="Placebo"){
      
      value<-pat.list[[i]]$mscore[selimp[var],
                                  pat.list[[i]]$pat$samples]
      
      x<-data.frame("time"=pat.list[[i]]$pat$Time,
               "response"=pat.list[[i]]$pat$response,
               "value"=as.numeric(value))
      meltM<-rbind(meltM,x)
      
    }
    
  }
  values<-c(mean(meltM[meltM$time=="baseline" & meltM$response=="YES","value"],na.rm=T),
            mean(meltM[meltM$time!="baseline" & meltM$response=="YES","value"],na.rm=T),
            mean(meltM[meltM$time=="baseline" & meltM$response=="NO","value"],na.rm=T),
            mean(meltM[meltM$time!="baseline" & meltM$response=="NO","value"],na.rm=T))
  
  Mpl<-data.frame("time"=c(0,1,0,1),
                  "response"=c("YES","YES","NO","NO"),
                  "value"=values)
  
  p1<-ggplot(Mpl,aes(x=time,y=value,color=response,group=response))+geom_line(size=1.2)+xlim(-0.1,1.1)+
    theme_bw()+ylim(min(values)-0.25,max(values)+0.25)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_color_manual(values=c("NO"="#7d8c93","YES"="#eebc4d"))
  plot.list[[var]]<-p1
  
  
}
plot.list

#Modules.ann[selimp,]
library("cowplot")
library("gridExtra")

grid.arrange(arrangeGrob(plot.list[[1]],plot.list[[2]],plot.list[[3]],plot.list[[4]],plot.list[[5]],ncol=5),
                arrangeGrob(plot.list[[6]],plot.list[[7]],plot.list[[8]],plot.list[[9]],plot.list[[10]],ncol=5),
                nrow = 2)  

## por signature

selsig<-rownames(pat.list[[1]]$msig)

plot.list2<-list()
for(var in 1:length(selsig)){
  meltM<-matrix(ncol=3,nrow=0)
  colnames(meltM)<-c("time","response","value")
  
  for(i in 1:length(pat.list)){
    
    if(pat.list[[i]]$pat$Treatment[2]!="Placebo"){
      
      value<-pat.list[[i]]$msig[selsig[var],
                                  pat.list[[i]]$pat$samples]
      
      x<-data.frame("time"=pat.list[[i]]$pat$Time,
                    "response"=pat.list[[i]]$pat$response,
                    "value"=as.numeric(value))
      meltM<-rbind(meltM,x)
      
    }
    
  }
  values<-c(mean(meltM[meltM$time=="baseline" & meltM$response=="YES","value"],na.rm=T),
            mean(meltM[meltM$time=="week16" & meltM$response=="YES","value"],na.rm=T),
            mean(meltM[meltM$time=="week52" & meltM$response=="YES","value"],na.rm=T),
            mean(meltM[meltM$time=="baseline" & meltM$response=="NO","value"],na.rm=T),
            mean(meltM[meltM$time=="week16" & meltM$response=="NO","value"],na.rm=T),
            mean(meltM[meltM$time=="week52" & meltM$response=="NO","value"],na.rm=T))
  
  Mpl<-data.frame("time"=c(0,1,2,0,1,2),
                  "response"=c("YES","YES","YES","NO","NO","NO"),
                  "value"=values)
  
  p1<-ggplot(Mpl,aes(x=time,y=value,color=response,group=response))+geom_line(size=1.2)+xlim(-0.1,2.1)+geom_point(size=1.3)+
    theme_bw()+ylim(min(values)-0.25,max(values)+0.25)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_color_manual(values=c("NO"="#7d8c93","YES"="#eebc4d"))+ggtitle(selsig[var])
  plot.list2[[var]]<-p1
  
  
}
plot.list2

#Modules.ann[selimp,]
library("cowplot")
library("gridExtra")

grid.arrange(arrangeGrob(plot.list2[[1]],plot.list2[[2]],plot.list2[[3]],plot.list2[[4]],plot.list2[[5]],ncol=5),
             arrangeGrob(plot.list2[[6]],plot.list2[[7]],plot.list2[[8]],plot.list2[[9]],ncol=5),
             nrow = 2)  






## ................................................................STEP 7
## Pie chart

#variable.list<-list() # "variable.type" "mscore" "msig" 

#vars<-c("MAJORCR","MAJRCRBG","MAJRCRSL","PARTLCR","PRTLCRBG","PRTLCRSL")
#vars<-c("SRI5RESP","SRI4RESP") #


  # 699 patients treated (50:50) 
  # 327 placebo
  
  tmp<-Response[,c("USUBJID","TRT","SRI5RESP")]
  tmp$TRT<-ifelse(tmp$TRT=="Placebo","Placebo","Tabalumab")
  tmp$SRI5RESP<-ifelse(tmp$SRI5RESP==0,"NO","YES") 
  

  ids.tmp<-ids[ids$Time=="week52",] # baseline
  ids.tmp<-ids.tmp[ids.tmp$Clinical_USUBJID %in% tmp$USUBJID,]
  
  response<-NULL
  for(j in 1:nrow(ids.tmp)){
    response<-c(response,
                tmp[tmp$USUBJID==ids.tmp$Clinical_USUBJID[j],"SRI5RESP"])
    
  }
  ids.tmp$response<-response
  
  trt<-NULL
  for(p in 1:nrow(ids.tmp)){
    trt<-c(trt,tmp[tmp$USUBJID==ids.tmp$Clinical_USUBJID[p],"TRT"])
  }
  ids.tmp<-cbind(ids.tmp,trt)
  
  mscores<-SLE.mscores[,ids.tmp$GEO_accession]
  msig<-SLE.msig[,ids.tmp$GEO_accession]
  
  
  mscores<-as.data.frame(t(mscores))
  msig<-as.data.frame(t(msig))
  
  mscores<-cbind(paste0(ids.tmp$response,"_",ids.tmp$trt,sep=""),mscores)
  msig<-cbind(paste0(ids.tmp$response,"_",ids.tmp$trt,sep=""),msig)
  colnames(mscores)[1]<-"group"
  colnames(msig)[1]<-"group"
  
  
Events<-NULL
for(coln in 2:ncol(msig)){
  active<-NULL
  for(r in 1:nrow(msig)){
    
    if(msig[r,coln]>=1.65){
      active<-c(active,as.character(msig[r,"group"]))
    }
  }
  active<-active[active!="NA_Placebo"]
  
  print(colnames(msig)[coln])
  print(table(active))
  
  x<-table(active)
  
  sum(x[c("NO_Tabalumab","YES_Tabalumab")])
  
  print(as.numeric(x["YES_Placebo"])/sum(as.numeric(x[c("YES_Placebo","NO_Placebo")])))
  print(as.numeric(x["YES_Tabalumab"])/sum(as.numeric(x[c("YES_Tabalumab","NO_Tabalumab")])))
  
  m<-data.frame("Placebo"=c(as.numeric(x["NO_Placebo"]),
                as.numeric(x["YES_Placebo"])),
                "Tabalumab"=c(as.numeric(x["NO_Tabalumab"]),
                              as.numeric(x["YES_Tabalumab"])),
                "Type"=c("Nonresponse","Response"),
                "Signature"=rep(colnames(msig)[coln],2)
                )
Events<-rbind(Events,m)
  
}

##------
## Get significance (Fisher test / Hypergeometric test)

signatures<-unique(Events$Signature)
pvals<-NULL

for(i in 1:length(signatures)){
  
  x<-ifelse(Events$Signature==signatures[i],T,F)
  test<-t(as.data.frame(cbind(Events[x,"Placebo"],Events[x,"Tabalumab"])))
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


# baseline:
E.bas<-Events
E.wk<-Events
#23/(16+23); 31/(22+31); 45/50+45; 38/(38+41); 25/(25+23); 104/(104+108);237/(237+248); 34/(34+34);16/(16+11)
#22/43; 21/44; 30/62; 25/(25+32); 5/6; 40/80; 55/(55+47); 28/(19+28); 11/23 

E.bas$Tabalumab










  ###########################################################
  

    tmp<-Response[,c("USUBJID","TRT","SRI5RESP")]
    tmp<-tmp[ifelse(!is.na(tmp[,"SRI5RESP"]),T,F),]
          
    tmp[,"SRI5RESP"]<-ifelse(tmp[,"SRI5RESP"]==0,"NO","YES") 
    tmp$TRT<-ifelse(tmp$TRT=="Placebo","Placebo","Tabalumab")

    tmp<-tmp[ifelse(tmp$TRT!="Placebo",T,F),]
    
    clin.pat.resp<-tmp
    
    
    ids.tmp<-ids[ids$Time!="week16",]
    ids.tmp<-ids[ids$Clinical_USUBJID %in% clin.pat.resp$USUBJID,]
    
    response<-NULL
    for(j in 1:nrow(ids.tmp)){
  
      response<-c(response,
                  tmp[tmp$USUBJID==ids.tmp$Clinical_USUBJID[j],"SRI5RESP"])
      
    }
    ids.tmp$response<-response
    rownames(ids.tmp)<-ids.tmp$GEO_accession
    ids.tmp<-ids.tmp[ids.tmp$Diagnosis!="SLE",]
    metadata<-ids.tmp[,c("Clinical_USUBJID","response","Time")]
    
    
    mscores<-SLE.mscores[,rownames(metadata)]
    msig<-SLE.msig[,rownames(metadata)]
    
    
## Comparaciones
    
## Differential expression between two categorical classes
    
dolimma<-function(data,    ## Matrix(genes in rows, samples in columns)
                  cl,      ## Numeric vector (0: reference sample, 1: case sample)
                  adjust="bonferroni"){  ## Pvalue multiple correction method, see details  
    require("limma")
      
    design<-matrix(data=c(ifelse(cl==0,1,0),ifelse(cl==1,1,0)),
                    ncol=2,nrow=ncol(data),
                    dimnames = list(c(colnames(data)),c("CONTROL","CASE")))
    fit <- lmFit(data, design)
    contrast.matrix <- makeContrasts(CASE-CONTROL, levels=design)
    fit2 <- contrasts.fit(fit, contrast.matrix); 
    fit2 <- eBayes(fit2)
    res<- topTable(fit2, coef=1, adjust.method =adjust,number = nrow(data))
    res<-cbind(rownames(res),res); colnames(res)[1]<-"Genes"
    return(res)
  }
  

## Baseline, resp vs no resp

sel<-ifelse(metadata$Time=="baseline",T,F)
sel<-metadata[sel,]

data<-mscores[,rownames(sel)]
cl<-ifelse(sel$response=="NO",0,1)

res<-dolimma(data = data,
             cl = cl,
            adjust = "bonferroni")

res<-res[ifelse(res$P.Value<=0.05,T,F),]

res.base<-res


results1<-data[res.base$Genes,]
results1<-data.frame("NO"=mean(as.numeric(results1[cl==0])),
                     "YES"=mean(as.numeric(results1[cl==1])))
rownames(results1)<-res.base$Genes

##.........................
## resp baseline wk52
sel<-ifelse(metadata$Time!="week16" & metadata$response=="YES",T,F)
sel<-metadata[sel,]

data<-mscores[,rownames(sel)]
cl<-ifelse(sel$Time=="baseline",0,1)

res<-dolimma(data = data,
             cl = cl,
             adjust = "bonferroni")

res<-res[ifelse(res$P.Value<=0.05,T,F),]
#Modules.ann[rownames(res),]

res.response<-res

results2<-data[res.response$Genes,]
r1<-results2[,cl==0]
r2<-results2[,cl==1]

r1<-apply(r1,1,mean)
r2<-apply(r2,1,mean)

results2<-data.frame("NO"=r1,
                     "YES"=r2)


rownames(results2)<-res.response$Genes

results2<-data.frame("NO"=mean(as.numeric(results1[cl==0])),
                     "YES"=mean(as.numeric(results1[cl==1])))

##

##.........................
## noresp baseline wk52
sel<-ifelse(metadata$Time!="week16" & metadata$response=="NO",T,F)
sel<-metadata[sel,]

data<-mscores[,rownames(sel)]
cl<-ifelse(sel$Time=="baseline",0,1)

res<-dolimma(data = data,
             cl = cl,
             adjust = "bonferroni")

res<-res[ifelse(res$P.Value<=0.05,T,F),]
#Modules.ann[rownames(res),]

res.noresponse<-res

results3<-data[res.noresponse$Genes,]
r1<-results3[,cl==0]
r2<-results3[,cl==1]

r1<-apply(r1,1,mean)
r2<-apply(r2,1,mean)

results3<-data.frame("NO"=r1,
                     "YES"=r2)

rownames(results3)<-res.noresponse$Genes


#library("pheatmap")

#m0<-data[rownames(res.response),ifelse(cl==0,T,F)]
#p0<-pheatmap(t(m0),show_rownames = F,
#             cluster_cols = F,cluster_rows = T,show_colnames = T,
#             breaks = seq(-2,2,length.out = 100))
#m1<-data[rownames(res.response),ifelse(cl==1,T,F)]
#p1<-pheatmap(t(m1),show_rownames = F,
#             cluster_cols = F,cluster_rows = T,show_colnames = T,
#             breaks = seq(-2,2,length.out = 100))

#m0<-m0[,p0$tree_row$order]
#m1<-m1[,p1$tree_row$order]


#m<-cbind(m0[,c(ncol(m0):1)],m1)
#cl<-c(rep("baseline",80),rep("week52",80))

#ann<-data.frame("time"=cl)
#rownames(ann)<-colnames(m)

#pheatmap(t(m),show_rownames = F,
#         cluster_cols = F,cluster_rows = T,show_colnames = T,
#         breaks = seq(-3,3,length.out = 100),annotation_row = ann,
#         color = colorRampPalette(c("turquoise4","black","yellow"))(100))


##.........................
## resp no resp semana 52
sel<-ifelse(metadata$Time=="week52",T,F)
sel<-metadata[sel,]

data<-mscores[,rownames(sel)]
cl<-ifelse(sel$response=="NO",0,1)

res<-dolimma(data = data,
             cl = cl,
             adjust = "bonferroni")

res<-res[ifelse(res$P.Value<=0.05,T,F),]
Modules.ann[rownames(res),]

res.wk<-res

results4<-data[res.wk$Genes,]
r1<-results4[,cl==0]
r2<-results4[,cl==1]

r1<-apply(r1,1,mean)
r2<-apply(r2,1,mean)

results4<-data.frame("NO"=r1,
                     "YES"=r2)

rownames(results4)<-res.wk$Genes

#
ann1<-Modules.ann[rownames(results1),c("category","Function","color")]
ann2<-Modules.ann[rownames(results2),c("category","Function","color")]
ann3<-Modules.ann[rownames(results3),c("category","Function","color")]
ann4<-Modules.ann[rownames(results4),c("category","Function","color")]

rownames(ann1)<-paste("base_",sep="",ann1$Function)
rownames(ann2)<-paste("tres_",sep="",ann2$Function)
rownames(ann3)<-paste("tnres_",sep="",ann3$Function)
rownames(ann4)<-paste("wk52_",sep="",ann4$Function)

#

anns<-ann1[,"color"]; names(anns)<-rownames(ann1)
for(i in 1:nrow(ann2)){
  if(!rownames(ann2)[i] %in% names(anns)){
    x<-ann2[i,"color"]; names(x)<-rownames(ann2)[i]
    anns<-c(anns,x)
  }
}
for(i in 1:nrow(ann3)){
  if(!rownames(ann3)[i] %in% names(anns)){
    x<-ann3[i,"color"]; names(x)<-rownames(ann3)[i]
    anns<-c(anns,x)
  }
}
for(i in 1:nrow(ann4)){
  if(!rownames(ann4)[i] %in% names(anns)){
    x<-ann4[i,"color"]; names(x)<-rownames(ann4)[i]
    anns<-c(anns,x)
  }
}

## PLOT ##

rownames(results1)<-rownames(ann1)
rownames(results2)<-rownames(ann2)
rownames(results3)<-rownames(ann3)
rownames(results4)<-rownames(ann4)

gaps<-c(1,14,23)

anns.df<-data.frame("color"=anns)
rownames(anns.df)<-names(anns)

all<-rbind(results1,results2,results3,results4)

pheatmap(t(all),show_rownames = F,
         cluster_cols = F,cluster_rows = F,show_colnames = T,
         gaps_col = gaps,annotation_col = anns.df,
         breaks = seq(-2,2,length.out = 100),border_color = "black",
         color = colorRampPalette(c("turquoise4","black","yellow"))(100)
         )










