##############################
## MyPROSLE 
## R version R version 4.0.4
## 
##############################
## Time-dependent analysis

load("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/RData/Clustering.RData")
set.seed(123456788)
source(paste0(opt$scriptPath,"/Resources.R",sep=""))

load.libraries(set="06")

## Load clinical data from the two longitudinal datasets
load(paste0(opt$dataPath,sep="","/LongitudinalSets_metadata.RData"))

## ................................................................STEP 1
## Formatting metadata

## Select columns of interest
Clin.pascual<-Clin.pascual[,c("PatientID","Visit","Days_sinceLastvisit","Cumulative_time","SLEDAI")]

Clin.petri<-Clin.petri.visit[,c("CohortID","Date","SLEDAI")]
colnames(Clin.petri)[1]<-"PatientID"
rm(Clin.petri.visit)

## Split date information (Petri, dataset 9)
days<-NULL
month<-NULL
year<-NULL
for(i in 1:nrow(Clin.petri)){
  a<-unlist(strsplit(x=as.character(Clin.petri$Date[i]), split="/", 
                     fixed = FALSE, perl = FALSE, useBytes = FALSE))
  days<-c(days,as.numeric(as.character(a[1])))
  month<-c(month,as.numeric(as.character(a[2])))
  year<-c(year,as.numeric(as.character(a[3])))
}
Clin.petri<-cbind(Clin.petri,days,month,year)

## Set Date into a cuantitative measurement
time<-NULL
for(i in 1:nrow(Clin.petri)){
  if(is.na(Clin.petri$days[i])==FALSE){
    daysTime<-0
    daysTime<-daysTime + (Clin.petri$days[i]-1)
    m<-(Clin.petri$month[i]-1); m<-m*30
    daysTime<-daysTime+m
    ## Using 2005 as year of reference
    y<-(Clin.petri$year[i]-2005);y<-y*360 
    daysTime<-daysTime+y
    time<-c(time,daysTime)
  }else{
    time<-c(time,NA)
  }
}

Clin.petri<-cbind(Clin.petri,time)
Clin.petri<-Clin.petri[ifelse(is.na(Clin.petri$time)==FALSE,T,F),]

##-------
## Join metadata from both datasets
Clin.pascual<-Clin.pascual[,c("PatientID","Cumulative_time","SLEDAI")]
Clin.petri<-Clin.petri[,c("PatientID","time","SLEDAI")]
colnames(Clin.petri)[2]<-"Cumulative_time"

clin<-rbind(Clin.pascual,Clin.petri)

rm(Clin.pascual,Clin.petri,a,days,daysTime,i,m,month,time,y,year)

## Join Mscores for gene-modules from both datasets
Ms.time<-cbind(DATA.Mscore$dataset1$SLE,DATA.Mscore$dataset9$SLE)
Ms.time<-Ms.time[rownames(Modules.ann),rownames(clin)]

## Summarize gene-modules into SLE-signatures
DATA.Msig<-GetSignatures(listM = DATA.Mscore,ann = Modules.ann)

## Join Mscores for SLE-signatures from both datasets
Msig.time<-cbind(DATA.Msig$dataset1,DATA.Msig$dataset9)
Msig.time<-Msig.time[,rownames(clin)]


## ................................................................STEP 2
## Extract relevant information to time-dependent analysis (with disease state)

## Set maximal SLEDAI value to consider patient in clinical remission
max.remission.value<-2


REM<-data.frame(matrix(ncol=5,nrow=0))
## Patient, sampleID, cummulative time, Time from Last active state, Time to next active state
colnames(REM)<-c("Pat","Sample","cumTime","TfLS","TfNS") 

patients<-unique(clin$PatientID)

## Extract information patient by patient
for(pats in patients){
  
  tmp<-clin[ifelse(clin$PatientID==pats,T,F),]
  tmp<-tmp[order(as.numeric(as.character(tmp$Cumulative_time)),decreasing = F),]
  tmp<-tmp[ifelse(is.na(tmp$SLEDAI),F,T),]
  
  ## At least 3 time points
  if(nrow(tmp)>=3){
    
    ## Getting times from/to events
    for(t in 2:nrow(tmp)){ ## First point is discarded due to we dont have previous information
      
      if(as.numeric(as.character(tmp$SLEDAI[t]))<=max.remission.value){ ## Visit in remission
        
        ## Keep information of the patient in remission
        res<-pats
        res<-c(res,rownames(tmp)[t],
               as.numeric(tmp$Cumulative_time[t])-as.numeric(tmp$Cumulative_time[1]))
        
        ## Past
        x0<-"NO"
        t0<-t-1
        while(t0>=1){
          if(as.numeric(as.character(tmp$SLEDAI[t0]))>max.remission.value){ ## Flare
            x0<-abs(as.numeric(tmp$Cumulative_time[t])-
                      as.numeric(tmp$Cumulative_time[t0]))
            t0<-0
          }else{
            t0<-t0-1
          }## else
        }## while
        
        ## Next
        x1<-"NO"
        t0<-t+1
        while(t0<=nrow(tmp)){ ##while
          if(as.numeric(as.character(tmp$SLEDAI[t0]))>max.remission.value){  ## Flare
            x1<-abs(as.numeric(tmp$Cumulative_time[t])-
                      as.numeric(tmp$Cumulative_time[t0]))
            t0<-1000
          }else{
            t0<-t0+1
          } ## else
        } ## while
        
        res<-c(res,x0,x1)
        names(res)<-c("Pat","Sample","cumTime","TfLS","TfNS") 
        REM<-rbind(REM,res)
        
        
      } ## In remission
      
    } ## t
  }
}
colnames(REM)<-c("Pat","Sample","cumTime","TfLS","TfNS") 

rm(clin,tmp,max.remission.value,patients,pats,res,t,t0,x0,x1)


## ................................................................STEP 3
## Comparison: Time from LAST and NEXT active state

##-------
## Get available point for: LAST

## 1. Filter points without previous information
REM.last<-REM[ifelse(as.character(REM$TfLS)=="NO",F,T),]

## 2. Filter points with a near flare (100 days)
Type<-rep(TRUE,nrow(REM.last))
for(i in 1:nrow(REM.last)){
  if(as.character(REM.last$TfNS[i])!="NO"){
    if(as.numeric(as.character(REM.last$TfNS[i]))<100){
      Type[i]<-F
    }}
}
REM.last<-REM.last[Type,]


##-------
## Get available point for: NEXT

## 1. Filter samples without information for previous and/or following time point  
Type<-rep(FALSE,nrow(REM))
for(i in 1:nrow(REM)){
  if(as.character(REM$TfNS[i])!="NO" & as.character(REM$TfLS[i])!="NO"){
    Type[i]<-TRUE
  }
}
REM.next<-REM[Type,]

## 2. Filter samples with a recent flare (less than 100 days)
REM.next<-REM.next[ifelse(as.numeric(as.character(REM.next$TfLS))<100,F,T),]

rm(i,Type)


##-------
## Get significance for 9 main SLE-signatures

Surv.l1<-SurvAnalysis(rem.matrix = REM.last,
                     mscores.matrix = Msig.time,
                     thresh=1.2,
                     pval.th = 0.05,
                     ann=NULL,
                     dir="Absolute",
                     Type="Last")

Surv.l2<-SurvAnalysis(rem.matrix = REM.next,
                     mscores.matrix = Msig.time,
                     thresh=1.2,
                     pval.th = 0.05,
                     ann=NULL,
                     dir="Absolute",
                     Type="Next")

##-------
## Get survival plots

dev.off()
tiff(filename=paste0(opt$resultsPath,sep="","/Time_fromLastFlare_.tiff"),res = 300,width = 6,height = 5.4,units="in")
arrange_ggsurvplots(Surv.l1, print=TRUE,
                    ncol = 3, nrow = 2, risk.table.height = NULL)
dev.off()

tiff(filename=paste0(opt$resultsPath,sep="","/Time_toNextFlare_.tiff"),res = 300,width = 6,height = 2.7,units="in")
arrange_ggsurvplots(Surv.l2, print=TRUE,
                    ncol = 3, nrow = 1, risk.table.height = NULL)
dev.off()


##-------
## Get significance for gene-modules

Surv.mL<-SurvAnalysis(rem.matrix = REM.last,
                      mscores.matrix = Ms.time,
                      thresh=1.2,
                      pval.th = 0.05,
                      ann=Modules.ann,
                      dir="Absolute",
                      Type="Last")

Surv.mN<-SurvAnalysis(rem.matrix = REM.next,
                      mscores.matrix = Ms.time,
                      thresh=1.2,
                      pval.th = 0.05,
                      ann=Modules.ann,
                      dir="Absolute",
                      Type="Next")

Surv.mL$Pvalues<-(-log10(Surv.mL$Pvalues))
Surv.mN$Pvalues<-(-log10(Surv.mN$Pvalues))

surv.vals<-cbind(Surv.mL$Pvalues,Surv.mN)
colnames(surv.vals)[1:2]<-c("Pvalue_Last","Pvalue_Next")



#maxi<-apply(surv[,1:2],1,max)
#sel<-ifelse(maxi<1.30103,T,F)
#names<-rownames(surv)
#names[sel]<-" "

dev.off()

tiff(filename=paste0(opt$resultsPath,sep="","/Time_NextLast_GeneModules_.tiff"),res = 300,width = 3.5,height = 3.5,units="in")
p1<-ggplot(surv.vals,aes(x=Pvalue_Last,y=Pvalue_Next)) + theme_classic()+
  geom_point(color=surv.vals$color,size=2.5) +
  geom_vline(xintercept=1.30103,color="grey") +
  geom_hline(yintercept=1.30103,color="grey")
plot(p1)
dev.off()


rm(fit,p,p1,REM.last,REM.next,Surv.l1,Surv.l2,Surv.mL,Surv.mN,surv.vals,survObject)

## ................................................................STEP 4
## Comparison: V and L points

##-------
## Get available points
## Filter 1: Remove samples without previous/posterior information
Type<-rep(F,nrow(REM))
for(i in 1:nrow(REM)){
  if(as.character(REM$TfNS[i])!="NO" & as.character(REM$TfLS[i])!="NO"){
    Type[i]<-T
  }
}
REM.VL<-REM[Type,]

## Filter 2: Remove patients that do not become from a near flare
Type<-rep(T,nrow(REM.VL))
for(i in 1:nrow(REM.VL)){
  if(as.numeric(as.character(REM.VL$TfLS[i]))>90){ ## 90
    Type[i]<-F
  }
}
REM.VL<-REM.VL[Type,]

## Divide V and L patient
Type<-rep(T,nrow(REM.VL))
for(i in 1:nrow(REM.VL)){
  if(as.numeric(as.character(REM.VL$TfNS[i]))>90){ ## 90
    Type[i]<-F
  }
}
REM_V<-REM.VL[Type,]
REM_L<-REM.VL[!Type,]

patsV<-as.character(REM_V$Sample); 
patsL<-as.character(REM_L$Sample); 

cat(paste0("\nV Patients: ",length(patsV),"\nL Patients: ",length(patsL),sep=""))

## Object to perform comparisons
Group<-c(rep(1,length(patsV)),rep(2,length(patsL)))
VL.mscore<-Ms.time[,c(patsV,patsL)]
VL.msig<-Msig.time[,c(patsV,patsL)]


## Get significant differences
annot=Modules.ann[,"category"]; names(annot)<-rownames(Modules.ann)

GroupComparison(mdata=VL.mscore,group=Group,annot=annot,namePlot = "Time_VL_modules",
                Width = 5.8,Height = 5,pval.th = 0.01,method = "wilcox",
                Ylim=c(-4,3.5),textSize = 2.6,colori=c("cadetblue2","darkseagreen4"),
                groupName=c("V","L"),pointPval = -4)

GroupComparison(mdata=VL.msig,group=Group,annot=NULL,namePlot = "Time_VL_signatures",
                Width = 3.5,Height = 5,pval.th = 0.05,method = "wilcox",
                Ylim=c(-4,3.5),textSize = 2.6,colori=c("cadetblue2","darkseagreen4"),
                groupName=c("V","L"),pointPval = -4)

##-------
## Check changes in Cell types

load(paste0(opt$dataPath,sep="","/Cells_pascual.RData"))

cells<-cells[intersect(colnames(VL.mscore),rownames(cells)),]

Groupscells<-Group; names(Groupscells)<-colnames(VL.mscore)
Groupscells<-Groupscells[rownames(cells)]
Groupscells<-as.numeric(Groupscells)

cells<-t(cells)

GroupComparison(mdata=cells,group=Groupscells,annot=NULL,namePlot = "Time_VL_cells",
                Width = 3.5,Height = 5,pval.th = 0.05,method = "wilcox",
                Ylim=c(-0.03,0.6),textSize = 2.6,colori=c("cadetblue2","darkseagreen4"),
                groupName=c("V","L"),pointPval = -0.035)

rm(cells,annot,Group,Groupscells,i,patsV,patsL,Type,REM_L,REM.VL,REM_V)


## ................................................................STEP 5
## ## Model to predict disease activation in X days

## ............................
## 100 Days

## 1. Filter points without posterior information
REM.tmp<-REM[ifelse(as.character(REM$TfNS)=="NO",F,T),]

## Optional (to remove samples with posible remain signals from past flare)
#REM.tmp<-REM.tmp[ifelse(as.character(REM.tmp$TfLS)=="NO",F,T),]
#REM.tmp<-REM.tmp[ifelse(as.numeric(REM.tmp$TfLS)<100,F,T),]

## 2. Divide patients with near flare or not
Type<-rep(TRUE,nrow(REM.tmp))
for(i in 1:nrow(REM.tmp)){
    if(as.numeric(as.character(REM.tmp$TfNS[i]))<100){
      Type[i]<-F
    }
}
REM.WS<-REM.tmp[!Type,]
REM.nW<-REM.tmp[Type,]

rm(REM.tmp)

patsW<-as.character(REM.WS$Sample); 
patsnW<-as.character(REM.nW$Sample); 

cat(paste0("\nWS Patients: ",length(patsW),"\nnWS Patients: ",length(patsnW),sep=""))

## Object to perform comparisons
Group<-c(rep("nW",length(patsnW)),rep("W",length(patsW)))
WS.mscore<-as.data.frame(t(Ms.time[,c(patsnW,patsW)]))
WS.msig<-as.data.frame(t(Msig.time[,c(patsnW,patsW)]))

WS.mscore<-cbind(Group,WS.mscore)
WS.msig<-cbind(Group,WS.msig)

WS.mscore$Group<-ifelse(WS.mscore$Group=="nW","NO","YES")
WS.msig$Group<-ifelse(WS.msig$Group=="nW","NO","YES")
colnames(WS.mscore)[1]<-"group"
colnames(WS.msig)[1]<-"group"

variable.type<-"categoric"

res100<-list(variable.type,WS.mscore,WS.msig)
names(res100)<-c("variable.type","mscore","msig")

## Get models

#model.ms3.worsening<-GetModel(data=WS.mscore,method="rf",bootstrap=0.8,kfold=10,repeatedCv=20,
#         perm.bias=10,featureSets=c(1,2,3,4,5,6,7,8,9,10,15,20,25,30),seed=12345)
        
#model.sig3.worsening<-GetModel(data=WS.msig,method="rf",bootstrap=0.8,kfold=10,repeatedCv=20,
#                          perm.bias=10,featureSets=c(1,2,3,4,5,6,7,8,9),seed=12345)              
                  
    
#rm(list=setdiff(ls(),"model.ms.worsening"))
#save.image("model_ws.RData")
  

## ............................
## 200 Days

## 1. Filter points without posterior information
REM.tmp<-REM[ifelse(as.character(REM$TfNS)=="NO",F,T),]

## Optional (to remove samples with posible remain signals from past flare)
#REM.tmp<-REM.tmp[ifelse(as.character(REM.tmp$TfLS)=="NO",F,T),]
#REM.tmp<-REM.tmp[ifelse(as.numeric(REM.tmp$TfLS)<100,F,T),]

## 2. Divide patients with near flare or not
Type<-rep(TRUE,nrow(REM.tmp))
for(i in 1:nrow(REM.tmp)){
  if(as.numeric(as.character(REM.tmp$TfNS[i]))<200){
    Type[i]<-F
  }
}
REM.WS<-REM.tmp[!Type,]
REM.nW<-REM.tmp[Type,]

rm(REM.tmp)

patsW<-as.character(REM.WS$Sample); 
patsnW<-as.character(REM.nW$Sample); 

cat(paste0("\nWS Patients: ",length(patsW),"\nnWS Patients: ",length(patsnW),sep=""))

## Object to perform comparisons
Group<-c(rep("nW",length(patsnW)),rep("W",length(patsW)))
WS.mscore<-as.data.frame(t(Ms.time[,c(patsnW,patsW)]))
WS.msig<-as.data.frame(t(Msig.time[,c(patsnW,patsW)]))

WS.mscore<-cbind(Group,WS.mscore)
WS.msig<-cbind(Group,WS.msig)


WS.mscore$Group<-ifelse(WS.mscore$Group=="nW","NO","YES")
WS.msig$Group<-ifelse(WS.msig$Group=="nW","NO","YES")
colnames(WS.mscore)[1]<-"group"
colnames(WS.msig)[1]<-"group"

variable.type<-"categoric"

res200<-list(variable.type,WS.mscore,WS.msig)
names(res200)<-c("variable.type","mscore","msig")

variable.list<-list()
variable.list[[1]]<-res100
variable.list[[2]]<-res200
names(variable.list)<-c("worsening3","worsening6")

saveRDS(variable.list,"D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/RData/ObjectWS.rds")

## Get models

#model.ms6.worsening<-GetModel(data=WS.mscore,method="rf",bootstrap=0.8,kfold=10,repeatedCv=20,
#                             perm.bias=10,featureSets=c(1,2,3,4,5,6,7,8,9,10,15,20,25,30),seed=12345)

#model.sig6.worsening<-GetModel(data=WS.msig,method="rf",bootstrap=0.8,kfold=10,repeatedCv=20,
#                              perm.bias=10,featureSets=c(1,2,3,4,5,6,7,8,9),seed=12345)              

#setwd(opt$modelPath)
#rm(list=setdiff(ls(),c("model.ms3.worsening","model.ms6.worsening")))

#MODELS<-list()
#MODELS[[1]]<-model.ms3.worsening; names(MODELS)[1]<-"Worsening.3months"
#MODELS[[2]]<-model.ms6.worsening; names(MODELS)[2]<-"Worsening.6months"

#rm(list=setdiff(ls(),c("MODELS")))

#save.image("Models.RData")

