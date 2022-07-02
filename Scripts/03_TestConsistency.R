##############################
## MyPROSLE 
## R version R version 4.0.4
## 
##############################
## Test consistency

setwd("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/RData")

load("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/RData/Clustering.RData")

Modules.list<-Modules.list[rownames(Modules.ann)]

## Generate random subsets 

RES<-list()

for(d in 1:length(DATA)){
  print(d)
  
  
  ## Create random subsets
  H.i<-DATA[[d]]$Healthy
  subsets<-list()
  for(i in 1:20){
    subsets[[i]]<-sample(1:ncol(H.i),size = round(ncol(H.i)*0.9,digits = 0))
  }
  for(i in 21:40){
    subsets[[i]]<-sample(1:ncol(H.i),size = round(ncol(H.i)*0.8,digits = 0))
  }
  for(i in 41:60){
    subsets[[i]]<-sample(1:ncol(H.i),size = round(ncol(H.i)*0.7,digits = 0))
  }
  for(i in 61:80){
    subsets[[i]]<-sample(1:ncol(H.i),size = round(ncol(H.i)*0.6,digits = 0))
  }
  for(i in 81:100){
    subsets[[i]]<-sample(1:ncol(H.i),size = round(ncol(H.i)*0.5,digits = 0))
  }
  
  sle.subs<-sample(1:ncol(DATA[[d]]$SLE),size = round(ncol(DATA[[d]]$SLE)*0.5,digits = 0))
  
  results.d<-list()
  ## Calculate M-scores for each subset()
  for(sub in 1:length(subsets)){
    SLE.sub <- Get.MyPROSLE(modules.list = Modules.list, Patient=DATA[[d]]$SLE[,sle.subs],
                      Healthy = H.i[,subsets[[sub]]], method="mean")
    results.d[[sub]]<-SLE.sub
  }

  RES[[d]]<-results.d
  
}

save.image("RandomSubsets.RData")
# load("RandomSubsets.RData")

## Calcular la desviación estandar dentro de cada dataset (por rangos de 20 o todo)
## Objetivo, que la desciacion estandar sea pequeña, o sino buscar correlaciones, que sean altas


boots<-c(90,80,70,60,50)
range0<-c(1,21,41,61,81)
range1<-c(20,40,60,80,100)

meltRES<-matrix(ncol=2,nrow=0)
colnames(meltRES)<-c("boots","value")

for(b in 1:length(boots)){
  
  M.s<-NULL
  for(d in 1:length(RES)){ #dataset
    pats<-unique(as.character(colnames(RES[[d]][[1]])))
    
    for(p in 1:length(pats)){ #patient
      pat.sub<-matrix(data=0,ncol=0,nrow=206)
      
      for(i in range0[b]:range1[b]){ #subset
        x<-RES[[d]][[i]][,pats[p]]
        pat.sub<-cbind(pat.sub,x)
      } #subset
      
      pat.sub<-apply(pat.sub,1,sd)
      
      M.s<-c(M.s,pat.sub)
      
    } #patient
    
  }
  
  x<-cbind(rep(boots[b],length(M.s)),M.s)
  colnames(x)<-c("boots","value")
  meltRES<-rbind(meltRES,x)
  
}

library("ggplot2")

meltRES<-as.data.frame(meltRES)
meltRES$boots<-paste0(sep="",meltRES$boots,"%")
meltRES$boots<-factor(meltRES$boots,levels = unique(meltRES$boots))



ggplot(meltRES,aes(x=boots,y=value))+theme_bw()+
  geom_boxplot(outlier.fill = NA) + ylim(0,0.085)+
  xlab("random subset size")+ylab("standard deviation")

  

## Mas subsets ###############################33

d=1
  
  ## Create random subsets
  H.i<-DATA[[d]]$Healthy
  subsets<-list()
  for(i in 1:20){
    subsets[[i]]<-sample(1:ncol(H.i),size = round(ncol(H.i)*0.4,digits = 0))
  }
  for(i in 21:40){
    subsets[[i]]<-sample(1:ncol(H.i),size = round(ncol(H.i)*0.3,digits = 0))
  }
  for(i in 41:60){
    subsets[[i]]<-sample(1:ncol(H.i),size = round(ncol(H.i)*0.2,digits = 0))
  }
  for(i in 61:80){
    subsets[[i]]<-sample(1:ncol(H.i),size = round(ncol(H.i)*0.1,digits = 0))
  }
  
  sle.subs<-sample(1:ncol(DATA[[d]]$SLE),size = round(ncol(DATA[[d]]$SLE)*0.5,digits = 0))
  
  results.d<-list()
  ## Calculate M-scores for each subset()
  for(sub in 1:length(subsets)){
    SLE.sub <- Get.MyPROSLE(modules.list = Modules.list, Patient=DATA[[d]]$SLE[,sle.subs],
                            Healthy = H.i[,subsets[[sub]]], method="mean")
    results.d[[sub]]<-SLE.sub
  }
  
  RES2<-results.d
  
  boots<-c(40,30,20,10)
  range0<-c(1,21,41,61)
  range1<-c(20,40,60,80)
  
  meltRES2<-matrix(ncol=2,nrow=0)
  colnames(meltRES2)<-c("boots","value")
  
  for(b in 1:length(boots)){
    
    M.s<-NULL
      pats<-unique(as.character(colnames(RES2[[1]])))
      
      for(p in 1:length(pats)){ #patient
        pat.sub<-matrix(data=0,ncol=0,nrow=206)
        
        for(i in range0[b]:range1[b]){ #subset
          x<-RES2[[i]][,pats[p]]
          pat.sub<-cbind(pat.sub,x)
        } #subset
        
        pat.sub<-apply(pat.sub,1,sd)
        
        M.s<-c(M.s,pat.sub)
        
      } #patient
      

    
    x<-cbind(rep(boots[b],length(M.s)),M.s)
    colnames(x)<-c("boots","value")
    meltRES2<-rbind(meltRES2,x)
    
  }

##
  meltRES2<-as.data.frame(meltRES2)
  meltRES2$boots<-paste0(sep="",meltRES2$boots,"%")
  meltRES2$boots<-factor(meltRES2$boots,levels = unique(meltRES2$boots))
  
  MELT<-rbind(meltRES,meltRES2)
  
  ## POR AQUI VOY
  
  setwd("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/RData")
  
  saveRDS(MELT,"consistency.rds")

  setwd("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/Results")
  
tiff(filename="Consistency.tiff",res=300,width=3,height = 2.2, unit="in")
  p1<-ggplot(MELT,aes(x=boots,y=value))+theme_classic()+
    geom_boxplot(outlier.fill = NA,outlier.alpha = 0) + ylim(0,0.35)+
    xlab("Random subset size")+ylab("Standard deviation")+
    theme(axis.text.x = element_text(color = "grey20", size = 7, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 7,  face = "plain"),  
          axis.title.x = element_text(color = "grey20", size = 9, face = "bold"),
          axis.title.y = element_text(color = "grey20", size = 9,  face = "bold"))
    
    
  plot(p1)
dev.off()


  
  
  
  
  
  