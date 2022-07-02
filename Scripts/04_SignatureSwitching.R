##############################
## MyPROSLE 
## R version R version 4.0.4
## 
##############################
## Prepare longitudinal data and shared Signatures between pathotypes 

## @@@@@@@@@@@@@@@@@@@@@@@@@@@ Input options
## Input general options

load("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/RData/Clustering.RData")
set.seed(123456788)
source(paste0(opt$scriptPath,"/Resources.R",sep=""))

load.libraries(set="04")

## Test normality

values<-NULL
for(i in 1:nrow(Mscore.matrix)){
  values<-c(values,as.numeric(Mscore.matrix[i,]))
}
plot(density(values),xlim=c(-5,5))

## ................................................................STEP 1
## Extract longitudinal datasets

##-------
## Pascual
clin1<-read.csv(file=paste0(opt$dataPath,sep="","/Pascual_allClin.csv"),sep="\t")
rownames(clin1)<-clin1[,1];  clin1<-clin1[,-1]
clin1<-as.data.frame(clin1); clin1<-t(clin1); clin1<-as.data.frame(clin1)
clin1<-clin1[,c("Patient","Visit")]
clin1<-clin1[colnames(DATA.Mscore$dataset1$SLE),c("Patient","Visit")]

## Petri
clin2<-read.csv(file=paste0(opt$dataPath,sep="","/Metadata.petri.csv"),sep=";")
rownames(clin2)<-clin2$GZ_Filenames
clin2<-clin2[colnames(DATA.Mscore$dataset9$SLE),c("Subject_ID","Visit_Date")]
colnames(clin2)<-c("Patient","Visit")

clin.longitudinal<-rbind(clin1,clin2)
data.longitudinal<-cbind(DATA.Mscore$dataset1$SLE,DATA.Mscore$dataset9$SLE)

rm(clin1,clin2)

pats<-unique(clin.longitudinal$Patient)
cat(paste0("\nNumber of patients: ",length(pats),sep=""))


## ................................................................STEP 2
## Get switching rate

##-------
## Switching rate
Sw.rate<-matrix(data=0,ncol=nrow(Modules.ann),nrow=nrow(Modules.ann))
colnames(Sw.rate)<-rownames(Modules.ann)
rownames(Sw.rate)<-rownames(Modules.ann)
npats<-0
for(i in 1:length(pats)){
  x<-rownames(clin.longitudinal[clin.longitudinal$Patient==pats[i],])
  
  if(length(x)>=3){ ## At least 3 time points
    npats<-npats+1
    tmp<-data.longitudinal[,x]
    sel<-ifelse(apply(tmp,1,max)>=1.65,T,F)
    
    if(!is.na(table(sel)["TRUE"])){
      tmp<-rownames(tmp[sel,])
      
      for(pair.i in 1:length(tmp)){
        for(pair.j in 1:length(tmp)){
          if(pair.i!=pair.j){
            Sw.rate[tmp[pair.i],tmp[pair.j]]<-Sw.rate[tmp[pair.i],tmp[pair.j]]+1
            #Sw.rate[tmp[pair.j],tmp[pair.i]]<-Sw.rate[tmp[pair.j],tmp[pair.i]]+1
          }
        }
      }
    }## all false
  }
}

Sw.rate<-Sw.rate/npats



##-------
## Filtering non high dysregulated pathways
Sw.rate<-Sw.rate[ifelse(apply(Sw.rate,1,max)>=0.05,T,F),ifelse(apply(Sw.rate,1,max)>=0.05,T,F)]

for(i in 1:nrow(Sw.rate)){
  Sw.rate[i,i]<-1
}


annrow<-data.frame("Signature"=Modules.ann[colnames(Sw.rate),"path10"])
rownames(annrow)<-colnames(Sw.rate)


dev.off()
tiff(filename=paste0(opt$resultsPath,sep="","/Switching_Figure1B__1.tiff"),res = 300,width = 4.2,height = 3.45,units="in")
heatmap<-pheatmap(Sw.rate,show_rownames = F,show_colnames = F, 
                  cluster_rows = T,cluster_cols = T,
                  breaks=seq(0,0.3,length.out = 100),
                  annotation_col = annrow,
                  annotation_row = annrow,
                  annotation_colors = Module.colors,
                  color = colorRampPalette(c("white","seagreen"))(100),
                  border_color = F)

dev.off()

saveRDS(Sw.rate,paste0(opt$rdataPath,sep="","Heatmapswitch.rds"))

## ................................................................STEP 3
## Significance between different signature dysregulation across time

Msig<-GetSignatures(listM = data.longitudinal,ann = Modules.ann)

pats<-unique(clin.longitudinal$Patient)

Sw.rate<-matrix(data=0,ncol=nrow(Msig),nrow=nrow(Msig))
colnames(Sw.rate)<-rownames(Msig)
rownames(Sw.rate)<-rownames(Msig)
npats<-0

for(i in 1:length(pats)){
  x<-rownames(clin.longitudinal[clin.longitudinal$Patient==pats[i],])
  
  if(length(x)>=3){ ## At least 3 time points
    npats<-npats+1
    tmp<-Msig[,x]
    sel<-ifelse(apply(tmp,1,max)>=1.65,T,F)
    
    if(!is.na(table(sel)["TRUE"])){
      tmp<-rownames(tmp[sel,])
      
      for(pair.i in 1:length(tmp)){
        for(pair.j in 1:length(tmp)){
          if(pair.i!=pair.j){
            Sw.rate[tmp[pair.i],tmp[pair.j]]<-Sw.rate[tmp[pair.i],tmp[pair.j]]+1
            #Sw.rate[tmp[pair.j],tmp[pair.i]]<-Sw.rate[tmp[pair.j],tmp[pair.i]]+1
          }
        }
      }
    }## all false
  }
}

Nactive<-rep(0,nrow(Msig))
names(Nactive)<-rownames(Msig)

  for(i in 1:length(pats)){
   
    x<-rownames(clin.longitudinal[clin.longitudinal$Patient==pats[i],])
    
    if(length(x)>=3){
    
      tmp<-Msig[,x]
      tmp<-apply(tmp,1,max)
      tmp<-ifelse(tmp>=1.65,1,0)
      Nactive<-Nactive+tmp
    }
  }

##

pvals<-matrix(data=1,ncol=ncol(Sw.rate),nrow=nrow(Sw.rate))
colnames(pvals)<-colnames(Sw.rate)
rownames(pvals)<-rownames(Sw.rate)

for(i in 1:length(Nactive)){
  
  for(j in 1:ncol(Sw.rate)){
    if(i!=j){
     
      if(as.numeric(Nactive[i])>0){
        x<-prop.test(x=Sw.rate[i,j],n = as.numeric(Nactive[i])) 
        pvals[i,j]<-x$p.value
      }
    }
  }
}

sel<-ifelse(apply(pvals,1,min)>=0.05,F,T)

pvals<-pvals[sel,sel]

Sw.rate<-Sw.rate[sel,sel]

pvaL<-(-(log10(pvals)))

Nactive<-Nactive[sel]

##

Nactive
Sw.rate
pvaL

RESULTS<-data.frame(matrix(data=0,ncol=3,nrow=12))
colnames(RESULTS)<-c("Signatures","With","Pvalue")

RESULTS[1:3,1]<-"Plasma_cell_Cell_cycle"
RESULTS[4:6,1]<-"Interferon"
RESULTS[7:9,1]<-"Neutrophil_Inflammation"
RESULTS[10:12,1]<-"Platelet"

RESULTS[,2]<-c("Interferon","Neutrophil_Inflammation","Platelet",
               "Plasma_cell_Cell_cycle","Neutrophil_Inflammation","Platelet",
               "Plasma_cell_Cell_cycle","Interferon","Platelet",
               "Plasma_cell_Cell_cycle","Interferon","Neutrophil_Inflammation")

RESULTS[,3]<-c(8.875,-2.021,-4.509,
               -7.198,-1.469,-6.450,
               -5.654,20.395,-5.763,
               -3.14,7.15,1.00)

RESULTS$variable<-paste0(RESULTS$Signatures,sep="_",RESULTS$With)

RESULTS$colors<-c("#a9182d","#f29a42","#b486cb",
                  "#009000","#f29a42","#b486cb",
                  "#009000","#a9182d","#b486cb",
                  "#009000","#a9182d","#f29a42")

colors<-list(variable=c(Plasma_cell_Cell_cycle_Interferon="#a9182d",
                        Plasma_cell_Cell_cycle_Neutrophil_Inflammation="#f29a42",
                        Plasma_cell_Cell_cycle_Platelet="#b486cb",
                        Interferon_Plasma_cell_Cell_cycle="#009000",
                        Interferon_Neutrophil_Inflammation="#f29a42",
                        Interferon_Platelet="#b486cb",
                        Neutrophil_Inflammation_Plasma_cell_Cell_cycle="#009000",
                        Neutrophil_Inflammation_Interferon="#a9182d",
                        Neutrophil_Inflammation_Platelet="#b486cb",
                        Platelet_Plasma_cell_Cell_cycle="#009000",
                        Platelet_Interferon="#a9182d",
                        Platelet_Neutrophil_Inflammation="#f29a42"))

colors<-colors[[1]]

#library("ggbreak")

RESULTS$variable<-as.factor(rev(RESULTS$variable))

dev.off()
tiff(filename=paste0("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/Results",sep="","/Pvalues1D.tiff"),res = 300,width = 8,height = 2.8,units="in")

ggplot(RESULTS,aes(x=variable,y=Pvalue,fill=variable))+
  geom_bar(stat = "identity",color="black",width = 0.8)+
  coord_flip()+theme_classic()+geom_hline(yintercept = 0,size=1.1) +
  scale_fill_manual(values=colors)+ scale_y_break(c(9.5, 19.5))
 #+
#geom_hline(yintercept = 1.30,size=0.8,linetype="dashed",alpha=0.5)+
#  geom_hline(yintercept = -1.30,size=0.8,linetype="dashed",alpha=0.5)

dev.off()

save.image(paste0(opt$rdataPath,sep="","/switchPlot.RData"))






