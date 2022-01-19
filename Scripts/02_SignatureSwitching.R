##############################
## MyPROSLE 
## R version R version 4.0.4
## 
##############################
## Prepare longitudinal data and shared Signatures between pathotypes 

## @@@@@@@@@@@@@@@@@@@@@@@@@@@ Input options
## Input general options

load("D:/DATA/WORK/Toro.et.al.MyPROSLE/RData/Clustering.RData")
set.seed(123456788)
source(paste0(opt$scriptPath,"/Resources.R",sep=""))

pkgs<-c("tmod","caret","matrixStats","biomaRt","NOISeq",
        "NbClust","SNFtool","tidyr","parallel","ggplot2",
        "ConsensusClusterPlus","raster","pheatmap")
check.packages(pkgs)
rm(pkgs)

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

annrow<-data.frame(Modules.ann$path10)
rownames(annrow)<-rownames(Modules.ann)
colnames(annrow)<-"Signature"

##-------
## Filtering non high dysregulated pathways
Sw.rate<-Sw.rate[ifelse(apply(Sw.rate,1,max)>=0.05,T,F),ifelse(apply(Sw.rate,1,max)>=0.05,T,F)]

for(i in 1:nrow(Sw.rate)){
  Sw.rate[i,i]<-1
}


dev.off()
tiff(filename=paste0(opt$resultsPath,sep="","/Switching_Figure1B__1.tiff"),res = 300,width = 4.9,height = 4.2,units="in")
heatmap<-pheatmap(Sw.rate,show_rownames = F,show_colnames = F, 
                  cluster_rows = T,cluster_cols = T,
                  breaks=seq(0,0.3,length.out = 100),
                  annotation_col = annrow,
                  annotation_row = annrow,
                  annotation_colors = Module.colors,
                  color = colorRampPalette(c("white","seagreen"))(100),
                  border_color = F)

dev.off()






