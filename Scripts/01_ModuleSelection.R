##############################
## MyPROSLE 
## R version R version 4.0.4
## 
##############################
## Module selecction and Clustering

## @@@@@@@@@@@@@@@@@@@@@@@@@@@ Input options
## Input general options

load("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/RData/Mscores.RData")
set.seed(123456788)
source(paste0(opt$scriptPath,"/Resources.R",sep=""))

load.libraries(set="01")

opt<-AddOpt(option.list = opt,name.option = "filterVars",variable.value = c(10,3))
## variable.value = c("Percentage of patients with high dysregulated path in each dataset",
##                     "Number of datasets where "Percentaje" is satisfied)




## ................................................................STEP 1
## Select SLE-important paths

## Filter 1: High dysregulated modules in at least 10% of the samples and in more than 3 datasets

HighDys.perc<-data.frame(matrix(ncol=0,nrow=nrow(DATA.Mscore$dataset1$SLE)))

for(i in 1:length(DATA.Mscore)){
  sle<-apply(DATA.Mscore[[i]]$SLE,1,function(x){
    values<-as.numeric(x[is.na(x)==F])
    values<-(length(values[abs(values)>=1.65])/length(x))*100})
  HighDys.perc<-cbind(HighDys.perc,sle)
}
colnames(HighDys.perc)<-paste0("dataset",sep="",1:ncol(HighDys.perc))

selected.path<-apply(HighDys.perc,1,function(x){
  nperc<-length(x[x>opt$filterVars[1]])
  ntimes<-ifelse(nperc>opt$filterVars[2],T,F)
  return(ntimes)
})
cat(paste0("\nNumber of selected gene-modules: ",sep="",table(selected.path)["TRUE"]))

Modules.ann<-Modules.ann[selected.path,]
Modules.list<-Modules.list[selected.path]

rm(HighDys.perc,sle,selected.path,i)

#-------
## Filter path from DATA.Mscore
for(i in 1:length(DATA.Mscore)){
  DATA.Mscore[[i]]$SLE<-DATA.Mscore[[i]]$SLE[rownames(Modules.ann),]
  DATA.Mscore[[i]]$H<-DATA.Mscore[[i]]$H[rownames(Modules.ann),]
}


#-------
## Filter 2: Remove wrong measure paths (0 values = NA)
tmp<-data.frame(matrix(ncol=0,nrow=nrow(DATA.Mscore$dataset1$SLE)))
rownames(tmp)<-rownames(DATA.Mscore$dataset1$SLE)
for(i in 1:length(DATA.Mscore)){
  tmp<-cbind(tmp,DATA.Mscore[[i]]$SLE)  
}
tmp[tmp==0]<-NA
tmp<-t(apply(tmp,1,is.na))
tmp<-apply(tmp,1,table)

wrong.path<-NULL
for(i in 1:length(tmp)){
  if(!is.na(tmp[[i]]["TRUE"])){
    wrong.path<-c(wrong.path,names(tmp)[i])
  }
}

Modules.list<-Modules.list[setdiff(rownames(Modules.ann),wrong.path)]
Modules.ann<-Modules.ann[setdiff(rownames(Modules.ann),wrong.path),]

#-------
## Filter path from DATA.Mscore
for(i in 1:length(DATA.Mscore)){
  DATA.Mscore[[i]]$SLE<-DATA.Mscore[[i]]$SLE[rownames(Modules.ann),]
  DATA.Mscore[[i]]$H<-DATA.Mscore[[i]]$H[rownames(Modules.ann),]
}

rm(tmp,i,wrong.path)

## ................................................................STEP 2
## Clustering stability 
##Many of the modules are related and act on similar pathways

save.image(paste0(opt$rdataPath,sep="","/For_Random_stability.RData"))

#data<-DATA.Mscore$dataset1$SLE
#for(i in 2:length(DATA.Mscore)){
#  data<-cbind(data,DATA.Mscore[[i]]$SLE)
#}

# Elbow method
#fviz_nbclust(data, kmeans, method = "wss") +
#  geom_vline(xintercept = 4, linetype = 2)+
#  labs(subtitle = "Elbow method")

# Silhouette method
#fviz_nbclust(data, kmeans, method = "silhouette")+
#  labs(subtitle = "Silhouette method")

#fviz_nbclust(data, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
#  labs(subtitle = "Gap statistic method")

#-------
## Get distances between modules in each dataset
Dist.list<-list()
for(i in 1:length(DATA.Mscore)){
  data<-DATA.Mscore[[i]]$SLE
  Dist = (dist2(as.matrix(data),as.matrix(data)))^(1/2)
  Dist.list[[i]]<-Dist
}
rm(data,Dist,i)

#cl<-makeCluster(detectCores())
#Dist.list<-parLapply(cl, DATA.Mscore,int=1, GetDist)
#stopCluster(cl)
#for(i in 1:length(Dist.list)){
#  rownames(Dist.list[[i]])<-rownames(DATA.Mscore[[i]]$SLE)
#  colnames(Dist.list[[i]])<-rownames(DATA.Mscore[[i]]$SLE)
#}

#-------
## Set clustering parameters to permutate
K<-seq(10,30,by=1) ## number of neighbors, usually (10~30)
alpha<-seq(0.3,0.8,by=0.1) # alpha, hyperparameter, usually (0.3~0.8)

## 500 permutations selecting different number of datasets
permutations<-sample(2:(length(Dist.list)-1),500,replace=T)
permutations<-c(permutations,rep(length(Dist.list),2)) # To select also all datasets

print(table(permutations))
permutations<-as.list(permutations)

## For each permutation, all parameter combinations must be tested
combinations<-as.matrix(crossing(K,alpha)) ## 126 parameter combinations

#-------
## Save RData object
save.image(paste0(opt$rdataPath,sep="","/Clustering.RData"))

## load("D:/DATA/WORK/Toro.et.al.MyPROSLE/RData/Clustering.RData")

cl<-makeCluster(detectCores()-1)

cat("\nMeasuring cluster stability, this process is slow...")
stability<-parLapply(cl,permutations,
                     Dist.list=Dist.list,
                     combinations=combinations,
                     ClusterStats)
melt.stability = as.data.frame(do.call("rbind", stability))

stopCluster(cl)
rm(combinations,permutations,K,alpha,Dist.list,cl)

#-------
## Save RData object
save.image(paste0(opt$rdataPath,sep="","/Clustering.RData"))

## load("D:/DATA/WORK/Toro.et.al.MyPROSLE/RData/Clustering.RData")

## Get best number of clusters
nbplot<-aggregate(melt.stability$nbindex,by=list(datasets=melt.stability$datasets,clusters=melt.stability$clusters),data=melt.stability,FUN=mean)
colnames(nbplot)[3]<-"frecuency"

dev.off()
tiff(filename=paste0(opt$resultsPath,sep="","/Path_clustering.tiff"),res = 300,width = 6.5,height = 3.2,units="in")
p1<-ggplot(data=nbplot, aes(x=clusters, y=frecuency, group=datasets,color=as.factor(datasets))) +
  geom_line(size=1.1)+xlim(c(2,14))+
  geom_point(size=2) + theme_bw()
plot(p1)
dev.off()

nclust.plot<-nbplot

saveRDS(nclust.plot,paste0(opt$rdataPath,sep="","/nclustPlot.rds"))

cat("\nSelect nuber of clusters")
## 2, 4, 6, ( 9 )

rm(nbplot,stability,p1)


## ................................................................STEP 3
## Clustering assignation: Get clusters of paths and pats

#-------
## Create Joint-matrix
Mscore.matrix<-data.frame(matrix(ncol=0,nrow=nrow(DATA.Mscore$dataset1$SLE)))
rownames(Mscore.matrix)<-rownames(DATA.Mscore$dataset1$SLE)
for(i in 1:length(DATA.Mscore)){
  Mscore.matrix<-cbind(Mscore.matrix,DATA.Mscore[[i]]$SLE)  
}

#-------
## Consensus clustering
if (!file.exists(paste0(opt$resultsPath,sep="","/ConsensusCluster"))){
  dir.create(paste0(opt$resultsPath,sep="","/ConsensusCluster"))
} 

x<-t(Mscore.matrix)
d = sweep(x,1, apply(x,1,median,na.rm=T))
title=paste0(opt$resultsPath,sep="","/ConsensusCluster")

results = ConsensusClusterPlus(as.matrix(d),maxK=12,reps=300,
                               pItem=0.8,pFeature=1,
                               title=title,
                               clusterAlg="km",distance="euclidean",
                               innerLinkage = "complete",
                               seed=1262118388.71279,plot="pdf")

rm(x,d,title)
path10<-results[[10]]$consensusClass ## Un cluster es muy pequeño,se una al similar

#-------
## Save RData object
save.image(paste0(opt$rdataPath,sep="","/Clustering.RData"))


# platelet activation (III)
# platelet activation and degranulation
Modules.ann<-cbind(path10,Modules.ann)
Modules.ann[order(Modules.ann$path10,decreasing = F),]

Modules.ann[Modules.ann$path10==8,"path10"]<-7
Modules.ann[Modules.ann$path10==9,"path10"]<-8
Modules.ann[Modules.ann$path10==10,"path10"]<-9
Modules.ann<-Modules.ann[order(Modules.ann$path10,decreasing = F),]

Modules.ann<-cbind(Modules.ann,
results[[4]]$consensusClass[rownames(Modules.ann)],
results[[6]]$consensusClass[rownames(Modules.ann)])

colnames(Modules.ann)[4]<-"path4"
colnames(Modules.ann)[5]<-"path6"

Modules.ann<-Modules.ann[,c("path10","path6","path4",
                            "ID","Function")]


#-------
## Signature characterization

category<-NULL
color<-NULL
for(i in 1:nrow(Modules.ann)){
  
  AA = paste0("p",sep="",as.character(Modules.ann$path10[i]))
  switch(AA, 
         p1={
           category<-c(category,"Plasma cell/ Cell cycle")
           #color<-c(color,"#a8d8ba") 
           color<-c(color,"#009000")
         },
         p2={
           category<-c(category,"Neutrophil/ Inflammation")
           #color<-c(color,"#fdd7ab")   
           color<-c(color,"#f29a42")
         },
         p3={
           category<-c(category,"T cell")
           #color<-c(color,"#afd8f4")      
           color<-c(color,"#0099b2") 
         },
         p4={
           category<-c(category,"Nk cell")
           #color<-c(color,"#aebde1")     
           color<-c(color,"#43639d")
         },
         p5={
           category<-c(category,"Interferon")
           #color<-c(color,"#fba191")      
           color<-c(color,"#a9182d")
         },
         p6={
           category<-c(category,"Mitochondrion")
           #color<-c(color,"#b8b8b8")     
           color<-c(color,"#b0b0b0") 
         },
         p7={
           category<-c(category,"Platelet")
           #color<-c(color,"#cfb0d4")   
           color<-c(color,"#b486cb")
         },
         p8={
           category<-c(category,"B cell/ Plasma cell")
           #color<-c(color,"#cbd8a2")  
           color<-c(color,"#4fbb4b")
         },
         p9={
           category<-c(category,"Inositol metabolism")
           #color<-c(color,"#ebeeb2")  
           color<-c(color,"#86ffd7") 
         })
}

Modules.ann<-cbind(category,Modules.ann)
Modules.ann<-cbind(Modules.ann,color)

#-------
## Clustering of patients
Nb<-NbClust(data = Mscore.matrix, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 15,
            method = "complete", index = "ch", alphaBeale = 0.1)

plot(names(Nb$All.index),Nb$All.index)
## 5 and 11 clusters

if (!file.exists(paste0(opt$resultsPath,sep="","/ConsensusClusterPats"))){
  dir.create(paste0(opt$resultsPath,sep="","/ConsensusClusterPats"))
} 

x<-Mscore.matrix
d = sweep(x,1, apply(x,1,median,na.rm=T))
title=paste0(opt$resultsPath,sep="","/ConsensusClusterPats")

resultsPats = ConsensusClusterPlus(as.matrix(d),maxK=12,reps=300,
                               pItem=0.8,pFeature=1,
                               title=title,
                               clusterAlg="km",distance="euclidean",
                               innerLinkage = "complete",
                               seed=1262118388.71279,plot="pdf")
#-------
## Figure 1A (Heatmap + Frecuency)
ann.pats<-cbind(resultsPats[[5]]$consensusClass,
                resultsPats[[11]]$consensusClass)
colnames(ann.pats)<-c("K5","K11")
ann.pats<-as.data.frame(ann.pats)
ann.pats<-ann.pats[resultsPats[[11]]$consensusTree$order,]

Modules.ann$path10<-paste0("pats",Modules.ann$path10,sep="")


## Set new order (based on consensus)
Modules.ann<-rbind(Modules.ann[Modules.ann$path10=="pats4",],
         Modules.ann[Modules.ann$path10=="pats3",],
         Modules.ann[Modules.ann$path10=="pats6",],
         Modules.ann[Modules.ann$path10=="pats9",],
         Modules.ann[Modules.ann$path10=="pats8",],
         Modules.ann[Modules.ann$path10=="pats1",],
         Modules.ann[Modules.ann$path10=="pats5",],
         Modules.ann[Modules.ann$path10=="pats2",],
         Modules.ann[Modules.ann$path10=="pats7",])


ann.pats<-rbind(ann.pats[ann.pats$K11==5,],
                ann.pats[ann.pats$K11==11,],
                ann.pats[ann.pats$K11==4,],
                ann.pats[ann.pats$K11==7,],
                ann.pats[ann.pats$K11==1,],
                ann.pats[ann.pats$K11==2,],
                ann.pats[ann.pats$K11==3,],
                ann.pats[ann.pats$K11==10,],
                ann.pats[ann.pats$K11==9,],
                ann.pats[ann.pats$K11==6,],
                ann.pats[ann.pats$K11==8,])

rowGaps<-GetGap(vect=Modules.ann$path10)
colGaps<-GetGap(vect=ann.pats$K11)

dataset<-NULL
for(i in 1:nrow(ann.pats)){
  for(l in 1:length(DATA.Mscore)){
    if(rownames(ann.pats)[i] %in% colnames(DATA.Mscore[[l]]$SLE)){
      dataset<-c(dataset,names(DATA.Mscore[l]))
      break
    }
  }
}
ann.pats<-cbind(ann.pats,dataset)


Module.colors<-list(Signature=c(pats1="#009000",pats2="#f29a42",pats3="#0099b2",
                                pats4="#43639d",pats5="#a9182d",pats6="#b0b0b0",
                                pats7="#b486cb",pats8="#4fbb4b",pats9="#86ffd7"),
                    Dataset=c(dataset1="#a8d8ba",dataset2="#aebde1",dataset3="#afd8f4",
                              dataset4="#fba191",dataset5="#f5744d",dataset6="#b8b8b8",
                              dataset7="#cfb0d4",dataset8="#cbd8a2",dataset9="#ebeeb2"))


annrow<-data.frame(Modules.ann$path10)
rownames(annrow)<-rownames(Modules.ann)
colnames(annrow)<-"Signature"
anncol<-data.frame(ann.pats$dataset)
rownames(anncol)<-rownames(ann.pats)
colnames(anncol)<-"Dataset"

#-------
## Save heatmap
Mscore.matrix<-Mscore.matrix[rownames(Modules.ann),]

dev.off()
tiff(filename=paste0(opt$resultsPath,sep="","/Clustering_Figure1A__1.tiff"),res = 300,width = 8,height = 3.5,units="in")
pheatmap(Mscore.matrix[,rownames(ann.pats)],
         show_rownames = F,show_colnames = F, 
         cluster_rows = F,cluster_cols = F,
         breaks=seq(-2.5,2.5,length.out = 100),
         annotation_row = annrow,
         annotation_col = anncol,
         annotation_color=Module.colors,
         gaps_col = colGaps,
         gaps_row = rowGaps,
         color = colorRampPalette(c("#3a5f81","white","#d36c3a"))(100))

dev.off()


#-------
## Imporant Signatures by module



frec<-NULL
for(i in 1:nrow(Mscore.matrix)){
  frec<-c(frec,sum(ifelse(Mscore.matrix[i,]>=1.65,1,0)))
}

FrecPlot<-cbind(length(frec):1,frec/ncol(Mscore.matrix))
FrecPlot<-as.data.frame(FrecPlot)
FrecPlot<-cbind(FrecPlot,Modules.ann$color)
colnames(FrecPlot)<-c("x","frecuency","color")


save.image(paste0(opt$rdataPath,sep="","/Heatmap1A.RData"))

dev.off()
tiff(filename=paste0(opt$resultsPath,sep="","/Frecuency_Figure1A__1.tiff"),res = 300,width = 3,height = 3,units="in")
p1<-ggplot(FrecPlot,aes(x=x,y=frecuency,fill=color))+theme_classic()+coord_flip()+
  geom_point(color=FrecPlot$color,size=1.2)
plot(p1)
dev.off()

#-------
## Cleaning work space

rm(anncol,annrow,d,heatmap,Nb,results,resultsPats,stability,x,AA,category,colGaps,
   color,dataset,i,l,path10,rowGaps,title,frec,FrecPlot,p1)


Modules.list<-Modules.list[rownames(Modules.ann)]

#-------
## Save RData object
save.image(paste0(opt$rdataPath,sep="","/Clustering.RData"))






