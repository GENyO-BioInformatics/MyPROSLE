##############################
## MyPROSLE 
## R version R version 4.0.4
## 
##############################
## Random stability

## @@@@@@@@@@@@@@@@@@@@@@@@@@@ Input options
## Input general options

load("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/RData/For_Random_stability.RData")


#-------
## Get distances between modules in each dataset (Randomized)
Dist.list<-list()
for(i in 1:length(DATA.Mscore)){
  
  data<-DATA.Mscore[[i]]$SLE
  data<-data[sample(1:nrow(data),replace = F),sample(1:ncol(data),replace = F)]
  
  Dist = (dist2(as.matrix(data),as.matrix(data)))^(1/2)
  Dist.list[[i]]<-Dist
}
rm(data,Dist,i)

K<-seq(10,30,by=1) ## number of neighbors, usually (10~30)
alpha<-seq(0.3,0.8,by=0.1) # alpha, hyperparameter, usually (0.3~0.8)

## 500 permutations selecting different number of datasets
permutations<-sample(2:(length(Dist.list)-1),500,replace=T)
permutations<-c(permutations,rep(length(Dist.list),2)) # To select also all datasets

print(table(permutations))
permutations<-as.list(permutations)

## For each permutation, all parameter combinations must be tested
combinations<-as.matrix(crossing(K,alpha)) ## 126 parameter combinations

cl<-makeCluster(detectCores()-1)

cat("\nMeasuring cluster stability, this process is slow...")
stability<-parLapply(cl,permutations,
                     Dist.list=Dist.list,
                     combinations=combinations,
                     ClusterStats)
melt.stability = as.data.frame(do.call("rbind", stability))

stopCluster(cl)
rm(combinations,permutations,K,alpha,Dist.list,cl)

## Get best number of clusters
nbplot<-aggregate(melt.stability$nbindex,by=list(datasets=melt.stability$datasets,clusters=melt.stability$clusters),data=melt.stability,FUN=mean)
colnames(nbplot)[3]<-"frecuency"


x<-aggregate(frecuency ~ clusters, data = nbplot, FUN = mean)

ggplot(data=x, aes(x=clusters, y=frecuency)) +
  geom_line(size=1,linetype="dashed",color="#333333")+xlim(c(2,14))+
   theme_classic()

setwd("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/RData")

saveRDS(x,"randomCurve.rds")

#dev.off()
#tiff(filename=paste0(opt$resultsPath,sep="","/Path_random_clustering.tiff"),res = 300,width = 6.5,height = 3.2,units="in")
#p1<-ggplot(data=nbplot, aes(x=clusters, y=frecuency, group=datasets,color=as.factor(datasets))) +
#  geom_line(size=1.1)+xlim(c(2,14))+
#  geom_point(size=2) + theme_classic()
#plot(p1)
#dev.off()


x<-readRDS("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/RData/randomCurve.rds")
nbplot<-readRDS("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/RData/nclustPlot.rds")


x$datasets<-rep(0,14)
x<-x[,c(3,1,2)]

M<-rbind(nbplot,x)

ltype<-ifelse(M$datasets==0,2,1)

M$ltype<-ltype

#setwd("D:/DATA/WORK/Toro.et.al.MyPROSLE")

dev.off()
tiff(filename=paste0(opt$resultsPath,sep="","/PathRandom.tiff"),res = 300,width = 5,height = 2.5,units="in")
p1<-ggplot(data=M, aes(x=clusters, y=frecuency, group=datasets,color=as.factor(datasets))) +
  geom_line(aes(group=datasets),linetype=M$ltype[nrow(M):1],size=0.9)+xlim(c(2,14))+ylim(0,60)+
  geom_point(size=1) + theme_classic()+
  scale_color_manual(values = c("0"="#666666","2"="#339900","3"="#336699","4"="#CC6666",
                                "5"="#CC99FF","6"="#99CCFF","7"="#99CCCC","8"="#66CC99",
                                "9"="#FF9900"))
plot(p1)
dev.off()



