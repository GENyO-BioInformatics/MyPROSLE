##########################################################
## MyPROSLE. R functions
## R version 4.0.3
## daniel.toro@genyo.es
##########################################################

## Chech installed packages
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)>0){
    print(paste0("Insalling ", paste(new.pkg,collapse = ", ")))
    if (!requireNamespace("BiocManager", quietly = TRUE)){
      install.packages("BiocManager")
    }
    bioconductor_packages <- BiocManager::available()
    bioconductor_packages <- bioconductor_packages[bioconductor_packages %in% new.pkg]
    if (length(bioconductor_packages) > 0){
      BiocManager::install(bioconductor_packages, dependencies = TRUE,ask=FALSE)
    }
    new.pkg = new.pkg[!(new.pkg %in% bioconductor_packages)]
    if (length(new.pkg)>0){
      install.packages(new.pkg, dependencies = TRUE)
    }
  }
  res <- lapply(pkg,load.packages)
}

## Load packages
load.packages <- function(pkg){
  print(paste0("Loading ",pkg))
  suppressMessages(require(pkg,character.only = T))
}

## Molecular dysregulated profiles measurement for individual patients against healthy control distribution
Get.MyPROSLE<-function(Modules.path,    ## Path of ModuleReference.RData object
                       Patient,         ## Numeric vector of gene expression values of a patient (with Gene_Symbol as names)
                       Healthy,         ## Matrix of healthy samples (reference for gene expression distributions). Genes as rownames
                       PreservationFilter, ## Integer to filter modules by Preservation median rank
                       show.plot=T,     ## Print plot for each patient in screen
                       save.plot=NULL,  ## Save plot in current work directory (NULL: no save, character: name of plot(for one patient), TRUE: save plots with column names as plot names (for matrix of patients))
                       save.as="pdf"){  ## Extension ("pdf", "jpeg", "tiff", "png", "bmp", "svg")
  # Resources
  load(Modules.path)

  # Preparing sets
  genes = Reduce(intersect, list(rownames(Healthy),names(Patient),rownames(modules)))
  Patient<-Patient[genes]
  Healthy<-Healthy[genes,]
  modules<-modules[genes,]

  RES<-list()
  res<-rep(0,length(colorPalette)); names(res)<-names(colorPalette)

  for(path in 1:length(colorPalette)){
    genes.path<-as.character(modules[modules$modules==names(res)[path],1])
    genes.path<-intersect(genes.path,names(Patient))
    tmpRef<-Healthy[genes.path,]
    tmpPat<-Patient[genes.path]

    # Zscore calculation by module (only in modules with more than two correctly meassured genes)
    sel<-ifelse(apply(tmpRef,1,sd)==0,F,T)
    if("TRUE" %in% sel){
      if(table(sel)[["TRUE"]]>2){
        tmpRef<-tmpRef[sel,]
        tmpPat<-tmpPat[sel]
        Z.Pat<-(tmpPat-apply(tmpRef,1,mean))/apply(tmpRef,1,sd)
        res[path]<-mean(Z.Pat,na.rm=T)
      }
    }
  }

  # Preparing output
  nmes<-modules[!duplicated(modules$modules),]; rownames(nmes)<-nmes$modules; nmes<-nmes[names(res),]
  names(res)<-nmes$ModuleNames
  res<-res[37:1]
  res = data.frame(Modules = names(res),mprosle=res)

  filter<-ifelse(preservation$PreservationMedianRank<=PreservationFilter,T,F)
  filter<-as.character(preservation[filter,1])
  sel<-rownames(res) %in% filter
  res<-res[sel,]

  ann<-res; ann<-as.data.frame(ann)
  ann[,2]<-as.numeric(as.character(ann[,2]))
  colnames(ann)<-c("Modules","Zscore")
  ann[ann$Zscore>3,"Zscore"]<- 3
  ann[ann$Zscore< -3,"Zscore"]= -3
  ann$Modules <- factor(ann$Modules, levels = ann$Modules)
  RES[[1]]<-ann; names(RES)[1]<-"plot"
  RES[[2]]<-res; names(RES)[2]<-"mprosle"

  gplot<-ggplot(ann,aes(x=Modules,y=Zscore,fill=Modules)) + coord_flip()+
    theme_bw()+ ylim(-3.05,3.05)+
    theme(legend.position = "none",text = element_text(size=10))+
    geom_bar(stat = "identity", colour="black") +
    scale_fill_manual(values=c(as.character(colorPalette[sel])))

  if(show.plot){
    print(gplot)
  }

  if(!is.null(save.plot)){
    ggsave(filename= paste0(save.plot,".",save.as),plot = gplot,units="in",dpi=300,width = 6, height = 5.2)
  }

  return(RES)

} ## end

launch.MyPROSLE = function(Modules.path,    ## Path of ModuleReference.RData object
                           Patient,         ## Numeric vector of gene expression values of a patient (with Gene_Symbol as names)
                           Healthy,         ## Matrix of healthy samples (reference for gene expression distributions). Genes as rownames
                           PreservationFilter=30, ## Integer to filter modules by Preservation median rank
                           show.plot=F,
                           save.plot=NULL,  ## Save plot in current work directory (NULL: no save, character: name of plot for one sample, TRUE: save plots for multiple samples)
                           save.as="pdf"){
  check.packages(c("ggplot2"))
  load(Modules.path)
  if (is.vector(Patient)){
    RES = Get.MyPROSLE(Modules.path,    ## Path of ModuleReference.RData object
                       Patient,         ## Numeric vector of gene expression values of a patient (with Gene_Symbol as names)
                       Healthy,         ## Matrix of healthy samples (reference for gene expression distributions). Genes as rownames
                       PreservationFilter=PreservationFilter,
                       show.plot=show.plot,
                       save.plot=save.plot,  ## Save plot in current work directory (NULL: no save, character: name of plot)
                       save.as="pdf")
    RES = RES$mprosle
    rownames(RES) = NULL

    if(is.null(save.plot)){
      write.table(RES,"Results_patient.txt",sep="\t",row.names = F)
    }else{
      write.table(RES,paste0(save.plot,sep="",".txt"),sep="\t",row.names = F)
    }


  } else{
    res = data.frame(matrix(nrow = length(unique(modules$modules)),ncol = ncol(Patient)+1))
    colnames(res) = c("Modules",colnames(Patient))
    rownames(res) = unique(modules$ModuleNames)
    res$Modules = unique(modules$ModuleNames)
    for (patient in colnames(Patient)){
      sPatient = Patient[,patient]
      names(sPatient) = rownames(Patient)
      if (!is.null(save.plot)){
        RES = Get.MyPROSLE(Modules.path,    ## Path of ModuleReference.RData object
                         sPatient,         ## Numeric vector of gene expression values of a patient (with Gene_Symbol as names)
                         Healthy,         ## Matrix of healthy samples (reference for gene expression distributions). Genes as rownames
                         PreservationFilter=PreservationFilter,
                         show.plot=show.plot,
                         save.plot=patient,  ## Save plot in current work directory (NULL: no save, character: name of plot)
                         save.as=save.as)
      } else{
        RES = Get.MyPROSLE(Modules.path,    ## Path of ModuleReference.RData object
                           sPatient,         ## Numeric vector of gene expression values of a patient (with Gene_Symbol as names)
                           Healthy,         ## Matrix of healthy samples (reference for gene expression distributions). Genes as rownames
                           PreservationFilter=PreservationFilter,
                           show.plot=show.plot,
                           save.plot=NULL,  ## Save plot in current work directory (NULL: no save, character: name of plot)
                           save.as=save.as)
      }
      mprosle = RES$mprosle
      rownames(mprosle) = mprosle$Modules
      mprosle = mprosle[rownames(res),]
      res[,patient] = mprosle$mprosle
    }
    RES = res
    rownames(RES) = NULL

    RES<-na.omit(RES)

    write.table(RES,"Results.txt",sep="\t",row.names = F)
  }




  return(RES)
}
