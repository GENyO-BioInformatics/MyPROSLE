## Log Normalization ···························································
#' @data Gene expression matrix
norm.log<-function(data){
  qx <- as.numeric(quantile(data, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) { data[which(data <= 0)] <- NaN
  data <- log2(data) }
  return(data)
}

## Calculate median of expression by gene ······································
#' @gene Gene to annotate
#' @genome Table with gene-probe set conections
#' @expressionMatrix Gene-expression data
calculate_medians = function(gene,
                             genome,
                             expressionMatrix){
  require("matrixStats")
  probe_ids = genome[genome$toGenes==gene,"fromGenes"]
  if (length(probe_ids)>1){
    res =  colMedians(as.matrix(expressionMatrix[probe_ids,]))
  } else{
    res <- as.numeric(expressionMatrix[probe_ids,])  
  }
  return(res)
}

## Function to annotate genes from a gene-expression matrix·····················
#' @data Gene expression matrix
#' @toGenes Gene ID
#' @fromGenes Probe ID
annotateGenes<-function(data,
                        toGenes='external_gene_name',
                        fromGenes='ensembl_gene_id'){
  require("parallel")
  require("biomaRt")
  
  mart = useMart("ensembl", dataset = paste0(casefold("hsapiens"),"_gene_ensembl"),host="https://jul2019.archive.ensembl.org")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host="https://jul2019.archive.ensembl.org")
  
  genome = getBM(attributes=c(toGenes,fromGenes),mart = human)
  colnames(genome)<-c("toGenes","fromGenes")
  genome <- genome[genome$fromGenes != "",]
  genome = genome[genome$fromGenes %in% rownames(data),]
  data = data[genome$fromGenes,]
  finalGenes = unique(genome$toGenes)
  
  
  temp = mclapply(finalGenes,calculate_medians,genome = genome,
                  expressionMatrix = data,mc.cores = 1)
  temp = as.data.frame(do.call("rbind", temp))
  rownames(temp) = finalGenes
  colnames(temp) = colnames(data)
  
  return(temp)
}


## Summary of variance of genes ················································
#' @data Gene expression matrix
summ.var<-function(data){
  require(caret)
  data<-t(data)
  nzv <- nearZeroVar(data, saveMetrics= TRUE)
  return(nzv)
}
