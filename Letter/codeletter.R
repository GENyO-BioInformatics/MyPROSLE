################################################################################
##
## R version 4.2.3
################################################################################
## Pre-procesing public LN datasets (GSE72326 and GSE99967)

# Files needed can be downloaded from:
# https://drive.google.com/drive/folders/1VscKg2blQVg7BJMF6nHaKpwX2R1tcOex?usp=sharing

## Set environment

set.seed(12345678)

library("GEOquery")
library("caret")
library("stringr")
library("stringi")
library("parallel")

source("utils.R")

Datasets.LN <- list()


## GSE72326 ····································································
gset <- getGEO("GSE99967", GSEMatrix = TRUE) ## Download. data from NCBI GEO
gset <- gset[[1]]

data <- exprs(gset) ## Get gene-expression matrix
data <- norm.log(data) ## Log2 transformation

nonVar.genes <- summ.var(data = data) ## Remove genes with near to zero variance
data <- data[!nonVar.genes$nzv, ]

## Annotate Probes to GeneSymbol
## Download GPL21970 db from NCBI GEO
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL21970
genome <- getGEO("GPL21970")@dataTable@table[, c("GeneSymbol", "ID")]
colnames(genome) <- c("toGenes", "fromGenes")

genome <- genome[genome$fromGenes != "", ] ## Filtering non anotategd genes
genome <- genome[genome$fromGenes %in% rownames(data), ]
genome <- na.omit(genome)
data <- data[genome$fromGenes, ]
finalGenes <- unique(genome$toGenes)

temp <- mclapply(finalGenes, calculate_medians,
                 genome = genome,
                 expressionMatrix = data, mc.cores = 1
)
temp <- as.data.frame(do.call("rbind", temp))
rownames(temp) <- finalGenes
colnames(temp) <- colnames(data)
data <- temp

metadata <- phenoData(gset) ## Download clinical data from NCBI GEO
metadata <- pData(metadata)

## Select relevant variables for the analysis
metadata <- metadata[, c(
  "title", "disease state:ch1",
  "nephritis state:ch1", "biopsy class:ch1"
)]
colnames(metadata) <- c("title", "diagnosis", "LN", "biopsy")

## Separe Healthy controls and SLE samples
sel <- ifelse(metadata$diagnosis == "control", F, T)
Healthy <- data[, rownames(metadata)[!sel]]
metadata <- metadata[sel, ]

SLE <- data[, rownames(metadata)]

write.table(SLE, file = "GSE99967_cases.tsv", sep = "\t", quote = F)
write.table(Healthy, file = "GSE99967_healthy.tsv", sep = "\t", quote = F)

metadata$LN <- ifelse(metadata$LN == "active LN (ALN)", "YES", "NO")
metadata$pLN <- ifelse(metadata$LN == "NO", "NO",
                       ifelse(metadata$biopsy == "V", "NO", "YES")
)

dataset <- list(SLE, Healthy, metadata)
names(dataset) <- c("Disease", "Healthy", "metadata")
Datasets.LN[["GSE99967"]] <- dataset

## Summary: GSE99967
## 29 samples with LN
## 23 with pLN
## 13 samples without LN
## 17 Healthy samples

rm(data, dataset, genome, gset, Healthy, metadata, nonVar.genes, SLE, temp, finalGenes, sel)



## GSE72326 ····································································
gset <- getGEO("GSE72326", GSEMatrix = TRUE) ## Download data from NCBI GEO
gset <- gset[[1]]

data <- exprs(gset) ## Get gene-expression matrix
data <- norm.log(data) ## Log2 transformation

nonVar.genes <- summ.var(data = data) ## Remove genes with near to zero variance
data <- data[!nonVar.genes$nzv, ]

## Annotate Probes to GeneSymbol using biomaRt
data <- annotateGenes(
  data = data, toGenes = "external_gene_name",
  fromGenes = "illumina_humanht_12_v4"
)

metadata <- phenoData(gset) ## Download clinical data from NCBI GEO
metadata <- pData(metadata)

## Separe Healthy controls and SLE samples
Healthy <- data[, rownames(metadata)[metadata$`group:ch1` == "Healthy control of SLE"]]

metadata <- metadata[ifelse(metadata$"group:ch1" == "SLE", T, F), ]

## Select relevant variables for the analysis
metadata <- metadata[, c(
  "title", "collection date:ch1", "group:ch1",
  "renal:ch1", "class kb:ch1", "ever renal:ch1", "condition:ch1", "sledai:ch1"
)]
colnames(metadata) <- c("title", "date", "diagnosis", "LN", "biopsy", "everRenal", "condition", "sledai")

patientInfo <- do.call("rbind", str_split(metadata$title, "V"))

metadata <- data.frame(
  "PatientID" = patientInfo[, 1],
  "Visit" = as.numeric(patientInfo[, 2]), metadata
)

metadata <- metadata[order(metadata$PatientID, metadata$Visit, decreasing = F), ]

## Date to numeric measurement
metadata$date <- stri_replace_all_regex(
  metadata$date,
  pattern = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"),
  replacement = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"),
  vectorize = FALSE
)

metadata$TimeFromVisits <- NA
patients <- unique(metadata$PatientID)
for (i in 1:length(patients)) {
  tmp <- metadata[metadata$PatientID == patients[i], ]
  tmp <- tmp[order(tmp$Visit, decreasing = F), ]
  ref <- tmp[1, "date"]
  
  for (j in 1:nrow(tmp)) {
    x <- as.numeric(difftime(strptime(tmp[j, "date"], format = "%d-%m-%Y"),
                             strptime(ref, format = "%d-%m-%Y"),
                             units = "days"
    ))
    metadata[metadata$PatientID == patients[i] &
               metadata$Visit == tmp$Visit[j], "TimeFromVisits"] <- x
  }
}

metadata$biopsyClass <- ifelse(is.na(metadata$biopsy), NA,
                               ifelse(str_detect(metadata$biopsy, "III") | str_detect(metadata$biopsy, "IV"), "Prolif", "Other")
)

SLE <- data[, rownames(metadata)]

write.table(SLE, file = "GSE72326_cases.tsv", sep = "\t", quote = F)
write.table(Healthy, file = "GSE72326_healthy.tsv", sep = "\t", quote = F)

dataset <- list(SLE, Healthy, metadata)
names(dataset) <- c("Disease", "Healthy", "metadata")
Datasets.LN[["GSE72326"]] <- dataset

saveRDS(Datasets.LN, "DatasetsLN.rds")


## GSE72326 is a longitudinal datasets
## Identify patients without Nephritis (and ever LN)
## 62 Longitudinal patients
## 25 Patients with Nephritis (15 with proliferative Nephritis)
## 37 Patients without LN in any visit (17 with ever Renal disease)


################################################################################
## Cohen's kappa calculation

## Use the gene expression matrix in MyPROSLE to obtain the prediction results ##

library("psych")
library("dplyr")

## Here we provide a function to calculate Cohen's kappa from the two datasets of interest

library(tidyverse)
prediction.with.healthy.GSE99967 <- read.delim("GSE99967_pred_with_healthy.tsv")
prediction.without.healthy.GSE99967 <- read.delim("GSE99967_pred_without_healthy.tsv")
prediction.with.healthy.GSE72326 <- read.delim("GSE72326_pred_with_healthy.tsv")
prediction.without.healthy.GSE72326 <- read.delim("GSE72326_pred_without_healthy.tsv")


cohenK_function <- function(prediction.with.healthy, prediction.without.healthy) {
  prediction.with.healthy <- prediction.with.healthy %>% filter(Variable == "Proliferative nephritis")
  prediction.with.healthy <- prediction.with.healthy[, seq(4, ncol(prediction.with.healthy))] %>% as.numeric()
  
  prediction.without.healthy <- prediction.without.healthy %>% filter(Variable == "Proliferative nephritis")
  prediction.without.healthy <- prediction.without.healthy[, seq(4, ncol(prediction.without.healthy))] %>% as.numeric()
  
  cohenK <- as.numeric(cohen.kappa(x = cbind(
    ifelse(prediction.with.healthy >= 0.5, "YES", "NO"),
    ifelse(prediction.without.healthy >= 0.5, "YES", "NO")
  ))$kappa)
  
  return(cohenK)
}

cohenK_function(prediction.with.healthy.GSE72326, prediction.without.healthy.GSE72326) # 0.652378
cohenK_function(prediction.with.healthy.GSE99967, prediction.without.healthy.GSE99967) # 0.285714

## Cohen's kappa was calculated from 16 Datasets downloaded from ADEx (
## GSE108497, GSE110174, GSE11907_SLE_GPL96, GSE11907_SLE_GPL97, GSE24706, GSE50772,GSE45291, GSE110169
## GSE51092,GSE61635, GSE65391, GSE72509, GSE80183.tsv, GSE82221_GPL10558, GSE84844, GSE90081),
## the Dataset from precisesads (witout SLE samples), an non public RA dataset (RNASeq)
## and 3 datasets downloaded and NCBI GEO(GSE72326, GSE99967, GSE17755)


################################################################################
## Select one sample for each patient (longitudinal dataset)

Datasets.LN <- readRDS("DatasetsLN.rds")
clin <- Datasets.LN$GSE72326$metadata

pats <- unique(clin$PatientID)
sel <- NULL
for (i in 1:length(pats)) {
  tmp <- clin[clin$PatientID == pats[i], ]
  #if ("Y" %in% tmp$LN) {
    tmp <- tmp[order(tmp$biopsy), ]
    sel <- c(sel, rownames(tmp)[1])
  # } else {
  #   sel <- c(sel, c(rownames(tmp)[sample(1:nrow(tmp), 1)]))
  # }
}
# First visit selection for each patient, from line 231 to 241
# For random selection of visit from patients without biopsy, remove comments from line 231 to 241

clin <- clin[sel, ]
clin$pLN <- ifelse(is.na(clin$biopsyClass), "NO",
                   ifelse(clin$biopsyClass == "Prolif", "YES", "NO")
)

Datasets.LN$GSE72326$metadata <- clin
Datasets.LN$GSE72326$Disease <- Datasets.LN$GSE72326$Disease[, rownames(clin)]

## Cohen first visit
selectColumns<-c("Variable","Category","Algorithm",sel)
cohenK_function(prediction.with.healthy.GSE72326[,selectColumns], 
                prediction.without.healthy.GSE72326[,selectColumns]) # 0.656

## aggretment
x<-prediction.with.healthy.GSE72326[prediction.with.healthy.GSE72326$Variable=="Proliferative nephritis",sel]
x<-ifelse(x>=0.5,"YES","NO")
y<-prediction.without.healthy.GSE72326[prediction.without.healthy.GSE72326$Variable=="Proliferative nephritis",sel]
y<-ifelse(y>=0.5,"YES","NO")

agreet <- (table(x == y)["TRUE"]) /
             (sum(table(x == y))) * 100
# 83.9

saveRDS(Datasets.LN, "DatasetsLN2.rds")





################################################################################
## Missing genes in data
library("dplyr")

# devtools::install_github("jordimartorell/pathMED") # Install wrapper of functions contained in MyPROSLE. The installation may take a long time

library("ggplot2")
library("pathMED")

## R Objects used in MyPROSLE web (used locally)
model.pLN <- readRDS("modelpLN.rds") ## Model to predict pLN
paths <- readRDS("SLEpaths.rds") ##  SLE-relevant paths


Datasets.LN <- readRDS("DatasetsLN.rds")
data <- Datasets.LN$GSE99967$Disease
healthy <- Datasets.LN$GSE99967$Healthy

sampl <- c(
  rep(0.9, 10), rep(0.8, 10), rep(0.7, 10), rep(0.6, 10), rep(0.5, 10),
  rep(0.4, 10), rep(0.3, 10), rep(0.2, 10), rep(0.1, 10)
)

m <- getMscores(
  Patient = data,
  Healthy = healthy,
  genesets = paths,
  cores = 1
)

originalPred <- predict(model.pLN, newdata = t(m), type = "prob")

res <- as.data.frame(matrix(ncol = 2, nrow = 0))

## This loop takes a long ##
for (i in 1:length(sampl)) {
  sampling <- sample(1:nrow(data), nrow(data) * sampl[i])
  dat <- data[sampling, ]
  
  tmp <- getMscores(
    Patient = dat,
    Healthy = healthy[rownames(dat), ],
    genesets = paths,
    cores = 10
  )
  
  pred <- predict(model.pLN, newdata = t(tmp), type = "prob")
  
  ## Agreement
  agreet <- (table(ifelse(originalPred$YES >= 0.5, "YES", "NO") == ifelse(pred$YES >= 0.5, "YES", "NO"))["TRUE"] /
               (sum(table(ifelse(originalPred$YES >= 0.5, "YES", "NO") == ifelse(pred$YES >= 0.5, "YES", "NO"))))) * 100
  
  res <- rbind(res, unname(c(agreet, sampl[i])))
  colnames(res) <- NULL
}
colnames(res) <- c("Agreement", "MissingGenes")
res$MissingGenes <- res$MissingGenes * 100

## Summarize results
res_agg <- res %>%
  mutate(MissingGenes = factor(MissingGenes)) %>%
  group_by(MissingGenes) %>%
  summarise(mean_agg = mean(Agreement))

res$MissingGenes<-paste0(res$MissingGenes,"%")
ggplot(res,aes(x=MissingGenes,y=Agreement)) + geom_jitter(alpha=0.7) + geom_boxplot()+theme_bw()


## Save results table
write.table(res_agg, "MissingGenes.txt", sep = "\t", quote = F)


################################################################################
