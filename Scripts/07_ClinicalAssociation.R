##############################
## MyPROSLE 
## R version R version 4.0.4
## 
##############################
## PredictiveModels - Get clinical information


##-------------------------------------------------- (STEP 0)
## Load Mscores and Msignatures

set.seed(123456788)

load("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/RData/Clustering.RData")
source(paste0(opt$scriptPath,"/Resources.R",sep=""))

DATA.Msig<-GetSignatures(listM = DATA.Mscore,ann = Modules.ann)

tmp<-list()
for(i in 1:length(DATA.Mscore)){
  tmp[[i]]<-DATA.Mscore[[i]]$SLE
}
names(tmp)<-names(DATA.Mscore)
DATA.Mscore<-tmp;

rm(list=setdiff(ls(),c("DATA.Mscore","DATA.Msig","Modules.ann","opt")))


##-------
## List to store variable information
variable.list<-list()
pnt<-1

##-------------------------------------------------- (STEP 1)
## Clinical manifestations

#source(paste0(opt$scriptPath,"/Resource_models.R",sep=""))

## Clinical data
load(paste0(opt$dataPath,sep="","/Clin_pascualSLE.RData"))
load(paste0(opt$dataPath,sep="","/Clin_precisesadsSLE.RData"))
load(paste0(opt$dataPath,sep="","/Clin_petriSLE.RData"))


Clin.pascual<-Clin.pascual[intersect(colnames(DATA.Mscore$dataset1),rownames(Clin.pascual)),]
Clin.petri<-Clin.petri.visit[intersect(colnames(DATA.Mscore$dataset9),rownames(Clin.petri.visit)),]
rm(Clin.petri.patient,Clin.petri.visit)
Clin.precisesads<-Clin.precisesads[intersect(colnames(DATA.Mscore$dataset8),rownames(Clin.precisesads)),]

clinLIST<-list(Clin.pascual,Clin.petri,Clin.precisesads)
names(clinLIST)<-c("pascual","petri","precisesads")


##--------
## Fit data values

val2vals<-c("0"="NO","8"="YES","4"="YES","1"="YES","2"="YES",
            "No"="NO","Past"="X","Present"="YES","Unknown"="X",
            "Moderate"="YES","Severe"="YES","Yes"="YES")

## VASCULAR
vars<-c("Vasculitis","Vasc.1") ##1
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Vascular_Vasculitis"; pnt<-pnt+1}



## MUSCLE/SKELETAL
vars<-c("Arthritis","MS.1","MUSCLE_AND_SKELETAL_ARTHRITIS") ##2
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Musculoskeletal_Arthritis"; pnt<-pnt+1}

vars<-c("Myositis","MS.2","MUSCLE_AND_SKELETAL_EVIDENCE_OF_INFLAMMATORY_MYOPATHY") ##3
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Musculoskeletal_Myositis"; pnt<-pnt+1}

vars<-c("MUSCLE_AND_SKELETAL_MUSCLE_WEAKNESS") ##4
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Musculoskeletal_MuscleWeakness"; pnt<-pnt+1}



## RENAL
vars<-c("UrinaryCast","Renal.1") ##5
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Renal_UrinaryCast"; pnt<-pnt+1}

vars<-c("Hematuria","Renal.2") ##6
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Renal_Haematuria"; pnt<-pnt+1}

vars<-c("Proteinuria","Renal.3","KIDNEY_PROTEINURIA") ##7
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Renal_Proteinuria"; pnt<-pnt+1}

vars<-c("Pyuria","Renal.4") ##8
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Renal_Pyuria"; pnt<-pnt+1}

vars<-c("KIDNEY_ABNORMAL_CREATININE") ##27
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Renal_Abnormal_Creatinine"; pnt<-pnt+1}

vars<-c("KIDNEY_ABNORMAL_LIPID_PROFILE") ##28
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Renal_Abnormal_Lipid_Profile"; pnt<-pnt+1}

vars<-c("KIDNEY_ABNORMAL_URINE_ANALYSIS") ##29
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Renal_Abnormal_Urine_Analysis"; pnt<-pnt+1}

vars<-c("KIDNEY_BIOPSY_PROVEN_NEPHRITIS") ##30
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Renal_Biosy_Proven_Nephritis"; pnt<-pnt+1}



## DERMAL
vars<-c("NewRash","Skin.1") ##9
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Dermal_Rash"; pnt<-pnt+1}

vars<-c("Alopecia","Skin.2") ##10
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Dermal_Alopecia"; pnt<-pnt+1}

vars<-c("Mucosal_ulcers","Skin.3") ##11
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Dermal_MucosalUlcers"; pnt<-pnt+1}

vars<-c("SKIN_CALCINOSIS_CUTIS") ##34
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Dermal_Calcinosis_Cutis"; pnt<-pnt+1}

vars<-c("SKIN_APHTOUS_ULCERS") ##35
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Dermal_Aphtous_Ulcers"; pnt<-pnt+1}

vars<-c("SKIN_TELANGECTASIA") ##36
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Dermal_telangectasia"; pnt<-pnt+1}



## HEART
vars<-c("Pleurisy","Sero.1") ##12
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Heart_Pleurisy"; pnt<-pnt+1}

vars<-c("Pericarditis","Sero.2","HEART_PERICARDITIS") ##13
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Heart_Pericarditis"; pnt<-pnt+1}

vars<-c("HEART_VALVE_LESIONS") ##26
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Heart_Valve_Lesions"; pnt<-pnt+1}



## IMMUNOLOGICAL
vars<-c("LowComplement","Imm.1") ##14
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Immunological_LowComplement"; pnt<-pnt+1}

vars<-c("LAB_REDUCED_C3_LEVELS") ##31
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Immunological_lowC3"; pnt<-pnt+1}

vars<-c("LAB_REDUCED_C4_LEVELS") ##32
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Immunological_lowC4"; pnt<-pnt+1}

vars<-c("Increased_DNAbinding","Imm2") ##15
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Immunological_Increased_DNAbinding"; pnt<-pnt+1}



## CONSTITUTIONAL
vars<-c("Fever","Cons.1","SYMPTOM_FEVER") ##16
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Constitutional_Fewer"; pnt<-pnt+1}



## HAEMATOLOGICAL
vars<-c("Thrombocytopenia","Heme.1") ##17
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = c("0"="NO","1"="YES"))
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Haematological_Trombocytopenia"; pnt<-pnt+1}

vars<-c("Leukopenia","Heme.2") ##18
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Haematological_Leukopenia"; pnt<-pnt+1}

vars<-c("SYMPTOM_HYPERGAMMABULINEMIA") ##37
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Haematological_Hypergammabulinemia"; pnt<-pnt+1}



## COMORBIDITY
vars<-c("COMORBIDITY_ABDOMINAL_PAIN") ##19
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Comorbidity_Abdominal_Pain"; pnt<-pnt+1}

vars<-c("COMORBIDITY_ASTHMA") ##20
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Comorbidity_Asthma"; pnt<-pnt+1}

vars<-c("COMORBIDITY_DIARRHEA__RECURRENT_") ##21
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Comorbidity_Diarrhea_Recurrent"; pnt<-pnt+1}

vars<-c("COMORBIDITY_DYSLIPIDEMIA") ##22
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Comorbidity_Dyslipidemia"; pnt<-pnt+1}

vars<-c("COMORBIDITY_OBESITY__BMI____30_") ##23
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Comorbidity_Obesity"; pnt<-pnt+1}

vars<-c("COMORBIDITY_STIPSIS_CONSTIPATION") ##24
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Comorbidity_Stipsis_Constipation"; pnt<-pnt+1}

vars<-c("COMORBIDITY_THYROIDITIS") ##25
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Comorbidity_Thyroiditis"; pnt<-pnt+1}


## LUNG
vars<-c("LUNG_FUNCTIONAL_VENTILATORY_RESTRICTION") ##33
x<-MergeClinical(variableIDs=vars,clin=clinLIST,minSize=10,
                 Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                 val2val = val2vals)
if(is.null(x)==F){variable.list[[pnt]]<-x; 
names(variable.list)[pnt]<-"Lung_Functional_Ventilatory_Restriction"; pnt<-pnt+1}


#for(i in 1:length(variable.list)){
#  print(table(variable.list[[i]]$mscore$group))
#}

##-------------------------------------------------- (STEP 2)
## Autoantibodies

library("stringr")

## Load data
clin.aa<-read.csv(paste0(opt$dataPath,sep="","//auto_antibodies.txt"),header=T,sep="\t",dec=",")
rownames(clin.aa)<-clin.aa$OMICID

clin.aa<-clin.aa[intersect(colnames(DATA.Mscore$dataset8),rownames(clin.aa)),
                 ifelse(str_detect(string = colnames(clin.aa),pattern = "CALL",negate = F),T,F)]
colnames(clin.aa)<-gsub(pattern = "_CALL",replacement = "",x = colnames(clin.aa))
clinLIST<-list(clin.aa)

for(i in 1:ncol(clin.aa)){
  
  x<-MergeClinical(variableIDs=colnames(clin.aa)[i],clin=clinLIST,minSize=10,
                   Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="categoric",
                   val2val = c("negative"="NO","positive"="YES"))
  if(is.null(x)==F){variable.list[[pnt]]<-x; 
  names(variable.list)[pnt]<-colnames(clin.aa)[i]; pnt<-pnt+1}
  
}


##-------------------------------------------------- (STEP 3)
## Cytokines

##-------
## Load data
clin.cyt<-read.csv(paste0(opt$dataPath,sep="","/cytokines.txt"),header=T,sep="\t",dec=".")
rownames(clin.cyt)<-clin.cyt$OMICID

clin.cyt<-clin.cyt[intersect(colnames(DATA.Mscore$dataset8),rownames(clin.cyt)),c(106:ncol(clin.cyt))]


##-------
## Select the most measured cytokines
varSel<-c("BAFF_ELISA","BLC","CRP_ELISA","FAS_LIGAND","GDF_15",
            "IL_1_RA","IL_1_RII","IL_6_ELISA","IP_10","MCP_2","MCP_4",
            "MIP_1_BETA","MMP_2_ELISA","MMP_8","TARC","TGF_BETA_ELISA",
            "TNF_ALPHA_ELISA","TNF_RI")
sel<-(colnames(clin.cyt) %in% varSel)
clin.cyt<-clin.cyt[,sel]
rm(varSel,sel)

colnames(clin.cyt)<-gsub(pattern = "_ELISA",replacement = "",x = colnames(clin.cyt))
clinLIST<-list(clin.cyt)

for(i in 1:ncol(clin.cyt)){
  
  x<-MergeClinical(variableIDs=colnames(clin.cyt)[i],clin=clinLIST,minSize=10,
                   Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="numeric",
                   val2val = NULL)
  if(is.null(x)==F){variable.list[[pnt]]<-x; 
  names(variable.list)[pnt]<-colnames(clin.cyt)[i]; pnt<-pnt+1}
  
}

##-------------------------------------------------- (STEP 4)
## cell populations

clin.cell<-read.csv(paste0(opt$dataPath,sep="","/FlowCytometry.csv"),header=T,sep=";",dec=".")
rownames(clin.cell)<-clin.cell$OMICID
clin.cell<-clin.cell[,-1]

clin.cell<-clin.cell[intersect(colnames(DATA.Mscore$dataset8),rownames(clin.cell)),]

clin.cell<-clin.cell[ifelse(is.na(clin.cell$PBMC)==FALSE &is.na(clin.cell$PMN)==FALSE,T,F),]

total<-clin.cell$PBMC + clin.cell$PMN
for(i in 1:ncol(clin.cell)){
  clin.cell[,i]<-clin.cell[,i]/total
}
clin.cell<-clin.cell[,-c(18,19,20,21,22)]

clinLIST<-list(clin.cell)

for(i in 1:ncol(clin.cell)){
  
  x<-MergeClinical(variableIDs=colnames(clin.cell)[i],clin=clinLIST,minSize=10,
                   Mscore.list=DATA.Mscore,Msig.list=DATA.Msig,variable.type="numeric",
                   val2val = NULL)
  if(is.null(x)==F){variable.list[[pnt]]<-x; 
  names(variable.list)[pnt]<-colnames(clin.cell)[i]; pnt<-pnt+1}
  
}

# "D:/DATA/WORK/Toro.et.al.MyPROSLE/RData"

rm(list=setdiff(ls(),c("variable.list")))

save.image("D:/DATA/WORK/Toro.et.al.MyPROSLE.2022/RData/ClinicalVariables.RData")






















