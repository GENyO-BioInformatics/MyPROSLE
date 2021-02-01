# Overview
Here all the necessary functions are stored (in the R programming language) to perform the measurements of the Molecular dYsregulated  PROfiles for Systemic Lupus Erythematosus patients (MyPROSLE), in which it is estimated how affected a series of molecular pathways are in a patient individually with respect to a distribution of healthy controls.

![ ](bin/F1.png)

# Requirements
You will get file structure of the project and main scripts by using:

```
git clone https://github.com/GENyO-BioInformatics/MyPROSLE.git
```

or open it with GitHub Desktop (download it [here](https://desktop.github.com/)).

Also, you need to have ```R``` installed. You will find how to do it in the R project web (https://www.r-project.org/).


# How to use it

MyPROSLE estimates how affected are a serie of molecular pathways in a patient with respect to a distributions of healthy controls. The functions needed are located at [bin/Resources.R](bin/Resources.R).


# Inputs and Outputs of MyPROSLE

The input for MyPROSLE consist of an expression matrix of more than one patients or a numeric vector with gene expression of one single patient, an expression matrix of healthy samples and the module reference provided [here](data/ModuleReference.RData). Gene expression matrices have to have gene symbols in rows.

With these three input, you can run the function launch.MyPROSLE and obtain the MyPROSLE of your samples.

The output of MyPROSLE is a table with the molecular profile score of each module and the corresponding figure for each SLE sample. All output results are saved by default in results folder.

# Run example

To show you how use the MyPROSLE software, there is a script called [bin/run_examples.R](bin/run_examples.R). Example data is obtained from NCBI GEO serie (GSE65391).

```R

### Load functions and data
source("Resources.R")
Modules.path = "../data/ModuleReference.RData"
load("../data/DataExample.RData")
dir.create("../results")

### Single Patient

Patient = SLE[,1] # Select one sample from the example dataset. For this example we are going to select the first one
names(Patient) = rownames(SLE)


RES = launch.MyPROSLE(Modules.path,
                      Patient,
                      Healthy,
                      show.plot = TRUE, #FALSE: do not show the plot or TRUE: show the plot
                      save.plot = "Example_patient", #NULL: do not save the plot or character with the name of the plot. For example ("Example_patient")
                      save.as = "pdf") #Output format of the plot (pdf, png, jpg, tiff)



### Multiple Patient

Patient = SLE[,c(1,2,3,4,5)] # Select some samples from the example dataset. For this example we are going to select the first five samples

RES = launch.MyPROSLE(Modules.path,
                      Patient,
                      Healthy,
                      show.plot = TRUE, #FALSE: do not show the plot or #TRUE: show the plot
                      save.plot = TRUE, #NULL: do not save the plots or TRUE: save the plots with the sample names
                      save.as = "pdf") #Output format of the plot (pdf, png, jpg, tiff)
```
