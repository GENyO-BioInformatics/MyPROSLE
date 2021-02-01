### Load useful functions, Modules and example dataset
source("Resources.R")
Modules.path = "ModuleReference.RData"
load("DataExample.RData")

### Single Patient

Patient = SLE[,1] # Select one sample from the example dataset. For this example we are going to select the first one
names(Patient) = rownames(SLE)


RES = launch.MyPROSLE(Modules.path,
                      Patient,
                      Healthy,
                      PreservationFilter=30,
                      show.plot = TRUE, #FALSE: do not show the plot or TRUE: show the plot
                      save.plot = "Example_patient", #NULL: do not save the plot or character with the name of the plot. For example ("Example_patient")
                      save.as = "pdf") #Output format of the plot (pdf, png, jpg, tiff)



### Multiple Patient

Patient = SLE[,c(1,2,3,4,5)] # Select some samples from the example dataset. For this example we are going to select the first five samples

RES = launch.MyPROSLE(Modules.path,
                      Patient,
                      Healthy,
                      PreservationFilter=30,
                      show.plot = TRUE, #FALSE: do not show the plot or #TRUE: show the plot
                      save.plot = TRUE, #NULL: do not save the plots or TRUE: save the plots with the sample names
                      save.as = "pdf") #Output format of the plot (pdf, png, jpg, tiff)



