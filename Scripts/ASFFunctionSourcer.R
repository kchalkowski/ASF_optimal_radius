#R functions
source(paste(getwd(), "/Scripts/areaOfinfection.R", sep = ''))
source(paste(getwd(), "/Scripts/CullingOneRun.R", sep = ''))
source(paste(getwd(), "/Scripts/FirstDetect.R", sep = ''))
source(paste(getwd(), "/Scripts/GetOutputs.R", sep = ''))
source(paste(getwd(), "/Scripts/GenerateFakeStateData.R", sep = ''))
source(paste(getwd(), "/Scripts/GetOutputs.R", sep = ''))
source(paste(getwd(), "/Scripts/InitializeSounders.R", sep = ''))
source(paste(getwd(), "/Scripts/FastMovement.R", sep = ''))
source(paste(getwd(), "/Scripts/SimulateOneRun.R", sep = ''))
source(paste(getwd(), "/Scripts/StateChanges.R", sep = ''))
source(paste(getwd(), "/Scripts/FOI.R", sep = ''))

#cpp Functions
Rcpp::sourceCpp("./Scripts/Movement_Parallel_Functionsmall.cpp", verbose=TRUE)
Rcpp::sourceCpp("./Scripts/Fast_FOI_Parallel.cpp", verbose=TRUE)

#Rcpp::sourceCpp("./Scratch/Fast_FOI.cpp", verbose=TRUE)
