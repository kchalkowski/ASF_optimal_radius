#R functions
source(paste(getwd(), "/Scripts/Functions/areaOfinfection.R", sep = ''))
source(paste(getwd(), "/Scripts/Functions/CullingOneRun.R", sep = ''))
source(paste(getwd(), "/Scripts/Functions/FirstDetect.R", sep = ''))
source(paste(getwd(), "/Scripts/Functions/GetOutputs.R", sep = ''))
source(paste(getwd(), "/Scripts/Functions/GetOutputs.R", sep = ''))
source(paste(getwd(), "/Scripts/Functions/InitializeSounders.R", sep = ''))
source(paste(getwd(), "/Scratch/Functions/Movement_ML_Exact.R", sep = ''))
source(paste(getwd(), "/Scripts/Functions/SimulateOneRun.R", sep = ''))
source(paste(getwd(), "/Scripts/Functions/Functions/StateChanges.R", sep = ''))
source(paste(getwd(), "/Scripts/Functions/Make_Grid.R", sep = ''))
source(paste(getwd(), "/Scripts/Functions/RunSimulationModel.R", sep = ''))
source(paste(getwd(), "/Scripts/Functions/FOI_R.R", sep = ''))

#cpp Functions
Rcpp::sourceCpp("./Scripts/Movement_Parallel_Functionsmall.cpp", verbose=TRUE)
Rcpp::sourceCpp("./Scripts/Fast_FOI_Parallel.cpp", verbose=TRUE)

