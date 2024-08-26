#R functions
source(paste(home, "/Scripts/Functions/areaOfinfection.R", sep = ''))
source(paste(home, "/Scripts/Functions/CullingOneRun.R", sep = ''))
source(paste(home, "/Scripts/Functions/FirstDetect.R", sep = ''))
source(paste(home, "/Scripts/Functions/GetOutputs.R", sep = ''))
source(paste(home, "/Scripts/Functions/GetOutputs.R", sep = ''))
source(paste(home, "/Scripts/Functions/InitializeSounders.R", sep = ''))
source(paste(home, "/Scripts/Functions/SimulateOneRun.R", sep = ''))
source(paste(home, "/Scripts/Functions/StateChanges.R", sep = ''))
source(paste(home, "/Scripts/Functions/Make_Grid.R", sep = ''))
source(paste(home, "/Scripts/Functions/RunSimulationModel.R", sep = ''))
source(paste(home, "/Scripts/Functions/FOI_R.R", sep = ''))

#cpp Functions
Rcpp::sourceCpp(paste0(home,"/Scripts/Functions/Movement_Parallel_Functionsmall.cpp"), verbose=TRUE)

