#R functions
source(paste(getwd(), "/Scripts/areaOfinfection.R", sep = ''))
source(paste(getwd(), "/Scripts/CullingOneRun.R", sep = ''))
source(paste(getwd(), "/Scripts/FirstDetect.R", sep = ''))
source(paste(getwd(), "/Scripts/GetOutputs.R", sep = ''))
#source(paste(getwd(), "/Scripts/GenerateFakeStateData.R", sep = ''))
source(paste(getwd(), "/Scripts/GetOutputs.R", sep = ''))
source(paste(getwd(), "/Scripts/InitializeSounders.R", sep = ''))
#source(paste(getwd(), "/Scripts/FastMovement.R", sep = ''))
source(paste(getwd(), "/Scratch/Movement_ML_Exact.R", sep = ''))
source(paste(getwd(), "/Scripts/SimulateOneRun.R", sep = ''))
source(paste(getwd(), "/Scripts/StateChanges.R", sep = ''))
source(paste(getwd(), "/Scripts/Make_Grid.R", sep = ''))
source(paste(getwd(), "/Scripts/RunSimulationModel.R", sep = ''))
source(paste(getwd(), "/Scripts/FOI_R.R", sep = ''))

#source(paste(getwd(), "/Scripts/FOI.R", sep = ''))

#for testing RCPP compilation
#home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/ASF_optimal_radius/"
#setwd(home)

#cpp Functions
Rcpp::sourceCpp("./Scripts/Movement_Parallel_Functionsmall.cpp", verbose=TRUE)
Rcpp::sourceCpp("./Scripts/Fast_FOI_Parallel.cpp", verbose=TRUE)
#Rcpp::sourceCpp("./Scripts/Fast_FOI_Parallel_Archive21MAY24.cpp", verbose=TRUE)

#June 3rd, 2024
#cpp FOI version results not same as R
#For now, replace with exact version to ML version, optimize later
#source(paste(getwd(), "/Scripts/FOI_R.R", sep = ''))
#source(paste(getwd(), "/Scratch/Movement_ML_Exact.R", sep = ''))
#Rcpp::sourceCpp("./Scratch/Fast_FOI.cpp", verbose=TRUE)
#pkgbuild::check_build_tools(debug = TRUE)
