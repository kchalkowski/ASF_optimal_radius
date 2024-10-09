#R functions

source(paste(home, "/Scripts/Functions/areaOfinfection.R", sep = ''))
source(paste(home, "/Scripts/Functions/CullingOneRun.R", sep = ''))
source(paste(home, "/Scripts/Functions/FirstDetect.R", sep = ''))
source(paste(home, "/Scripts/Functions/FastMovement.R", sep = ''))
source(paste(home, "/Scripts/Functions/GetOutputs.R", sep = ''))
source(paste(home, "/Scripts/Functions/GetOutputs.R", sep = ''))
source(paste(home, "/Scripts/Functions/InitializeSounders.R", sep = ''))
source(paste(home, "/Scripts/Functions/SimulateOneRun.R", sep = ''))
source(paste(home, "/Scripts/Functions/StateChanges.R", sep = ''))
source(paste(home, "/Scripts/Functions/Make_Grid.R", sep = ''))
source(paste(home, "/Scripts/Functions/RunSimulationModel.R", sep = ''))
source(paste(home, "/Scripts/Functions/FOI_R.R", sep = ''))
source(paste(home, "/Scripts/Functions/sounderlocsZone.R", sep = ''))
source(paste(home, "/Scripts/Functions/sounderlocsSummarize.R", sep = ''))
source(paste(home, "/Scripts/Functions/ExtApparentPrev.R", sep = ''))

#cpp Functions
Rcpp::sourceCpp(paste0(home,"/Scripts/Functions/Movement_Fast_Generalized.cpp"), verbose=TRUE)
Rcpp::sourceCpp(paste0(home,"/Scripts/Functions/Movement_Fast_RSFavail.cpp"), verbose=TRUE)
Rcpp::sourceCpp(paste0(home,"/Scripts/Functions/Fast_FOI_Matrix.cpp"), verbose=TRUE)
Rcpp::sourceCpp(paste0(home,"/Scripts/Functions/SpatialZones_fast.cpp"), verbose=TRUE)
Rcpp::sourceCpp(paste0(home,"/Scripts/Functions/FindCellfromCentroid.cpp"), verbose=TRUE)

