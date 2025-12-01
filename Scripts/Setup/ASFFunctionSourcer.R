#R functions

source(paste(home, "/Scripts/R_Functions/areaOfinfection.R", sep = ''))
source(paste(home, "/Scripts/R_Functions/CullingOneRun.R", sep = ''))
source(paste(home, "/Scripts/R_Functions/FirstDetect.R", sep = ''))
source(paste(home, "/Scripts/R_Functions/FastMovement.R", sep = ''))
source(paste(home, "/Scripts/R_Functions/GetOutputs.R", sep = ''))
source(paste(home, "/Scripts/R_Functions/GetOutputs.R", sep = ''))
source(paste(home, "/Scripts/R_Functions/InitializeSounders.R", sep = ''))
source(paste(home, "/Scripts/R_Functions/SimulateOneRun.R", sep = ''))
source(paste(home, "/Scripts/R_Functions/StateChanges.R", sep = ''))
source(paste(home, "/Scripts/R_Functions/Make_Grid.R", sep = ''))
source(paste(home, "/Scripts/R_Functions/RunSimulationModel.R", sep = ''))
source(paste(home, "/Scripts/R_Functions/FOI_R.R", sep = ''))
source(paste(home, "/Scripts/R_Functions/sounderlocsZone.R", sep = ''))
source(paste(home, "/Scripts/R_Functions/sounderlocsSummarize.R", sep = ''))
source(paste(home, "/Scripts/R_Functions/ExtApparentPrev.R", sep = ''))
source(paste(home, "/Scripts/R_Functions/Surveillance.R", sep = ''))


#cpp Functions
Rcpp::sourceCpp(paste0(home,"/Scripts/cpp_Functions/Movement_Fast_Generalized.cpp"), verbose=TRUE)
Rcpp::sourceCpp(paste0(home,"/Scripts/cpp_Functions/Movement_Fast_RSFavail.cpp"), verbose=TRUE)
Rcpp::sourceCpp(paste0(home,"/Scripts/cpp_Functions/Fast_FOI_Matrix.cpp"), verbose=TRUE)
Rcpp::sourceCpp(paste0(home,"/Scripts/cpp_Functions/SpatialZones_fast.cpp"), verbose=TRUE)
Rcpp::sourceCpp(paste0(home,"/Scripts/cpp_Functions/FindCellfromCentroid.cpp"), verbose=TRUE)

