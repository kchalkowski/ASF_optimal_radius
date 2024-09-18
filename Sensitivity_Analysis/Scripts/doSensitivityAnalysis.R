
#############################
########## Purpose ##########
#############################

#The purpose of this script is to run replicates of ASF simulation model, and match output to ASF outbreak surveillance data.

#Output: sl.summaries

#########################################
########## Set output location ##########
#########################################

#set home directory, location of ASF_Optimal_Rdius
#home="local/path/to/ASF_Optimal_Radius/repo"
home="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/ASF_optimal_radius"
setwd(home)


########################################
########## Parameter settings ##########
########################################

#Parameter settings to vary:
#fast/slow movement- state switch 1,2
#density- 1.5,3,5
#B1- 0.001, 0.01, 0.1, 1
states=c(1,2)
densities=c(1.5,3,5)
B1s=c(0.001,0.01,0.1,1)

##################################
########## Script setup ##########
##################################

#load libraries
library(sf)
library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(ggplot2)

#######################################
########## Toggle parameters ##########
#######################################

#need sounderlocs in out.opts
out.opts=c("sounderlocs")
N=100 #number of samples to take for each week
nrep=10 #number of reps to run
buffer=0.4

#sl.summary.opts-- different options to summarize sounderlocs output
sl.summary.opts=c("SEIRCZ_total_apparent",
                  "SEIRCZ_zone_apparent")

#originally had some cell-leve summaries, but not informative
#cell-level output options:
#SEIRCZ_zone_cells, SEIRCZ_zone_cells_apparent

###########################
########## Setup ##########
###########################

#outer loop would start here
#set manually for test run
state=1
density=1.5
B1=0.001
  
#Run InitializeASFModel.R once per session

#Run SetParameters once per session, or as needed to update parameters
source(paste0(home,"/Sensitivity_Analysis/Scripts/SetParmsSensAnalysis.R"))

#loads all functions, variables needed, etc.
source(paste0(home,"/Sensitivity_Analysis/Scripts/InitializeSensAnalysis.R"))

#set name of run folder and create folder for run output
run_folder=gsub("\\..*", "", gsub("-|:| ","_",paste0("TestRun_",Sys.time())),perl=TRUE)
if(!dir.exists(paste0("Sensitivity_Analysis/Output/",run_folder))){dir.create(paste0("Sensitivity_Analysis/Output/",run_folder))}

#####################################
########## Loop replicates ##########
#####################################

for(rep in 1:nrep){
  
  print(rep)
  
  #run function, generate list of output based on out.opts
  out.list=RunSimulationModel(rep)
  #invalid time argument, rep S
  sl.summaries_rep=sounderlocsSummarize(out.list$sounderlocs,rep,sl.summary.opts,DetP,N)
    
  if(rep==1){
      sl.summaries=sl.summaries_rep
    } else{
      sl.summaries=lapply(seq_along(sl.summaries), function(x) rbind(sl.summaries[[x]], sl.summaries_rep[[x]]))
    }
  
}


#save output
saveRDS(sl.summaries,
        paste0("Sensitivity_Analysis/Output/",
               run_folder,
               "/sl.summaries.rds"))

