
#############################
########## Purpose ##########
#############################

#The purpose of this script is to run replicates of ASF simulation model

#########################################
########## Set output location ##########
#########################################

#set home directory, location of ASF_Optimal_Rdius
#home="local/path/to/ASF_Optimal_Radius/repo"
home="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/ASF_optimal_radius"
setwd(home)

###########################
########## Setup ##########
###########################

#Run SetParameters once per session, or as needed to update parameters
source(paste0(home,"/Scripts/Setup/SetParameters.R"))

#loads all functions, variables needed, etc.
source(paste0(home,"/Scripts/Setup/InitializeASFModel.R"))

#set number of replicates to run
nrep=100

#set name of run folder and create folder for run output
run_folder=gsub("\\..*", "", gsub("-|:| ","_",paste0("TestRun_",Sys.time())),perl=TRUE)
if(!dir.exists(paste0("Output/",run_folder))){dir.create(paste0("Output/",run_folder))}

#####################################
########## Loop replicates ##########
#####################################
rep=1
#run loop
for(rep in 1:nrep){

print(rep)
  
#run function, generate list of output based on out.opts
out.list=RunSimulationModel(rep)

#this is summarization method for out.opts with 'sounderlocs'
#produces summaries of numbers of SEIRCZ for each time step
#useful for sensitivity analyses
if("sounderlocs"%in%out.opts){
  SEIRCZ.only=out.list$sounderlocs[,3:9]
  SEIRCZ.rep=SEIRCZ.only %>% 
    dplyr::group_by(timestep) %>% 
    dplyr::summarise(S=sum(S),E=sum(E),I=sum(I),R=sum(R),C=sum(C),Z=sum(Z)) %>% as.data.frame()
  SEIRCZ.rep$rep=rep
  SEIRCZ.rep=SEIRCZ.rep[,c(8,1:7)] #want rep in front
  if(rep==1){
    SEIRCZ.summary=SEIRCZ.rep
  } else{
    SEIRCZ.summary=rbind(SEIRCZ.summary,SEIRCZ.rep)
  }
}
}


##################################
########## Save outputs ##########
##################################

saveRDS(out.list,paste0("/Output/",run_folder,"/out_list.rds"))
saveRDS(SEIRCZ.summary,paste0("/Output/",run_folder,"/SEIRCZ.summary.rds"))



