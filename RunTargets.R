
#set directories
setwd(this.path::this.dir())
outdir<-file.path("Output")

#source targets file
source("_targets.R")

#Get pipeline
tar_manifest()

#Make pipeline
tar_make()

tar_visnetwork()

#list("SurvDetProb"=0.9,"Detp"=1)
