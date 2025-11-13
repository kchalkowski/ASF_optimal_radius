
## for troubleshooting
options(warn = 1,
        show.error.locations = TRUE,
        show.error.messages = TRUE)

#set directories
setwd(this.path::this.dir())
outdir<-file.path("Output")

#source targets file
source("_targets.R")

#Get pipeline
tar_manifest()

#Make pipeline
tar_make(callr_function = NULL, use_crew=FALSE, as_job=FALSE) # for troubleshooting
# tar_make()
# tar_make_clustermq(workers = 6)
