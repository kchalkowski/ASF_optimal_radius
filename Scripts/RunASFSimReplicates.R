###Structure to run replicates of ASF simulation model

#try keep same seed set, catch issue in loop
#addTaskCallback(function(...) {set.seed(123);TRUE})

#loop parameters
reps=500

densities=c(1.5,3,5)
s.sizes=c(2,4,6)
#innermost loop, num replicates
#matrows=reps*length(densities)
#matrows=reps
#outmat=matrix(nrow=0,ncol=11)
#colnames(outmat)=c("density","rep","Tinc","sumTculled","idT","Mspread","IConDD","ICatEnd","TincToDD","TincFromDD","DET")

#outmat.rep=matrix(nrow=matrows,ncol=11)
#colnames(outmat.rep)=c("density","rep","Tinc","sumTculled","idT","Mspread","IConDD","ICatEnd","TincToDD","TincFromDD","DET")

#for(d in 1:length(densities)){
  d=1
  density=densities[d]
  ss=s.sizes[d]
  
  source(paste0(home,"/Scripts/SetParameters.R")) #uncomment this to set parms if not running in RunASFSimReplicates.R
  incmat=matrix(nrow=72,ncol=reps)
for(rep in 1:reps){
  r=rep+(reps*d-reps)
  cat("rep ", rep, "\n")
  #out.list=RunSimulationModel()
  incmat[,rep]=RunSimulationModel()
  #outmat.rep[r,1]=densities[d]
  #outmat.rep[r,2]=rep
  #outmat.rep[r,3]=out.list$Tinc
  #outmat.rep[r,4]=out.list$sumTculled
  #outmat.rep[r,5]=out.list$idT
  #outmat.rep[r,6]=out.list$Mspread
  #outmat.rep[r,7]=out.list$IConDD
  #outmat.rep[r,8]=out.list$ICatEnd
  #outmat.rep[r,9]=out.list$TincToDD
  #outmat.rep[r,10]=out.list$TincFromDD
  #outmat.rep[r,11]=out.list$DET
  
}
  
  #outmat=rbind(outmat,outmat.rep)

#}

#run slow (FL) first
#slow.outmat=outmat.rep
#Tinc

#remove task callback when done!
removeTaskCallback(1)
#1a-3a
#start running 300 reps at 235AM, finished at 249AM
#~5min/100reps

#Next steps:
#1-write notes on what output mat cols are
#2-get parameter settings from manuscript
#3-run manuscript parameter settings

#Tinc #sum of all exposures over simulation 
#sum(Tculled)
#idT #last timestep there is an infectious individual
#Mspread #max spread of infection
#IConDD #number of I, C, and E on detection day
#ICatEnd #number of I,C,E on last day 
#TincFromDD #sum of all exposures starting day after detection day
#TincToDD #sum of all exposures up until detection day
#DET #total number of detections

#Note: seems like more infected individuals at end of sim for full models.. need compare w matlab script output
#double check ALL parms, incl fixed
