
#############################
########## Purpose ##########
#############################

#The purpose of this script is to run the ASF simulation model


#state=1
#density=1.5
#ss=2

source(paste0(home,"/Scripts/SetParameters.R"))

for(rep in 1:100){

print(rep)
out.list=RunSimulationModel(rep)
#BB.rep=out.list$BB
if("sounderlocs"%in%out.opts){
  SEIRCZ.only=out.list$sounderlocs[,3:9]
  SEIRCZ.rep=SEIRCZ.only %>% 
    dplyr::group_by(timestep) %>% 
    dplyr::summarise(S=sum(S),E=sum(E),I=sum(I),R=sum(R),C=sum(C),Z=sum(Z)) %>% as.data.frame()
  SEIRCZ.rep$rep=rep
  SEIRCZ.rep=SEIRCZ.rep[,c(8,1:7)] #want rep in front
  if(rep==1){
    SEIRCZ.summary_214pm=SEIRCZ.rep
    #BB=BB.rep
  } else{
    SEIRCZ.summary_214pm=rbind(SEIRCZ.summary_214pm,SEIRCZ.rep)
    #BB=cbind(BB,BB.rep)
  }
}
}

path=paste0("/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/ASF_optimal_radius/Output/Pse_Compare_Trblsht/final/")
saveRDS(SEIRCZ.summary_214pm,paste0(path,"SEIRCZsummary_04JUN24_214pm.rds"))


#use this for testing
#outdir="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/ASF_optimal_radius/Output/Test_Run/"

outdir="/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/ASF_optimal_radius/Output/APR1124_Run/"

#loop parameters
reps=500

densities=c(1.5,3,5)
s.sizes=c(2,4,6)
states=c(1)

#innermost loop, num replicates
#matrows=reps*length(densities)
#matrows=reps
for(st in 1:length(states)){
  
  state=states[st]
  
for(d in 1:length(densities)){
  #d=1
  density=densities[d]
  ss=s.sizes[d]
  
  source(paste0(home,"/Scripts/SetParameters.R")) #uncomment this to set parms if not running in RunASFSimReplicates.R
  
  #Initialize general output mat, always gets returned
  gen.outmat=matrix(nrow=reps,ncol=12)
  colnames(gen.outmat)=c("state","density","rep","Tinc","sumTculled","idT","Mspread","IConDD","ICatEnd","TincToDD","TincFromDD","DET")
  
  for(rep in 1:reps){
  #r=rep+(reps*d-reps)
  cat("rep ", rep, "\n")
  out.list=RunSimulationModel(rep)
  #out.list2=RunSimulationModel()
  
  #general outmat summaries
  gen.outmat[rep,1]=states[st]
  gen.outmat[rep,2]=densities[d]
  gen.outmat[rep,3]=rep
  gen.outmat[rep,4]=out.list$Tinc
  gen.outmat[rep,5]=out.list$sumTculled
  gen.outmat[rep,6]=out.list$idT
  gen.outmat[rep,7]=out.list$Mspread
  gen.outmat[rep,8]=out.list$IConDD
  gen.outmat[rep,9]=out.list$ICatEnd
  gen.outmat[rep,10]=out.list$TincToDD
  gen.outmat[rep,11]=out.list$TincFromDD
  gen.outmat[rep,12]=out.list$DET
  
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
  
  ###test
  data=out.list$sounderlocs
  data.sum=data %>% group_by(x,y,timestep) %>% dplyr::summarise(locs=n()) %>% as.data.frame()
  data.sum[(data.sum$locs>1),]
  
  
  if("alldetections"%in%out.opts){
    alldet.rep=out.list$alldetections
    alldet.rep<-as.data.frame(alldet.rep)
    colnames(alldet.rep)=c("timestep","type","number","loc")
    alldet.rep=unique(alldet.rep[,c(1,2,3)])
    alldet.rep$type[alldet.rep$type==1]<-"live"
    alldet.rep$type[alldet.rep$type==0]<-"dead"
    alldet.rep$rep=rep
    alldet.rep<-alldet.rep[,c(4,1,2,3)]
    if(rep==1){
      alldet=alldet.rep
    } else{
      alldet=rbind(alldet,alldet.rep)
    }
  }
  
  if("incidence"%in%out.opts){
    incidence=out.list$incidence
    incidence<-incidence[,c(1,2)] #replicate summary, ignore locs
    inc.rep=incidence %>% dplyr::group_by(timestep,state) %>% dplyr::summarise(num=n()) %>% as.data.frame()
    inc.rep$rep=rep
    inc.rep=inc.rep[,c(4,1,2,3)]
    if(rep==1){
      inc.total=inc.rep
    } else{
      inc.total=rbind(inc.total,inc.rep)
    }
  }
  
}
  
  
  
  write.csv(gen.outmat,paste0(outdir,"gen.outmat_st",states[st],"_d",densities[d],".csv"))
  #if(d==1&st==1){
  #gen.outmat.compile=gen.outmat
  #} else{
  #  gen.outmat.compile=rbind(gen.outmat.compile,gen.outmat)
  #}
    
  if("sounderlocs"%in%out.opts){
  write.csv(SEIRCZ.summary,paste0(outdir,"SEIRCZ.summary_st",states[st],"_d",densities[d],".csv"))
  #if(d==1&st==1){
  #  SEIRCZ.summary.compile=SEIRCZ.summary
  #} else{
  #  SEIRCZ.summary.compile=rbind(SEIRCZ.summary.compile,SEIRCZ.summary)
  #}
  }
  
  if("alldetections"%in%out.opts){  
  write.csv(alldet,paste0(outdir,"alldet_st",states[st],"_d",densities[d],".csv"))
  #  if(d==1&st==1){
  #    alldet.compile=alldet
  #  } else{
  #    alldet.compile=rbind(alldet.compile,alldet)
  #  }
  }
  
  if("incidence"%in%out.opts){  
  write.csv(inc.total,paste0(outdir,"inc_st",states[st],"_d",densities[d],".csv"))
   # if(d==1&st==1){
  #    inc.total.compile=inc.total
  #  } else{
  #    inc.total.compile=rbind(inc.total.compile,inc.total)
  #  }
  }
  
}
  
}

#To output at end...
#write.csv(gen.outmat.compile,paste0(outdir,"Genoutmat.csv"))
#write.csv(SEIRCZ.summary.compile,paste0(outdir,"SEIRCZct.csv"))
#write.csv(alldet.compile,paste0(outdir,"Alldet.csv"))
#write.csv(inc.total.compile,paste0(outdir,"Inc.csv"))



#check output


#remove task callback when done! 
#removeTaskCallback(1)
#1a-3a
#start running 300 reps at 235AM, finished at 249AM
#~5min/100reps for slow/low density

