
#############################
########## Purpose ##########
#############################

#The purpose of this script is to run replicates of ASF simulation model, and match output to ASF outbreak surveillance data.

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

#Run InitializeASFModel.R once per session

#Run SetParameters once per session, or as needed to update parameters
source(paste0(home,"/Scripts/Setup/SetParameters.R"))

#loads all functions, variables needed, etc.
source(paste0(home,"/Scripts/Setup/InitializeASFModel.R"))

#load additional sens. analysis function
source(paste0(home,"/Sensitivity_Analysis/Scripts/ExtApparentPrev.R"))

#set number of replicates to run
nrep=3

#set name of run folder and create folder for run output
run_folder=gsub("\\..*", "", gsub("-|:| ","_",paste0("TestRun_",Sys.time())),perl=TRUE)
if(!dir.exists(paste0("Sensitivity_Analysis/Output/",run_folder))){dir.create(paste0("Sensitivity_Analysis/Output/",run_folder))}

##################################################
########## Toggle parameters externally ##########
##################################################

#need sounderlocs in out.opts
out.opts=c("sounderlocs")
DetP=0.9
N=100

#####################################
########## Loop replicates ##########
#####################################
#rep=1
#run loop
for(rep in 1:nrep){
  
  print(rep)
  
  #run function, generate list of output based on out.opts
  out.list=RunSimulationModel(rep)
  
    SEIRCZ.only=out.list$sounderlocs[,3:9]
    SEIRCZ.rep=SEIRCZ.only %>% 
      dplyr::group_by(timestep) %>% 
      dplyr::summarise(S=sum(S),E=sum(E),I=sum(I),R=sum(R),C=sum(C),Z=sum(Z)) %>% as.data.frame()
    SEIRCZ.rep$rep=rep
    SEIRCZ.rep=SEIRCZ.rep[,c(8,1:7)] #want rep in front
    SI_apparent.rep=ExtApparentPrev(SEIRCZ.rep,DetP,N,rep)
    if(rep==1){
      SEIRCZ.summary=SEIRCZ.rep
      SI_apparent.summary=SI_apparent.rep
    } else{
      SEIRCZ.summary=rbind(SEIRCZ.summary,SEIRCZ.rep)
      SI_apparent.summary=rbind(SI_apparent.summary,SI_apparent.rep)
    }
}


#save output
write.csv(SEIRCZ.summary,paste0("Sensitivity_Analysis/Output/",run_folder,"/SEIRCZ.summary.csv"),row.names=FALSE)
write.csv(SI_apparent.summary,paste0("Sensitivity_Analysis/Output/",run_folder,"/SI_apparent.summary.csv"),row.names=FALSE)

####################################################################
########## Match format to ASF outbreak surveillance data ##########
####################################################################
#SI_apparent.summary=read.csv(paste0("Sensitivity_Analysis/Output/",run_folder,"/SI_apparent.summary.csv"))

#summarize SI_apparent.summary across replicates
SI_apparent.summary[is.na(SI_apparent.summary)]<-0
SI_apparent.summary=SI_apparent.summary %>% dplyr::group_by(src,week) %>% dplyr::summarise(freq_prev=mean(freq_prev),freq_prev_q25=quantile(freq_prev,0.25),freq_prev_q75=quantile(freq_prev,0.75),appar_prev=mean(appar_prev),appar_prev_q25=quantile(appar_prev,0.25),appar_prev_q75=quantile(appar_prev,0.75))

#read ASF outbreak data
edat=readRDS(paste0("Sensitivity_Analysis/Input/datE80km_weekly.rds"))
wdat=readRDS(paste0("Sensitivity_Analysis/Input/datW80km_weekly.rds"))

edat2=edat[,c(1,2,3,9)]
wdat2=wdat[,c(1,2,3,9)]

#rename cols to match
colnames(edat2)[4]<-"appar_prev"
colnames(wdat2)[4]<-"appar_prev"
#colnames(SI_apparent.summary)[2]<-"week"

add_vars <- function(df,vars,val){
  df[vars] <- val
  df
}

edat2$type="observed_east"
wdat2$type="observed_west"
SI_apparent.summary$type="simulated"

new_vars=colnames(SI_apparent.summary)[!(colnames(SI_apparent.summary)%in%colnames(edat2))]

edat2=edat2 |>
  add_vars(df = _,vars = new_vars,val = NA)
wdat2=wdat2 |>
  add_vars(df = _,vars = new_vars,val = NA)

new_vars=colnames(edat2)[!(colnames(edat2)%in%colnames(SI_apparent.summary))]
SI_apparent.summary=SI_apparent.summary |>
  add_vars(df = _,vars = new_vars,val = NA)

l<-list(observed_east=edat2,
        observed_west=wdat2,
        simulated=SI_apparent.summary)

sadf=do.call(rbind, lapply(l, function(x) x[match(names(l[[1]]), names(x))]))

################################
########## Make plots ##########
################################

library(ggplot2)
live=sadf[sadf$src=="alive"|sadf$src=="live",]
dead=sadf[sadf$src=="carcass",]
ggplot(live, aes(x=week,y=appar_prev,color=type))+geom_line()
ggplot(dead, aes(x=week,y=appar_prev,color=type))+geom_line()









##################################################
########## Add method to get zone calcs ##########
##################################################


#out.list$sounderlocs
#centroids
#rl=out.list$sounderlocs
findZone_sim<-function(rl,type){
rl=st_as_sf(rl,coords=c("x","y"))
case1=rl[rl$timestep==1&rl$I>0,]
xw <- raster(ncol=200, nrow=200, xmn=0, xmx=80, ymn=0, ymx=80,crs="")
values(xw)=1
cns=extract(xw,case1,cellnumbers=TRUE)
xy_first=xyFromCell(xw,cns[,1])
case1

weeks=unique(rl$timestep)[order(unique(rl$timestep))]
rl$inZone=0
rl$seq=1:nrow(rl)
for(w in 1:length(weeks)){
  print(paste0("starting timestep ",w))
  rlw=rl[rl$timestep==w,]
  if(type=="I"){
  rl.pos=rlw[rlw$I>0,]
  }
  if(type=="C"){
  rl.pos=rlw[rlw$C>0]
  }
  if(w==1){
    pts=c(st_geometry(rl.pos),st_geometry(case1))
  } else{
    pts=c(pts,st_geometry(rl.pos))
  }
  
  wkpos_chull=st_convex_hull(st_union(pts)) %>% st_sf %>% st_cast()
  wkpos_chull=st_buffer(wkpos_chull,0.4)
  #cn_chull=terra::extract(x,wkpos_chull,cellnumbers=TRUE)
  #cn_inzone=cn_chull[[1]][,1]
  ints=st_intersection(rl[rl$timestep==w,],wkpos_chull)
  #rl[(rl$seq%in%ints$seq),]$inZone=1
  
  if(nrow(ints)>0){
  rl[(rl$seq%in%ints$seq),]$inZone=1
  }

  area.vec.w=st_area(wkpos_chull)
  if(w==1){
    area.vec=area.vec.w
  } else{
    area.vec=c(area.vec,area.vec.w)
  }
  
}

rl=rl[,1:9]
rl$x=st_coordinates(rl)[,1]
rl$y=st_coordinates(rl)[,1]
rl=st_drop_geometry(rl)

rl$cell=NA
for(i in 1:nrow(rl)){
  if(i%%1000==0){
  print(paste0(i," of ",nrow(rl)))
  }
  rl$cell[i]=which(centroids[,1]==rl$x[i]&centroids[,2]==rl$y[i])
}

area.vec.all=data.frame("week"=weeks,"area"=as.numeric(area.vec))
return(list("area"=area.vec.all,"dat"=rl[,1:9]))

}

zone_out_I=findZone_sim(rl,"I")
Izone<-zone_out_I$dat
Iarea<-zone_out_I$area


E_wk.cell.summary=
  Izone[Izone$inZone==1,] %>% 
  dplyr::group_by(sample_source,week,grid_80km) %>% 
  dplyr::summarise(startdt=min(DATE),
                   Nsampled=sum(number.of.animals),
                   Npos=sum(PCR_Result),
                   outbreak=first(region))





