##The purpose of this script is to run a single rep of the ASF control optimization model
#out.opts=c("sounderlocs","idzone","alldetections","incidence")

#SimulateOneRun<-function(Pcr,Pir,Pbd,death,F1,F2_int,F2_B,F2i_int,F2i_B,B1,B2,thyme,cells,N0,K,detectday,Rad,Intensity,alphaC,shift,centroids,cullstyle,inc,ss,gridlen,midpoint,pop,out.opts,grid.opts,rep,DetP){
SimulateOneRun<-function(outputs,pop,centroids,grid,parameters,cpp_functions,K){
require(dplyr)

	for(i in 1:length(cpp_functions)){
		print(paste0("sourcing ",cpp_functions[[i]]))
		Rcpp::sourceCpp(cpp_functions[[i]])
		}

############################################################
######## Release parameters to function environment ######## 
############################################################
list2env(parameters, .GlobalEnv)
	
###########################################
######## Initialize Output Objects ######## 
###########################################
list2env(outputs, .GlobalEnv)
	
#track first infection in Incidence matrix
Incidence[1]<-num_inf_0
	
##################################
######## Start simulation ######## 
##################################
  
for(i in 1:thyme){

print(paste0("timestep: ",i))
if(any(pop[,9,drop=FALSE]!=0|pop[,10,drop=FALSE]!=0|pop[,12,drop=FALSE]!=0)){

if("sounderlocs"%in%out.opts){
  loc.list[[i]]=pop[,c(3,8:13)]
  }

#####################################
######## Track I/C locations ######## 
#####################################
	print("Track I/C locs")

if(nrow(pop[pop[,10]>0,,drop=FALSE])>0){
Isums[i]<-nrow(pop[pop[,10]>0,,drop=FALSE])
} else{Isums[i]=0}

if(nrow(pop[pop[,12]>0,,drop=FALSE])>0){
Csums[i]<-nrow(pop[pop[,12]>0,,drop=FALSE])
} else{Csums[i]=0}

I_locs[[i]]<-pop[pop[,10]>0,3]
C_locs[[i]]<-pop[pop[,12]>0,3]



##########################
######## Movement ######## 
##########################


pop<-FastMovement(pop,centroids,shape,rate,inc,mv_pref)


###############################
######## State Changes ########
###############################
#births, natural deaths, disease state changes (exposure, infection, recovery, death), carcass decay

## many of these parameters are missing; looks like they are leftover from the pre-targets era
st.list<-StateChanges(pop,centroids,nrow(centroids),parameters, i)
# st.list<-StateChanges(pop,centroids,nrow(centroids),Pbd,B1,B2,F1,F2_int,F2_B,F2i_int,F2i_B,K,death,Pcr,Pir,Incidence,BB,i)

Eep_mat[i,]=st.list$Eep
Sdpb_mat[i,]=st.list$Sdpb
Sdpd_mat[i,]=st.list$Sdpd
Rep_mat[i,]=st.list$Rep
Cep_mat[i,]=st.list$Cep
Rdpd_mat[i,]=st.list$Rdpd
Ccd_mat[i,]=st.list$Ccd
Zcd_mat[i,]=st.list$Zcd
Iep_mat[i,]=st.list$Iep

if("incidence"%in%out.opts){
## is there a reason this is counting from the pre-state changes list instead of st.list[[1]]?
  inf.locs=rep(pop[(pop[,10]>0),3], pop[(pop[,10]>0),10]) #locs infected
#   inf.locs=pop[,3][pop[,10]>0] #locs infected
#   inf.num=sum(pop[(pop[,10]>0),10]) #num infected
  inf.num=sum(pop[,10]) #num infected
#   exp.locs=rep(pop[st.list[[4]]>0,3],st.list[[4]][st.list[[4]]>0]) #locs exposed
  ## is that what we want? pop[st.list[[4]]>0,3] gives all the rows because the row arg = TRUE
  ## so we end up repeating the entire pop's cell numbers a zillion times?
  ## I think it might be old, because it seems to expect a vector for st.list[[4]] where it is now a sum()
  exp.locs=rep(st.list[[1]][st.list[[1]][,9]>0,3], st.list[[1]][,9][st.list[[1]][,9]>0]) #locs exposed
  ## st.list[[4]] is not what we want (not a vector), but pop and st.list[[1]] have different number of rows, so have to use st.list[[1]] instead of pop I guess
#   exp.locs=st.list[[1]][,3][st.list[[1]][,9]>0] #locs exposed (E)
#   exp.num=sum(st.list[[4]][st.list[[4]]>0]) #num exposed
  exp.num=sum(st.list[[1]][,9]) #num exposed
#   c.locs=rep(pop[(pop[,12]>0),3],pop[(pop[,12]>0),12]) #locs infected
#   c.locs=st.list[[1]][,3][st.list[[1]]][,12]>0] #locs infected (C)
  c.locs=rep(st.list[[1]][st.list[[1]][,12]>0,3], st.list[[1]][,12][st.list[[1]][,12]>0]) #locs exposed
#   c.num=sum(pop[(pop[,12]>0),12]) #num infected
  c.num= sum(st.list[[1]][,12])#num infected

  if(inf.num+exp.num>0){
    inc.mat.i=matrix(nrow=(inf.num+exp.num+c.num),ncol=3)
    inc.mat.i[,1]=i

    if(length(inf.locs)!=0){
      inc.mat.i[1:inf.num,3]=inf.locs
      inc.mat.i[1:inf.num,2]=10 #code 10 for inf, same as pop colnum
    }

    if(length(exp.locs)!=0){

      inc.mat.i[(inf.num+1):((inf.num)+exp.num),3]=exp.locs
      inc.mat.i[(inf.num+1):((inf.num)+exp.num),2]=9 #code 9 for exp, same as pop colnum
    }

    if(length(c.locs)!=0){
      inc.mat.i[(inf.num+exp.num+1):nrow(inc.mat.i),3]=c.locs
      inc.mat.i[(inf.num+exp.num+1):nrow(inc.mat.i),2]=12 #code 12 for exp, same as pop colnum
    }

    if(i==1){
      inc.mat=inc.mat.i
    } else{
      inc.mat=rbind(inc.mat, inc.mat.i)
    }
  } else{
    if(i==1){
      inc.mat=matrix(nrow=0,ncol=3)
    }
  }

}

pop<-st.list[[1]]
Incidence<-st.list[[2]]
BB<-st.list[[3]]

###**start on day of first detection
###################################
######## Initiate Response ########
###################################

# If sampling turned on, run surveillance scheme based on user's input
# for the current week
if(sample == 1){ ## I think this is supposed to be the 'sample' parameter?
# if(sampling == 1){
  sample.design <- PrepSurveillance(sample) ## this is the only place sample.design is defined (sample doesn't do anything,but it's in the function definition)
  surv.list<-Surveillance(pop,i,sample.design,grid.list,inc,POSlive,POSdead,POSlive_locs,POSdead_locs,pigs_sampled_timestep) # Madison
  pop=surv.list[[1]]
  POSlive=surv.list[[2]]
  POSdead=surv.list[[3]]
  POSlive_locs=surv.list[[4]]
  POSdead_locs=surv.list[[5]]
  pigs_sampled_timestep=surv.list[[6]]
}

# If sampling turned off and it's detect day based on user input,
# run FirstDetect because there are infected pigs to detect,  and Rad>0
if(sample != 1 & i==detectday&sum(pop[, c(9, 10, 12)])>0&Rad>0){
# if(sampling != 1 & i==detectday&sum(pop[, c(9, 10, 12)])>0&Rad>0){
  fd.list<-FirstDetect(pop, i, POSlive, POSdead, POSlive_locs, POSdead_locs)
  pop=fd.list[[1]]
  POSlive=fd.list[[2]]
  POSdead=fd.list[[3]]
  POSlive_locs=fd.list[[4]]
  POSdead_locs=fd.list[[5]]
}



###**start day after day of first detection
#######################################
######## Initiate Culling Zone ######## 
#######################################

#if it is at least day after detect day,  and Rad>0
if(sample != 1 & i > detectday & Rad > 0) {
# if(sampling != 1 & i > detectday & Rad > 0) {

	#new detections from last step, bc day lag 
	#(either from initial detection or last culling period)
	#get locations in grid for detections
	idNEW=c(POSlive_locs[[i-1]],POSdead_locs[[i-1]])
	#remove NA/0 (may get NAs/zeroes if no live/dead detected)
	idNEW<-idNEW[idNEW>0&!is.na(idNEW)]
	idZONE_t=idZONE
	#if there were detections in previous time steps, only get newly detected infected grid cells
	#"infected grid cell"=grid cell where there was an infected pig or carcass
	if(length(idZONE_t[,1])>0){
	uniqueidNEW<-which(!(idNEW %in% idZONE_t))
	idNEW<-idNEW[uniqueidNEW]
	} else{idNEW=idNEW}

	#Culling process
	output.list<-CullingOneRun(pop,idNEW,idZONE,Intensity,alphaC,centroids,Rad,inc,i,detected,POSlive,POSdead,POSlive_locs,POSdead_locs,NEGlive,NEGdead,DetP,cullstyle)
# 	output.list<-CullingOneRun(pop,idNEW,idZONE,Intensity,alphaC,centroids,Rad,inc,i,detected,POSlive,POSdead,POSlive_locs,POSdead_locs,NEGlive,NEGdead,DetP)
	POSlive[[i]]<-output.list[[1]]
	POSdead[[i]]<-output.list[[2]]
	POSlive_locs[[i]]<-output.list[[3]]
	POSdead_locs[[i]]<-output.list[[4]]
	NEGlive[[i]]<-output.list[[5]]
	NEGdead[[i]]<-output.list[[6]]
	idZONE<-output.list[[7]]
	removalcells[[i]]<-output.list[[8]]
	culled<-output.list[[9]]
	ZONEkm2[i,]<-output.list[[10]]
	pop<-output.list[[11]]
	Ct[i,1]<-output.list[[12]]
	#Total number culled at each timestep
	Tculled[i]=culled
	
	#compile optional outputs
	if("idzone"%in%out.opts){

	  #get list index
	  idz=(detectday-i)
	  if(idz==1){
	    idzone.mat.idz=idzone[,2]
	    idzone.mat=cbind(idzone.mat.idz,rep(i,length=length(idzone.mat.idz)))
	    colnames(idzone.mat)=c("cell","timestep")
	  } else{
	    #only store new locations
	    idzone.mat.idz=idzone[,2](which(!(idzone[,2])%in%idzone.mat[,1]))
	    idzone.mat.idz=cbind(idzone.mat.idz,rep(i,length=length(idzone.mat.idz)))
	    colnames(idzone.mat.idz)=c("cell","timestep")
	  }
	  
	}
	
	
} else{ #if greater than detectday closing bracket
	Ct[i,1]=0
	
	} #else if not greater than detect day


#############################
####Track true spatial spread
#############################
#if any infected individuals
if(nrow(pop[pop[,9,drop=FALSE]>0|pop[,10,drop=FALSE]>0|pop[,12,drop=FALSE]>0,,drop=FALSE])>0){
  out[i,]<-areaOfinfection(pop,centroids,inc)
} else{out[i,]=c(0,0,0)}

#############################
####Summarize infections
#############################
#sum all infectious cases (I,C,E) at each timestep
#ICtrue = sum(I + C,2); sum of all infectious cases over time
if(i==detectday){
ICtrue[i]<-(sum(colSums(pop)[c(9,10,12)])+1) #account for having removed that first detected
} else{
	ICtrue[i]<-sum(colSums(pop)[c(9,10,12)])
}

#} #for manual testing of loop

####Update population matrix
#Remove rows in pop with 0 pigs
pigcols=c(1,8:13)
pop=pop[which(rowSums(pop[,pigcols,drop=FALSE])!=0),,drop=FALSE]

}}
# } else{
#   #print("Exiting loop, no infections")
#   } #if any infected closing bracket/else
# 	} #for timestep closing bracket

#############################
#############################
if(length(out.opts)>0){
input.opts=vector(mode="list",length=1)
input.opts[[1]]=out.opts
names(input.opts)[1]="out.opts"
if("sounderlocs"%in%out.opts){
  templist=vector(mode="list",length=1)
  templist[[1]]=loc.list
  input.opts=append(input.opts,templist)
  names(input.opts)[length(input.opts)]="loc.list"
}

if("idzone"%in%out.opts){
  templist=vector(mode="list",length=1)
  templist[[1]]=idzone.mat
  input.opts=append(input.opts,templist)
  names(input.opts)[length(input.opts)]="idzone.mat"
}

if("alldetections"%in%out.opts){
  templist = list(POSlive)  # directly create a list with POSlive
  input.opts = append(input.opts, templist)
  names(input.opts)[length(input.opts)] = "POSlive"
  
  templist = list(POSdead)
  input.opts = append(input.opts, templist)
  names(input.opts)[length(input.opts)] = "POSdead"
  
  templist = list(POSlive_locs)
  input.opts = append(input.opts, templist)
  names(input.opts)[length(input.opts)] = "POSlive_locs"
  
  templist = list(POSdead_locs)
  input.opts = append(input.opts, templist)
  names(input.opts)[length(input.opts)] = "POSdead_locs"
  
#   if(sampling == 1){
  if(sample == 1){
    templist = list(pigs_sampled_timestep)  # directly create a list with pigs_sampled_timestep
    input.opts=append(input.opts,templist)
    names(input.opts)[length(input.opts)]="pigs_sampled_timestep"
  }
}

if("incidence"%in%out.opts){
  templist=vector(mode="list",length=1)
  templist[[1]]=inc.mat
  input.opts=append(input.opts,templist)
  names(input.opts)[length(input.opts)]="incidence"
  
  
}
}

list.all<-GetOutputs(pop,centroids,BB,Incidence,Tculled,ICtrue,out,detectday,Ct,out.opts,input.opts)


return(list.all)

	} #function closing bracket
