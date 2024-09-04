
GetOutputs<-function(pop,BB,Incidence,Tculled,ICtrue,out,detectday,out.opts,input.opts){
  #print("Entering GetOutputs")
  
  #List of outputs created here:	
	#Tinc #sum of all exposures over simulation 
	#sum(Tculled)
	#idT #last day there is an infectious individual
	#Mspread #max spread of infection
	#IConDD #number of I, C, and E on detection day
	#ICatEnd #number of I,C,E on last day 
	#TincFromDD #sum of all exposures starting day after detection day
	#TincToDD #sum of all exposures up until detection day
	#DET #total number of detections

#to add:
  #toggle option: locs, dataframe of x/y locs for each sounder at each timestep, with SEIRCZ status
  #toggle option: cell numbers of all in zone at each timestep 
  #toggle options: POSlive, POSdead, POSlive_locs, POSdead_locs, but as dataframes
  #toggle options: num new infections at each step
  #toggle options: R0 at each step
  
#############################################################
  
##############################
### Summarize main outputs ###
##############################  

#Tinc, this just sums all of the exposures
Tinc=sum(Incidence)

#sum all culled
sumTculled=sum(Tculled)

if(any(ICtrue!=0)){
#Find last day there was an infectious individual
idT=which(ICtrue!=0)[length(which(ICtrue!=0))]
#IConDD #number of I, C, and E on detection day
IConDD=ICtrue[detectday]
#ICatEnd #number of I,C,E on last day
ICatEnd=ICtrue[idT]
} else{
  idT=1
  IConDD=0
  ICatEnd=0
}

#Find max spread of infection
Mspread<-max(out[,2])

#TincToDD #sum of all exposures up until detection day
TincToDD<-sum(Incidence[idT:detectday])

#TincFromDD #sum of all exposures starting day after detection day
TincFromDD<-sum(Incidence[detectday:idT])

#DET #total number of detections
#DET=sum(unlist(POSlive),unlist(POSdead))
DET=NA
#############################################################

############################
### Compile main outputs ###
############################

#send all main outputs to list
list.all=list("Tinc"=Tinc,
              "sumTculled"=sumTculled,
              "BB"=BB,
              "Mspread"=Mspread,
              "IConDD"=IConDD,
              "ICatEnd"=ICatEnd,
              "TincToDD"=TincToDD,
              "TincFromDD"=TincFromDD,
              "DET"=DET)

#############################################################

##################################
### Summarize optional outputs ###
##################################

#loc.list: is pop[,c(3,8:13)] for each timestep in list of length thyme
#cell num, S, E, I, R, C, Z
#convert loc.list into dataframe
if("sounderlocs"%in%out.opts){
  
  #print("Summarizing sounderlocs")
  #print(input.opts$loc.list)
  loc.list=input.opts$loc.list
for(i in 1:length(input.opts$loc.list)){
  #print(i)
  #print(length(loc.list[[i]]))
  #print(length(loc.list))
  if(length(loc.list[[i]]!=0)){
  #print(loc.list[[i]][,1])
  locs.i=as.data.frame(centroids[loc.list[[i]][,1],])
  #print(locs.i)
  colnames(locs.i)=c("x","y")
  locs.i$timestep=i
  locs.i$S=loc.list[[i]][,2]
  locs.i$E=loc.list[[i]][,3]
  locs.i$I=loc.list[[i]][,4]
  locs.i$R=loc.list[[i]][,5]
  locs.i$C=loc.list[[i]][,6]
  locs.i$Z=loc.list[[i]][,7]
  
  #print(locs.i)
  if(i==1){
    locs.df=locs.i
  } else{
    locs.df=rbind(locs.df,locs.i)
  }
  }
}
  #convert to integer for plotting locations with geom_tile
  #locs.df$x=as.integer(locs.df$x)
  #locs.df$y=as.integer(locs.df$y)
}


#Get POSlive, dead, and all locs in datafame format from list
if("alldetections"%in%out.opts){
  #print("Summarizing alldetections")
  
  #all should be same length, thyme
  POSlive=input.opts$POSlive
  POSlive_locs=input.opts$POSlive_locs
  POSdead=input.opts$POSdead
  POSdead_locs=input.opts$POSdead_locs
  
  detections=matrix(nrow=0,ncol=4)
  for(i in 1:length(POSlive)){
    live.detections.i=matrix(nrow=length(c(POSlive_locs[[i]])),ncol=4)
    live.detections.i[,1]=i
    live.detections.i[,2]=1 #code for live
    live.detections.i[,3]=POSlive[[i]]
    live.detections.i[,4]=c(POSlive_locs[[i]])
    dead.detections.i=matrix(nrow=length(c(POSdead_locs[[i]])),ncol=4)
    dead.detections.i[,1]=i
    dead.detections.i[,2]=0 #code for dead
    dead.detections.i[,3]=POSdead[[i]]
    dead.detections.i[,4]=c(POSdead_locs[[i]])
    
    detections.i=rbind(live.detections.i,dead.detections.i)
    detections=rbind(detections,detections.i)
    }
  
}


#Get Incidence and R0 vals summarized in data frame
if("incidence"%in%out.opts){
  #print("Summarizing incidence")
  
  #Timestep
  #Num new exposures from infected (incidence)
  #Num infected (prev timestep)
  #Rt (total, homog) at each time step
  inc.mat=input.opts$incidence
  inc.df<-as.data.frame(inc.mat)
  colnames(inc.df)=c("timestep","state","loc")
  inc.df[inc.df$state==10,2]<-"infected" #any infected that could have exposed susceptibles
  inc.df[inc.df$state==9,2]<-"exposed" #exposed this timestep
  inc.df[inc.df$state==12,2]<-"carcass" #exposed this timestep
  
  #inc.summary=inc.df %>% dplyr::group_by(timestep) %>% dplyr::summarise(num_inf=n(state))
}

#############################################################

###############################
### Append optional outputs ###
###############################

if("sounderlocs"%in%out.opts){
  templist=vector(mode="list",length=1)
  templist[[1]]=locs.df
  list.all=append(list.all,templist)
  names(list.all)[length(list.all)]="sounderlocs"
}

if("idzone"%in%out.opts){
  templist=vector(mode="list",length=1)
  templist[[1]]=input.opts$idzone.mat
  list.all=append(list.all,templist)
  names(list.all)[length(list.all)]="idzone"
}

if("alldetections"%in%out.opts){
  templist=vector(mode="list",length=1)
  templist[[1]]=detections
  list.all=append(list.all,templist)
  names(list.all)[length(list.all)]="alldetections"
}

if("incidence"%in%out.opts){
  templist=vector(mode="list",length=1)
  templist[[1]]=inc.df
  list.all=append(list.all,templist)
  names(list.all)[length(list.all)]="incidence"
}

#############################################################

###########################
### Return final output ###
###########################

return(list.all)

}
