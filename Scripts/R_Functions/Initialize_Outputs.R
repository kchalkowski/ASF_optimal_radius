
###########################################
######## Initialize Output Objects ######## 
###########################################

Initialize_Outputs<-function(parameters){

#pull needed items from parameters 	
thyme=parameters$thyme	
out.opts=parameters$out.opts

Nall=matrix(nrow=thyme) #track total abundance
BB=matrix(nrow=thyme) #track births

POSlive=as.list(rep(0,thyme)) #Positive cases observed and removed from landscape
POSdead=as.list(rep(0,thyme))#Positive carcasses observed and removed from landscape
NEGlive=as.list(rep(0,thyme)) #Negative tests of detected carcasses that are removed from landscape
NEGdead=as.list(rep(0,thyme)) #Negative tests of carcasses that are removed from landscape
pigs_sampled_timestep=as.list(rep(0,thyme))  # Initialize an empty list to hold the number of pigs sampled at each timestep
Ct=matrix(nrow=thyme,ncol=1) #effective removal rate

POSlive_locs<-as.list(rep(0,thyme))
POSdead_locs<-as.list(rep(0,thyme))

idZONE=matrix(nrow=1,ncol=3) #grid cell ids that had a positive detection, grid cell ids that are within the zone, distance
#idZONE<-as.list(rep(NA,thyme)) #comment out of list mar 28
Tculled=matrix(0,nrow=thyme) #total number culled at each time step
ZONEkm2=matrix(0,nrow=thyme) 
Carea=matrix(0,nrow=thyme) #area of culling zone at each time step
Spread=matrix(0,nrow=thyme, ncol=3) #number of infectious individuals, area of infection, max distance between any two cases
Incidence=matrix(0,nrow=thyme) #store new cases for each time step
I_locs=vector("list",thyme)
C_locs=vector("list",thyme)
removalcells=vector("list",thyme)
I_locs[1:thyme]<-0
C_locs[1:thyme]<-0
Isums<-matrix(0,nrow=thyme)
Csums<-matrix(0,nrow=thyme)
out=matrix(c(0,0,0),nrow=thyme,ncol=3)
ICtrue=matrix(0,nrow=thyme,ncol=1)

#State change outputs when needed
#list(pop,Incidence,BB,"Eep"=Eep,"Sdpb"=Sdpb,"Sdpd"=Sdpd,"Iep"=Iep,"Edp"=Edpd,"Rep"=Rep,"Cep"=Cep,"Rdp"=Rdpd,"Ccd"=Ccd,"Zcd"=Zcd)
Eep_mat=matrix(0,nrow=thyme,ncol=1)
Sdpb_mat=matrix(0,nrow=thyme,ncol=1)
Sdpd_mat=matrix(0,nrow=thyme,ncol=1)
Iep_mat=matrix(0,nrow=thyme,ncol=1)
Rep_mat=matrix(0,nrow=thyme,ncol=1)
Cep_mat=matrix(0,nrow=thyme,ncol=1)
Rdpd_mat=matrix(0,nrow=thyme,ncol=1)
Ccd_mat=matrix(0,nrow=thyme,ncol=1)
Zcd_mat=matrix(0,nrow=thyme,ncol=1)

##Initialize out.opts objects as needed
if("sounderlocs"%in%out.opts){
  #Initialize list to track locations
  loc.list=vector(mode="list",length=thyme)
}

if("idzone"%in%out.opts){
  #Initialize list of idzones
  idzone.mat=matrix(nrow=0,ncol=2)
}


values = as.list(environment())
values = values[setdiff(names(values), names(formals()))]
values = values[setdiff(names(values),names(parameters))]

return(values)

}




