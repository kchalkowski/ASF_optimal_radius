#Modeling state data for ASF simulation contact parameters

#################################################
## Below is a summary of the variables and other info for the files needed to Model State Data

#Order readout from list.files:
#[1] "direct_marec_all1.csv" 
#[2] "direct_srel_all1.csv"  
#[3] "indirect_marec_all.csv"
#[4] "indirect_srel_all.csv" 

#File use:
#Pig contact     = "direct"
#Carcass contact = "indirect"
#FL State        = "marec"
#SC State        = "srel"

#Variables inside each data set: 
  
#Pig Contact Data (both direct) (direct_marec_all1.csv, direct_srel_all1.csv)
#- "id1"          = id of pig # 1      
#- "id2"          = id of pig # 2
#- "datetime"     = date and time when contact happened, date
#- "numPair"      = Not needed in script-- not sure what this is
#- "numbers"      = Number of contacts for each date, for given id pairing
#- "HRdist"       = Home range distance in meters, numeric
#- "HRover"       = Not needed in script, not sure what this is
#- "sex1"         = sex of pig # 1, male/female
#- "sex2"         = sex of pig # 2, male/female
#- "dynamicgroup" = whether the contact represents a within-group or between-group contact
    #1- between-group contact
    #0- within-group contact
  
#Carcass Contact Data (FL - marec)
#- "id"          = ID of pig contacting carcass, number-factor
#- "baitid"      =id of carcass (bait?), factor
#- "datetime"    = date and time when contact happened, date
#- "numPair"     = Not needed in script-- not sure what this is
#- "pairID"      = unique identifier for each pairing for each day
#- "numbers"     = Number of contacts for each pairing/day
#- "DistBait"    = Distance from bait in meters, numeric
#- "type...8"    = bait type, character **
#- "type...9"    = bait type, numeric identifier **
#- "...10"       = type of attractant menu **

#** Carcass contact data is currently underway
#** For now, just using bait contact rates as a proxy for carcass contact

#Carcass Contact Data (SC - srel)
#- "id1"         = ID of pig contacting bait, number-factor    
#- "id2"         = unique ID of bait, number-factor
#- "datetime"    = Date and time when contact happened, date
#- "numbers"     = Number of contacts between pig/bait pairing for that date, numeric
#- "DistBait"    = Distance from bait in meters, numeric
#- "sex1"        = sex of pig contacting bait, male/female
#- "numPair"     = Not needed in script-- not sure what this is
#- "type"        = type of attractant, character-factor

######

#1-read in data
contact.dat=list.files(paste0(home,"/Input"),full.names=TRUE)
contact.dat=contact.dat[grep("csv",contact.dat)]

pig.contact=contact.dat[grep("/direct",contact.dat)]
carcass.contact=contact.dat[grep("/indirect",contact.dat)]
contact.names=c("FL","SC")

F1.list=vector(mode="list",length=length(contact.names))
F2.list=vector(mode="list",length=length(contact.names))
names(F1.list)=contact.names
names(F2.list)=contact.names

#do direct files first, pig contact
for(i in 1:length(pig.contact)){
  
  #read data
  contact=read.csv(pig.contact[i])

  #format datetime
  if(i==1){
    contact$datetime=as.POSIXct(contact$datetime,format="%Y-%m-%d",tz="UTC")
  } else{
    contact$datetime=as.POSIXct(contact$datetime,format="%m/%d/%Y",tz="UTC")
  }
  contact$numdat=as.numeric(contact$datetime) #in seconds
  contact$numdat=contact$numdat/60/60/24
  contact$numdat=contact$numdat-min(contact$numdat)+1 #get all in reference to earliest date
  
  ####reflect same matlab work
  Z=matrix(nrow=nrow(contact),ncol=2)
  Z[,1]=contact$numbers
  Z[,2]=contact$HRdist
  
  G=matrix(nrow=nrow(contact),ncol=1)
  G[,1]=contact$dynamicgroup

  ID=matrix(nrow=nrow(contact),ncol=2)
  ID[,1]=contact$id1
  ID[,2]=contact$id2
  
  date=contact$numdat
  
  uid1=unique(contact$id1)
  uid2=unique(contact$id2)
  
  steps=seq(0,max(date),7)
  
  #initialize weekly matrix
  weekly=matrix(nrow=0,ncol=5)
  for(kk in 2:length(steps)){
    #kk=2
    #which date is less than or equal to 2, AND greater than 1
    #ask Kim-- maybe didn't want to include week 1 for some reason
    index=which(date<=steps[kk]&date>steps[(kk-1)])
    ID.kk=ID[index,]
    Z.kk=Z[index,]
    G.kk=G[index,]
    
    uid=unique(ID.kk)
    #for each unique pair in week kk
    for(ii in 1:nrow(uid)){
      index2=which(paste(ID.kk[,1],ID.kk[,2],sep="_")==paste(uid[ii,1],uid[ii,2],sep="_"))
      
      ID.kk.ii=t(as.matrix(ID.kk[index2,]))
      Z.kk.ii=c(NA,NA)
      Z.kk.ii[1]=sum(Z.kk[index2,1]) #sum of contacts
      Z.kk.ii[2]=mean(Z.kk[index2,2]) #mean of hr dist
      G.kk.ii=G.kk[index2]
      #just take first value for dynamic group, even tho changes
      weekly=rbind(weekly,c(ID.kk.ii[1,1],ID.kk.ii[1,2],Z.kk.ii,G.kk.ii[1]))
    
      }
    
  }
  
  #Get within-group contact
  wt=weekly[which(weekly[,5]==0),]
  wt1=wt[which(wt[,3]>0),3]
  wt0=wt[which(wt[,3]==0),3]
  W=length(wt1)/(length(wt0)+length(wt1))
  
  #get between-group contact
  bt=weekly[which(weekly[,5]>0),]
  bt1=bt[which(bt[,3]>0),3]
  bt0=bt[which(bt[,3]==0),3]
  
  X=bt[,4]/1000
  Y=bt[,3]
  Y[Y>0]=1
  M=glm(Y~X,family="binomial")

  #assign to lists
  F1.list[[i]]=W
  F2.list[[i]]=M


  } #FL/SC loop


#FL=
#results spot on for both types of contact

#SC=
#W should be 0.9375, is 0.9375
#M model should be y=0.63955-1.16x, is (after fixed date formatting)


########################################
##### Carcass contact modeling #########
########################################

F2i.list=vector(mode="list",length=length(contact.names))
names(F2i.list)=contact.names
for(i in 1:length(carcass.contact)){
  contact=read.csv(carcass.contact[i])
  
  if(i==1){
    #for FL pigs, subset just molasses ones
    contact=contact[contact$type=="molasses",]
  } 
  
  colnames(contact)[2]<-"baitid"
  
  #format datetime
  contact$datetime=as.POSIXct(contact$datetime,format="%m/%d/%Y",tz="UTC")
  contact$numdat=as.numeric(contact$datetime) #in seconds
  contact$numdat=contact$numdat/60/60/24
  contact$numdat=contact$numdat-min(contact$numdat)+1 #get all in reference to earliest date
  
  ####reflect same matlab work
  Z=matrix(nrow=nrow(contact),ncol=2)
  Z[,1]=contact$numbers
  Z[,2]=contact$DistBait
  
  #G=matrix(nrow=nrow(contact),ncol=1)
  #G[,1]=contact$DistBait
  
  ID=matrix(nrow=nrow(contact),ncol=2)
  ID[,1]=as.numeric(contact$id)
  ID[,2]=as.numeric(contact$baitid)
  
  date=contact$numdat
  
  uid1=unique(contact$id1)
  uid2=unique(contact$id2)
  
  steps=seq(0,max(date),7)
  
  #initialize weekly matrix
  weekly=matrix(nrow=0,ncol=4)
  for(kk in 2:length(steps)){
    #kk=2
    #which date is less than or equal to 2, AND greater than 1
    #ask Kim-- maybe didn't want to include week 1 for some reason
    index=which(date<=steps[kk]&date>steps[(kk-1)])
    if(length(index)>0){
    ID.kk=ID[index,]
    Z.kk=Z[index,]
    G.kk=G[index,]

    uid=unique(ID.kk)
    #for each unique pair in week kk
    for(ii in 1:nrow(uid)){
      index2=which(paste(ID.kk[,1],ID.kk[,2],sep="_")==paste(uid[ii,1],uid[ii,2],sep="_"))
      
      ID.kk.ii=ID.kk[index2,,drop=FALSE]
      Z.kk.ii=c(NA,NA)
      Z.kk.ii[1]=sum(Z.kk[index2,1])
      Z.kk.ii[2]=mean(Z.kk[index2,2])
      #just take first value for dynamic group, even tho changes
      weekly=rbind(weekly,c(ID.kk.ii[1,1],ID.kk.ii[1,2],Z.kk.ii))
      
    }
    } #if index length > 0 closing bracket
  }
  
  Z2=weekly[,3:4]
  Bt=Z2
  Bt1=Bt[Bt[,1]>0,]  
  Bt0=Bt[Bt[,1]==0,]
  
  X=weekly[,4]/1000
  Y=Bt[,1]
  Y[Y>0]=1
  
  M=glm(Y~X,family="binomial")
  #FL: 2.322-6.372, should be 2.3215-6.3723
  #SC: 9.820-4.107, is same as Matlab version
  F2i.list[[i]]=M

}



