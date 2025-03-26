
#############################
########## Purpose ##########
#############################
#Series of functions to sample pigs in given summarized sounderlocs data

###############################
########## Functions ##########
###############################

#SampleIndividuals: 
  #outputs matrix with actual freq and apparent freq
  #helper function used in ExtApparentPrev
#SEIRCZ- data frame, single timestep summary of all states in spatial area
#N- integer, number of samples to draw
#DetP- double, probability that infected individuals will be positive
#timestep- integer, current timestep
#prev- boolean, whether to get prevalence or just raw counts
SampleIndividuals <-function(SEIRCZ,N,DetP,timestep,prev){
  
  S<-rep("S",SEIRCZ[1,which(colnames(SEIRCZ)=="S")]) 
  E<-rep("S",SEIRCZ[1,which(colnames(SEIRCZ)=="E")]) #Would E have antibodies? assign as S for now
  I<-rep("I",SEIRCZ[1,which(colnames(SEIRCZ)=="I")]) 
  R<-rep("S",SEIRCZ[1,which(colnames(SEIRCZ)=="R")]) #Would R have antibodies? assign as S for now
  C<-rep("C",SEIRCZ[1,which(colnames(SEIRCZ)=="C")])
  Z<-rep("Z",SEIRCZ[1,which(colnames(SEIRCZ)=="Z")])
  
  allpig=c(S,E,I,R)
  samples=sample(allpig,min(N,length(allpig)))
  counts.r=plyr::count(samples)
  
  allcarcass=c(C,Z)
  csamples=sample(allcarcass,min(N,length(allcarcass)))
  counts.c=plyr::count(csamples)
  
  if(nrow(counts.r)>0){
    counts.r$appar<-0
    
    if(length(counts.r[counts.r$x=="I",]$freq)>0){
      if(length(counts.r[counts.r$x=="S",]$freq)==0){
        counts.r=rbind(counts.r,data.frame(x="S","freq"=0,"appar"=0))
      }
      num.appar.inf<-sum(rbinom(1,counts.r[counts.r$x=="I",]$freq,prob=c(DetP,(1-DetP))))
      num.add.sus=counts.r[counts.r$x=="I",]$freq-num.appar.inf
      counts.r[counts.r$x=="S",]$appar=counts.r[counts.r$x=="S",]$freq+num.add.sus
      counts.r[counts.r$x=="I",]$appar=num.appar.inf
    } else{
      counts.r$appar=counts.r$freq
      counts.r=rbind(counts.r,data.frame(x="I","freq"=0,"appar"=0))
    }
    
  } else{
    counts.r=data.frame(x="S",freq=0,appar=0)
    counts.r=rbind(counts.r,data.frame(x="I","freq"=0,"appar"=0))
  }
  
  if(nrow(counts.c)>0){
    counts.c$appar=NA
    if(length(counts.c[counts.c$x=="C",]$freq)>0){
      if(length(counts.c[counts.c$x=="Z",]$freq)==0){
        counts.c=rbind(counts.c,data.frame(x="Z","freq"=0,"appar"=0))
        
      }
      
      num.appar.inf<-sum(rbinom(1,counts.c[counts.c$x=="C",]$freq,prob=c(DetP,(1-DetP))))
      num.add.sus=counts.c[counts.c$x=="C",]$freq-num.appar.inf
      counts.c[counts.c$x=="Z",]$appar=counts.c[counts.c$x=="Z",]$freq+num.add.sus
      counts.c[counts.c$x=="C",]$appar=num.appar.inf
      
      } else{
      counts.c$appar=counts.c$freq
      counts.c=rbind(counts.c,data.frame(x="C","freq"=0,"appar"=0))
    }
    
  } else{
    counts.c=data.frame(x="Z",freq=0,appar=0)
    counts.c=rbind(counts.c,data.frame(x="C","freq"=0,"appar"=0))
  }
  
  counts.out=rbind(counts.r,counts.c)
  
  counts.out=as.data.frame(pivot_wider(counts.out,names_from=x,values_from=c(freq,appar)))
  counts.out[is.na(counts.out)]<-0
  
  
  if(prev){
    prev_out=data.frame(matrix(nrow=1,ncol=4))
    colnames(prev_out)=c("Ipt","Ipa","Cpt","Cpa")
    prev_out$Ipt=counts.out$freq_I/(counts.out$freq_I+counts.out$freq_S)
    prev_out$Ipa=counts.out$appar_I/(counts.out$appar_I+counts.out$appar_S)
    prev_out$Cpt=counts.out$freq_C/(counts.out$freq_C+counts.out$freq_Z)
    prev_out$Cpa=counts.out$appar_C/(counts.out$appar_C+counts.out$appar_Z)
    
    #if(!is.na(prev_out$Cpa>1)){
    #if(prev_out$Cpa>1){
    #  print(counts.out)
    #}
    #}
    
    #if(!is.na(prev_out$Ipa>1)){
    #  if(prev_out$Ipa>1){
    #    print(counts.out)
    #  }
    #}
    
    prev_out[is.na(prev_out)]=0
    
    if(!missing(timestep)){
      prev_out$timestep=timestep
    }
    
    return(prev_out)
    
  } else{
    if(!missing(timestep)){
      counts.out$timestep=timestep
    }
    return(as.matrix(counts.out))
    
  }
  
  
  
  
}

#doQuartileSummary: helper function used in ExtApparentPrev, used in cell-level summary option
  #gets quantile summaries (median, 0.25, 0.75) for a given matrix, not including the cols.exclude indices
#mat-input matrix to summarize
#cols.exclude-vector of column id numbers to exclude from summary. just gets first value.
doQuartileSummary<-function(mat,cols.exclude){
  cols=1:ncol(mat)
  cols=cols[-cols.exclude]
  strings=colnames(mat)[cols]
  
  #do summaries
  #sum.mat=matrix(nrow=1,ncol=3*length(cols))
  for(c in 1:length(cols)){
    sub.mat=matrix(nrow=1,ncol=3)
    colnames(sub.mat)=c(paste(strings[c],"med",sep="_"),
                        paste(strings[c],"q25",sep="_"),
                        paste(strings[c],"q75",sep="_"))
    sub.mat[1,1]=median(mat[,cols[c]])
    sub.mat[1,2]=quantile(mat[,cols[c]],0.25)
    sub.mat[1,3]=quantile(mat[,cols[c]],0.75)
    
    if(c==1){
      out.mat=sub.mat
    } else{
      out.mat=cbind(out.mat,sub.mat)
    }
  }
  
  
  out.mat=cbind(out.mat,mat[1,cols.exclude])
  colnames(out.mat)[ncol(out.mat)]=colnames(mat)[cols.exclude]
  
  return(out.mat)
}

#ExtApparentPrev- gets apparent prevalence across each timestep from summarized count data of each state
  #by sampling N individuals per timestep, and calculating infection prevalence from samples
  #used in sounderlocsSummarize.R to get apparent prev output
#SEIRCZ.rep-data frame summarized in sounderlocsSummarize.R
#DetP-double, probability that infected individuals will show positive (i.e., sensitivity)
#N-integer, number of individuals to sample each timestep
#rep-replicate number
ExtApparentPrev<-function(SEIRCZ.rep,DetP,N,rep){
  require(plyr)
  
  #inputs from sounderlocsSummarize
  #SEIRCZ.rep
  #SEIRCZ.rep_zone
  #SEIRCZ.rep_zone_cells
  #View(SEIRCZ.rep_zone_cells)
  
  timesteps=unique(SEIRCZ.rep$timestep)
  is_cells="cell"%in%colnames(SEIRCZ.rep)
  
  for(t in 1:length(timesteps)){
    
    SEIRCZ.t=SEIRCZ.rep[SEIRCZ.rep$timestep==t,,drop=FALSE]
    
    if(nrow(SEIRCZ.t)>0&is_cells){
      #sample 13 rows (cells) for given timestep
      cells_sampled=sample(1:nrow(SEIRCZ.t),min(13,nrow(SEIRCZ.t)))
      
      for(c in 1:(length(cells_sampled))){
        counts.t.c=SampleIndividuals(SEIRCZ.t[cells_sampled[c],],1,DetP,t,TRUE)
        #throws error here for sample where no S individuals-- only 1 Z
        
        if(c==1){
          counts.t=counts.t.c
        } else{
          counts.t=rbind(counts.t,counts.t.c)
        }
        
        
        
        
      }
      
      counts.t=doQuartileSummary(counts.t,which(colnames(counts.t)=="timestep"))
      
    } else{
      if(nrow(SEIRCZ.t)>0){
      
      counts.t=SampleIndividuals(SEIRCZ.t,N,DetP,t,TRUE)
      #counts.t=doQuartileSummary(counts.t,which(colnames(counts.t)=="timestep"))
      } else{
        if(is_cells){
          counts.t=matrix(0,nrow=1,ncol=13)
          colnames(counts.t)=c("Ipt_med","Ipt_q25","Ipt_q75","Ipa_med","Ipa_q25",
                               "Ipa_q75","Cpt_med","Cpt_q25","Cpt_q75","Cpa_med",
                               "Cpa_q25","Cpa_q75","timestep")
          counts.t[,13]=t
        } else{
        counts.t=matrix(0,nrow=1,ncol=5)
        colnames(counts.t)=c("Ipt","Ipa","Cpt","Cpa","timestep")
        counts.t[,5]=t
        }
      }
    }
    
    if(t==1){
      counts=counts.t
    } else{
      counts=rbind(counts,counts.t)
    }
  }
  
  counts=cbind(rep,counts)
  
  return(counts)
  
}




