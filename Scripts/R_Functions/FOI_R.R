FOI_R<-function(pop,centroids,cells,B1,B2,F1,F2_int,F2_B,F2i_int,F2i_B){
  #Create I/C matrices. Just want vector with locations of I/C, with nums of I/C in right place
  I=matrix(0,nrow=cells,ncol=1)
  C=matrix(0,nrow=cells,ncol=1)
  I[pop[which(pop[,10]!=0),3],1]=pop[which(pop[,10]!=0),10]
  C[pop[which(pop[,12]!=0),3],1]=pop[which(pop[,12]!=0),12]
  
  Pse=matrix(0,nrow=cells,ncol=1)
  temp=I+C
  id=which(temp>0)
  W=matrix(0,nrow=cells,ncol=2)

  W[,1]=F1*I
  W[,2]=F1*C
  
  #B=zeros(length(id),cells)
  
  if(length(id)==1){
  X1=matrix(nrow=length(id),ncol=1)
  Y1=matrix(nrow=length(id),ncol=1)
  
  X1=centroids[id,1]
  Y1=centroids[id,2]
  
  dist=sqrt((centroids[,1]-X1)^2+(centroids[,2]-Y1)^2)
  dist=as.data.frame(dist)
  colnames(dist)="X"
  #prob=predict(F2,dist,type="response")
  prob=F2_int+F2_B*dist
  #probi=predict(F2i,dist,type="response")
  probi=F2i_int+F2i_B*dist
  prob[dist==0]=0
  probi[dist==0]=0
  
  B=B1*I[id[1]]*prob + B2*C[id[1]]*probi
  }
  
  if(length(id)>1){
    #B=Fast_FOI_function((id-1),centroids[,1,drop=FALSE],centroids[,2,drop=FALSE],
    #                    cells,F2_int,F2_B,F2i_int,F2i_B,t(repmat(I[id],cells,1)),t(repmat(C[id],cells,1)),B1,B2)
  
    B=Fast_FOI_function((id-1),centroids[,1,drop=FALSE],centroids[,2,drop=FALSE],
                        cells,F2_int,F2_B,F2i_int,F2i_B,I[id,,drop=FALSE],C[id,,drop=FALSE],B1,B2)
    
    }
  
  #Any B val < machine epsilon should just be zero
  #imprecise calculations below this value
  #B[B<.Machine$double.eps]<-0
  
  if(length(id)==1){
  Pse=1-exp(-W[,1]-W[,2]-t(B))
  } else{
    if(length(id)>1){
      Pse=1-exp(-W[,1]-W[,2]-t(colSums(B)))
    }
  }
  
  
  return(Pse)
  
}
