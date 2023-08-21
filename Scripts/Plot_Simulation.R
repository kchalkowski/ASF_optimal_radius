#function to plot simulation output
#Input: I_locs, C_locs, POSlive_locs, POSdead_locs
#Output: plot of simulation progress at each timestep

#PlotSimulation(I_locs, C_locs, POSlive_locs, POSdead_locs)

PlotSimulation<-function(I_locs, C_locs, POSlive_locs, POSdead_locs){

I_locs.t=transform(reshape2::melt(list.all$I_locs), id = sequence(lengths(list.all$I_locs)))
C_locs.t=transform(reshape2::melt(list.all$C_locs), id = sequence(lengths(list.all$C_locs)))
POSlive_locs.t=transform(reshape2::melt(list.all$POSlive_locs), id = sequence(lengths(list.all$POSlive_locs)))
POSdead_locs.t=transform(reshape2::melt(list.all$POSdead_locs), id = sequence(lengths(list.all$POSdead_locs)))

colnames(I_locs.t)[1]<-"I_locs"
colnames(C_locs.t)[1]<-"C_locs"
colnames(POSlive_locs.t)[1]<-"POSlive_locs"
colnames(POSdead_locs.t)[1]<-"POSdead_locs"

I_locs.t<-I_locs.t[,1:2]
C_locs.t<-C_locs.t[,1:2]
POSlive_locs.t<-POSlive_locs.t[,1:2]
POSdead_locs.t<-POSdead_locs.t[,1:2]

IC_locs=left_join(I_locs.t,C_locs.t,key="L1")
IC_locs2=left_join(IC_locs,POSlive_locs.t,key="L1")
IC_locs3=left_join(IC_locs2,POSdead_locs.t,key="L1")

IC_locs3[IC_locs3==0]<-NA
#which(IC_locs3[,c(1,3,4,5)]
IC_locs3=IC_locs3[rowSums(is.na(IC_locs3)) != ncol(IC_locs3),]

IC_locs3=as.data.frame(pivot_longer(IC_locs3,cols=c(1,3,4,5),names_to="type",values_to="cellnum"))
IC_locs3=IC_locs3[complete.cases(IC_locs3),]

IC_locs3$x=centroids[IC_locs3$cellnum,1]
IC_locs3$y=centroids[IC_locs3$cellnum,2]

#2-make radius data frame
	#current timestep
	#POSlive cellnum from last time step, is center of each radius
	#x and y coords

IC_locs4=IC_locs3[!duplicated(IC_locs3),]

#2- make rad data frame, join to IC locs
#Need to loop through IC_locs4
#starting day after detection day
#needs be data frame with timestep, and rad
#for each timestep, all POSlive/dead pts become rad centers

POS=IC_locs4[which(IC_locs4$type=="POSlive_locs"|IC_locs4$type=="POSdead_locs"),]

radd=data.frame("L1"=integer(),"x.rad"=numeric(),"y.rad"=numeric())
for(i in 1:nrow(POS)){
L1=POS[i,1]+1
x.rad=POS[i,4]
y.rad=POS[i,5]

rad1=data.frame("L1"=L1,"x.rad"=x.rad,"y.rad"=y.rad)
radd=rbind(radd,rad1)

if(L1!=max(IC_locs4$L1)+1){
for(j in 1:(max(IC_locs4$L1)-L1)){
radj=data.frame("L1"=L1+j,"x.rad"=x.rad,"y.rad"=y.rad)
radd=rbind(radd,radj)
}
}	
}

IC_locsj=full_join(IC_locs4,radd,by="L1")
test=IC_locsj[,2:5]
IC_locs.plot=IC_locsj[-which(rowSums(is.na(test))==4),]

#3- Plot all
ggplot(IC_locs.plot, aes(x, y,color=as.factor(type)))+
	geom_circle(aes(x0 =x.rad, y0 =y.rad, r = 10),fill="#D2EBAC",linetype=0, inherit.aes = FALSE)+
geom_point(size=0.05)+#color=as.factor(ICdf$I_C))+
		facet_wrap(IC_locs.plot$L1)+
	theme(#axis.text.x=element_blank(),
          #axis.text.y=element_blank(),
					#axis.ticks=element_blank(),
          #axis.title.x=element_blank(),
					panel.grid.minor = element_blank(),
					panel.grid.major = element_blank(),
					panel.background = element_blank(),
					panel.border = element_rect(colour = "black", fill=NA, size=0.5),
					strip.background = element_rect(color="black", fill="black", size=1.5, linetype="solid"),
				  strip.text.x = element_text(size = 8, color = "white",face="bold"),
      		strip.text.y = element_text( size = 8, color = "white",face="bold"))+
          #axis.title.y=element_blank())+
					xlim(0,80)+ylim(0,80)	
	}