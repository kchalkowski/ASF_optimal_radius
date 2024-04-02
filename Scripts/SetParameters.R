#########################
####Load fixed parameters
#########################

Pir = 0.05 #; %proportion of individuals that recover instead of die
death = 7/(365*3) #assume pop growth rate of 1.5 so make death rate = birth rate*(1/1.5); 1/(365*3); % natural death rate for S and R
alphaC = 1.1 #; % scaling parameter on relationship of effort to capture success

######################
####Define Variables
#####################

state=1 #state switch
loc=c("FL","SC") #state strings
#density=5 #density per X
area=80^2 #total area of grid
N0=density*area #initial population size
K=N0*1.5 #carrying capacity for whole population
Rad=0 #culling radius
#detectday = 72
detectday=73
#thyme=52+detectday #using 'time' sometimes gives weird errors with different packages, hence the diff spelling
thyme=72
intensity=0.05
cullstyle="startOUT"
inc=0.4
Intensity = 0.05 #proportion of population targeted for removal per day based on capabilities
#ss=4
num_inf_0=1 #how many pigs to infect starting off
#stplot=1 #start plot at timestep i
#enplot=time+1 #end plot at timestep j


if(state == 1){
  #specify which set if gamma fit parameters to use: xSC or xFL
  shift=xFL
  #specify which contact data to use (SC vs FL)
  F1=F1.list$FL
  F2=F2.list$FL
  F2i=F2i.list$FL
  if(density == 1.5){
    B1=0.9
  } else if(density == 3){
    B1=0.4
  } else if(density == 5){
    B1=0.2
  }
} else if(state == 2){
  #specify which set if gamma fit parameters to use: xSC or xFL
  shift=xSC
  #specify which contact data to use (SC vs FL)
  F1=F1.list$SC
  F2=F2.list$SC
  F2i=F2i.list$SC
  if(density == 1.5){
    B1=0.009
  } else if(density == 3){
    B1=0.004
  } else if(density == 5){
    B1=0.002
  }
}

B2 = B1*0.5;

#########################
####Create Seasonally Varying Parameters
#########################
#*For now they are constant

mc_time=0.0027
Pcr = 7/30 #fix at 30 days - average for the year; %zeta.*1; % transition probability from carcass to removal from landscape
#Pbd = 7*mc_time #; %repmat(mean(c_time(1:364/7)),time,1).*1; % constant birth rate for S; rescale trend as needed to produce realistic pop dynamics
Pbd=0.0192 #From opt rad ms

#Set F2/F2i manually here
#need as explicit parameters to run FOI in cpp
#will later be incorporating this once finish U.S. contact prediction mapping
#infectious pig contact prob based on distance
F2_int=F2$coef[[1]]
F2_B=F2$coef[[2]]
#infected carcass contact prob based on distance
F2i_int=F2i$coef[[1]]
F2i_B=F2i$coef[[2]]


