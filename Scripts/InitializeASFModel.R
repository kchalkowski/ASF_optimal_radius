#The purpose of this script is to set all the parameters 
#and initialize all state variables needed to run the ASF Model
#Run this before running SimulateOneRun.R

######################
####Set directories
#####################
#home<-"/Users/kchalkowski/Desktop/USDA_Pigs/Projects/ASF_simulation_model/Rcpp_ASF/ASF_optimal_radius"
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/ASF_optimal_radius"
setwd(home)

######################
####Set libraries
#####################
library(Rcpp)
library(profvis)
library(R.matlab)
library(pracma)
library(rdist)
library(tidyverse)
library(microbenchmark)
library(RcppArmadillo)
library(RcppParallel)


######################
####Source Functions
#####################

source(paste(getwd(), "/Scripts/ASFFunctionSourcer.R", sep = ''))

######################
####Define Variables
#####################

state=2 #state switch
loc=c("FL","SC") #state strings
density=3 #density per X
area=80^2 #total area of grid
N0=density*area #initial population size
K=N0*1.5 #carrying capacity for whole population
Rad=10 #culling radius
detectday = 20
thyme=72+detectday #using 'time' sometimes gives weird errors with different packages, hence the diff spelling
#thyme=72
#detectday=73
intensity=0.05
cullstyle="startOUT"
inc=0.4
Intensity = 0.05 #proportion of population targeted for removal per day based on capabilities
#stplot=1 #start plot at timestep i
#enplot=time+1 #end plot at timestep j

######################
####Import grid
######################

#grid is a matrix, with nrow=ncell of landscape grid
#col 1- sequence 1:ncell
#col 2-5: X1,Y1; X2,Y2
#col 6-7: centroid x and y
grid<-readMat(paste0(home,"/Matlab_ASF_Model/Grid_80x80_0pt4km.mat"))
grid<-grid$grid

#Need to look at code for making grid in matlab, make function to generate grid in R
#This will likely be something more complex in future, since grid will be incorporating other landscape/env elements
#think about settings that would be good to toggle

#function to control sounder size to grid size/resolution...
#maybe could save grids as a matrix within a list (in an RDS file), and include the resolution in the list
#then when bringing it in, use that part of the list to get the sounder size/etc. 
#for now, just going to set sounder size manually because I only have the one grid anyways

inc=0.4
ss=4
centroids=grid[,6:7]
cells=nrow(grid)
area<-grid[cells,4]*grid[cells,5]
midpoint=c(max(centroids[,1]/2),max(centroids[,2]/2))

#########################
####Import State contact rules
#########################

#Need get data and script for data 
#generate FL and SC contact rules with a glm in R
#for now, generate fake data for glm for indirect contact rules

xFL=c(0.7515,0.3550)
xSC=c(0.5657,1.9082)

F1_FL=0.7381
F1_SC=0.9375

#source F2_FL and F2_SC
#source("GenerateFakeStateData.R",local=TRUE)

#########################
####Load fixed parameters
#########################

Pir = 0.05 #; %proportion of individuals that recover instead of die
death = 7/(365*3) #assume pop growth rate of 1.5 so make death rate = birth rate*(1/1.5); 1/(365*3); % natural death rate for S and R
alphaC = 1.1 #; % scaling parameter on relationship of effort to capture success

#########################
####Choose State Rules
#########################
#define movement characteristics of the population (slow=FL, fast=SC)

if(state == 1){
	#specify which set if gamma fit parameters to use: xSC or xFL
	shift=xFL
	#specify which contact data to use (SC vs FL)
	F1=F1_FL 
	F2=F2_FL
	F2i=F2i_FL
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
	F1=F1_SC
	F2=F2_SC
	F2i=F2i_SC
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
Pbd = 7*mc_time #; %repmat(mean(c_time(1:364/7)),time,1).*1; % constant birth rate for S; rescale trend as needed to produce realistic pop dynamics

######################
####Initialize Population
#####################
pop<-InitializeSounders(N0,ss,cells,centroids,0,0,0)

######################
####RunModel
#####################

num_inf_0=1 #how many pigs to infect starting off

                                                                                                                                                                                                                                                                                                                                       
