
#Loads libraries, sources functions, loads spatial objects 
#all things not dependent on specific parms
#Run this before running RunSimulationModel.R

######################
####Load libraries
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
library(stringr)
library(dplyr, warn.conflicts = FALSE)

# Suppress summarise info from dplyr
options(dplyr.summarise.inform = FALSE)

######################
####Source Functions
#####################

source(paste(home, "/Scripts/Setup/ASFFunctionSourcer.R", sep = ''))

######################
####Import grid
######################

#grid is a matrix, with nrow=ncell of landscape grid
#col 1- sequence 1:ncell
#col 2-5: X1,Y1; X2,Y2
#col 6-7: centroid x and y
if(grid.type=="ML"){
grid<-readMat(paste0(home,"/Input/Grid_80x80_0pt4km.mat"))
grid<-grid$grid

#scale land class values by RSFs
#if(grid.opts!="homogenous"){
#  prefs=grid[,8]
#  prefs[prefs==1]=RSF_mat[RSF_mat[,1]==1,2] #convert to RSF_prefs
#  prefs[prefs==0]=RSF_mat[RSF_mat[,1]==0,2] #convert to RSF_prefs
#  grid[,8]=prefs
#}

#set up objects
centroids=grid[,6:7]
cells=nrow(grid)
area<-cells*inc
midpoint=c(max(centroids[,1]/2),max(centroids[,2]/2))
} 

if(grid.type=="DIY"){
  grid.list=Make_Grid(len,inc,"homogenous")
  cells=grid.list$cells
  grid=grid.list$grid
  centroids=grid.list$centroids
  area=cells*inc
  midpoint=c(max(centroids[,1]/2),max(centroids[,2]/2))

#grid-dependent parms
area=len^2*inc #total area of grid
N0=density*area #initial population size
K=N0*1.5 #carrying capacity for whole population
}
#########################
####Import State contact rules
#########################

#Need get data and script for data
#generate FL and SC contact rules with a glm in R
#for now, generate fake data for glm for indirect contact rules

#F1- pig contact probability if in same cell
#F2- pig contact probability with distance
#F2i- pig contact probability with carcass

#Movement parms
xFL=c(0.7515,0.3550)
xSC=c(0.5657,1.9082)

#source F2_FL and F2_SC
source(paste0(home,"/Scripts/Setup/Model_State_Data.R"))

if(state == 1){
  #specify which set if gamma fit parameters to use: xSC or xFL
  shift=xFL
  #specify which contact data to use (SC vs FL)
  F1=F1.list$FL
  F2=F2.list$FL
  F2i=F2i.list$FL
  if(density == 1.5){
    B1=0.9
    ss=2
  } else if(density == 3){
    B1=0.4
    ss=4
  } else if(density == 5){
    B1=0.2
    ss=6
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
    ss=2
  } else if(density == 3){
    B1=0.004
    ss=4
  } else if(density == 5){
    B1=0.002
    ss=6
  }
}

B2 = B1*0.5;

#Set F2/F2i here
#need as explicit parameters to run FOI in cpp
#will later be incorporating this once finish U.S. contact prediction mapping
#infectious pig contact prob based on distance
F2_int=F2$coef[[1]]
F2_B=F2$coef[[2]]
#infected carcass contact prob based on distance
F2i_int=F2i$coef[[1]]
F2i_B=F2i$coef[[2]]

##### Determine cutoff for simulation model
#runs gamma for approximate number of distances that will be generated
#This will be used to determine cutoff of distances for which FOI will be calculated
#helps to cut down on FOI calculation time, and reduces inaccurate calculations with numbers < machine.epsilon
FOI_cutoff=round(max(rgamma(2*nrow(pop)*thyme,shape=shift[1],scale=shift[2])),2)

