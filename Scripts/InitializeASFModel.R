#The purpose of this script is to set all the parameters
#Loads libraries, sources functions, all things not dependent on spec. parms
#Run this before running RunSimulationModel.R

######################
####Set directories
#####################
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/ASF_optimal_radius/"
#home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/ASF_optimal_radius"
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
library(stringr)
library(dplyr)


######################
####Source Functions
#####################

#To test compilation ability
#old_path <- Sys.getenv("PATH")
#Sys.setenv(PATH = paste(old_path, "/Library/Developer/CommandLineTools", sep = ":"))
#"/Applications/Xcode.app/Contents/Developer"
#Sys.which("make")

source(paste(getwd(), "/Scripts/ASFFunctionSourcer.R", sep = ''))


######################
####Import grid
######################

#grid is a matrix, with nrow=ncell of landscape grid
#col 1- sequence 1:ncell
#col 2-5: X1,Y1; X2,Y2
#col 6-7: centroid x and y
grid<-readMat(paste0(home,"Matlab_ASF_model/Grid_80x80_0pt4km.mat"))
grid<-grid$grid

#Need to look at code for making grid in matlab, make function to generate grid in R
#This will likely be something more complex in future, since grid will be incorporating other landscape/env elements
#think about settings that would be good to toggle

#function to control sounder size to grid size/resolution...
#maybe could save grids as a matrix within a list (in an RDS file), and include the resolution in the list
#then when bringing it in, use that part of the list to get the sounder size/etc.
#for now, just going to set sounder size manually because I only have the one grid anyways

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

#F1- pig contact probability if in same cell
#F2- pig contact probability with distance
#F2i- pig contact probability with carcass

#Movement parms
xFL=c(0.7515,0.3550)
xSC=c(0.5657,1.9082)

#source F2_FL and F2_SC
#source("GenerateFakeStateData.R",local=TRUE)
source(paste0(home,"/Scripts/Model_State_Data.R"))


