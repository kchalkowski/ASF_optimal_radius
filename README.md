# ASF optimal radius model 

The purpose of this Readme is to describe the pipeline for our African Swine Fever simulation model. This is an adaptation of the ASF meta-population model in Pepin et al. 2022, converted from Matlab to R, with changes made to optimize model speed and allow flexible incorporation of environmental and management-relevant parameters.

### Table of Contents:    
1. Scripts summaries
2. Function summaries
3. Data descriptions    
4. Pipeline 
5. Troubleshooting needed

### Scripts Overview

**ASFFunctionSourcer.R**    
Sources all the functions used in the model. No inputs/outputs. 

**InitializeASFModel.R**    
Sets parameters, runs ASFFunctionSourcer.R to source functions,  initializes state variables used in the simulation, loads data and grid, sets number to infect at first time step. This is run before starting the simulation.

**SimulateOneRun.R**    
Runs the simulation. Outputs sum of all exposures, total number culled at each time step, last day there is an infectious individual, max spread of infection, number of I,C,E (infected, infected carcass, exposed) individuals on detection day, number of I,C,E on last day, sum of all exposures starting day after detection day, total number of detections, locations of infected at each time step, locations of carcasses at each time step, Isums, Csums.

### Functions    

**InitializeSounders (InitializeSounders.R)**     
Purpose is to initialize the starting population matrix in InitializeASFModel or add new pigs/sounders to existing population

**FastMovement (FastMovement.R)**        
Assigns distances using a gamma distribution, parameterized by collar data, then runs Rcpp function parallelMovementRcpp_portion, which conducts the movement process and outputs the new locations to pop[,3] (present location cell numbers). Outputs an updated population matrix with the new locations.

**parallelMovementRcpp_portion (Movement_Parallel_Functionsmall.cpp)**        
This is the Rcpp function that conducts the movement process. As inputs, it takes the population matrix, two columns subsetted from that population matrix (cols 1 and 3) and the centroids matrix. 
Function loops through each row of pop (each sounder), and if assigned movement distance is greater than 0, and if there is at least one pig in sounder, then loops through each cell in centroids to get difference between all the distances between that sounder and all other cells in the grid, and the assigned distance. Then, select set of cells that are closest to assigned distance. From this set, select cell to move to with minimum abundance. If multiple cells in set have same abundance, select a cell from those at random.

**StateChanges (StateChanges.R)**    
This function conducts all state changes including births, natural deaths, exposure (via the force of infection function), recovery, disease mortalities, and carcass decay. It outputs the pop matrix, and output vectors Incidence (tracks exposures over time) and BB (births over time). 

**FirstDetect (FirstDetect.R)**    
This function is only run in SimulateOneRun if the timestep is equal to detectday. The function detects an infected pig or carcass at random, records the infection in either POSlive or POSdead for that time step, removes the detected live pig or carcass, and returns the updated population matrix and POSlive/POSdead vectors

### In progress
1. more optimize on movement func
2. more optimize on FOI
3. add formal sanity checks to funcs

