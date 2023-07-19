# ASF optimal radius model 

This is an adaptation of the ASF meta-population model in Pepin et al. 2022, converted from Matlab to R, with changes made to optimize model speed and allow flexible incorporation of environmental and management-relevant parameters.

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

Objects: 
Nall- matrix of nrow=number of time steps, tracks number of (live) pigs in population at each time step
BB- matrix, nrow=number of time steps, tracks number of births at each time step
POSlive- list, length=number of time steps, positive (live) cases observed and removed from landscape
POSdead- list, length=number of time steps, positive carcasses observed and removed from landscape
NEGlive- list, length=number of time steps, negative (live) cases observed and removed from landscape
NEGdead- list, length=number of time steps, negative carcasses observed and removed from landscape
POSlive_locs- list, length=number of time steps, locations of all positive detected cases
POSdead_locs- list, length=number of time steps, locations of all positive detected carcasses
idZONE- list, length=number of time steps, each item in list contains a matrix of _______
Tculled- matrix, nrow=number of time steps, total number culled at each time step
ZONEkm2- matrix, nrow=number of time steps, total area of control zone at each time step in km2
Spread- matrix, nrow=number of time steps, ncol=3, col 1=number of infectious individuals, col 2=true area of infection (using convex hull around infected pigs), col 3=max distance between any two cases
Incidence- matrix, nrow=number of time steps, new (true) cases at each time step
I_locs- list, length=number of time steps, locations of all infected pigs
C_locs- list, length=number of time steps, locations of all infected carcasses
removalcells- list, length=number of time steps, cells where pig removals occurred during culling
Isums- matrix, nrow=number of time steps, number of rows in population matrix with infected pigs at each time step
Csums-matrix, nrow=number of time steps, number of rows in population matrix with infected carcasses at each time step
out-matrix with 3 cols, nrow=number of time steps, col 1=number of infected pigs/carcasses, col 2=area of infection, col 3=max distance between infected pigs/carcasses
ICtrue-vector containing total number of infections (I,C,E) at each time step
IConDD-number of ICE on detect day
ICatEnd-number of infection on last day
idNEW-vector containing newly detected infected pigs/carcasses, begins on day after detectday, reflecting a 1 day lag between initiating surveillance (initiated on detectday) and new detections 
idZONE-list with length=timestep, each item in list containing matrices with 3 cols, containing all grid cells in zone per row. col 1=grid cell, col2=each grid cell in zone paired with grid cell in col 1, col 3=distance between paired grid cells

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

**CullingOneRun (CullingOneRun.R)**  
Finds cells within set radius from cells where an infected pig/carcass was found in last time step, then combines with cells in zone from previous time steps. If there are pigs within the culling zone, remove some randomly based on set culling intensity. Output new population matrix, updated control zone cells (idZONE) and detection numbers and locations.

**areaOfinfection (areaOfinfection.R)**
Identifies where infected pigs/carcasses are, outputs convex hull area around infected pigs/carcasses, and max distance between them

### In progress
-add all current scripts to readme
-finish neatening/commenting scripts

