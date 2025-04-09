# Spatial meta-population model   

Based on ASF simulation model from Pepin et al. 2022, Optimizing response to an introduction of African Swine Fever in wild pigs, converted from Matlab to R/C++.    

Readme/scripts updated periodically.    

## Steps to run simulation
1. Clone repo to your machine
2. Store contact data files in Input folder (email kim.m.pepin@usda.gov and kayleigh.chalkowski@usda.gov for these files)
3. Set home directory in RunASFSimReplicates.R to your ASF_optimal_radius path
4. If needed, change parameters in SetParameters.R

## Recent Updates
**October 9, 2024:**  
Added capability in Make_Grid to input external raster with lc values.

**October 9, 2024:**  
Incorporated new movement function, Movement_Fast_RSFavail.cpp (parm setting mv_pref=3 to use), which allows for flexibility in RSF preferences according to availability. Requires matrix of RSF probabilities with different combinations, with corresponding dummy coded RSF availability matrix.

**September 18, 2024:**   
Added options to initialize (neutral landscape model) heterogeneous landscape, initialize population according to RSF probabilities assigned to landscape, and allow for RSF-driven movement preference. Made movement function more general to allow for different movement preference types (distance-only, abundance-avoidant, and rsf land class preference).    

**September 18, 2024:**    
Added functionality for additional output/aggregation options for sounderlocs output, including some Rcpp-optimized spatial summary options (currently demonstrated in Sensitivity_Analysis/Scripts/doSensitivityAnalysis.R). Added comments to new spatial functions and finalized sensitivity analysis pipeline for comparing simulation output to observed, de-identified ASF outbreak data.    
    
**August 28, 2024:**    
Optimized FOI function-- created FOI_Fast_Matrix.cpp to speed up FOI calcs via Rcpp armadillo, loop unrolling, and added distance cutoff for calculations.

## Pending/Upcoming Updates
* Transitioning model to targets pipeline
* Add external input option to FOI_Fast_Matrix.cpp for FOI_cutoff to allow more flexiblity in determining calculation cutoffs-- currently just hardcoded to 5km in code on line 67
* Develop and test options for sounder initialization on grid when RSF varies with availability



