# Spatial meta-population model   

Based on ASF simulation model from Pepin et al. 2022, Optimizing response to an introduction of African Swine Fever in wild pigs, converted from Matlab to R/C++.    

Readme/scripts updated periodically.    

## Steps to run simulation
1. clone repo to your machine
2. store contact data files in Input folder (email kim.m.pepin@usda.gov or kayleigh.chalkowski@usda.gov for these files)
3. set home directory in RunASFSimReplicates.R to your ASF_optimal_radius path
4. if needed, change parameters in SetParameters.R

## Recent Updates
**September 18, 2024** Added functionality for additional output/aggregation options for sounderlocs output, including some Rcpp-optimized spatial summary options (currently demonstrated in Sensitivity_Analysis/Scripts/doSensitivityAnalysis.R). Added comments to new spatial functions and finalized sensitivity analysis pipeline for comparing simulation output to observed, de-identified ASF outbreak data.    
**August 28, 2024** Optimized FOI function-- created FOI_Fast_Matrix.cpp to speed up FOI calcs via Rcpp armadillo, loop unrolling, and added distance cutoff for calculations.

## Pending/Upcoming Updates
* Improve generalizeability of input options for InitializeSounders.R.   
* Add RSF-preference option to InitializeSounders.R.   
* Improve generalizeability of Movement_Parallel_Functionsmall.cpp, make abundance-driven movement an option rather than required
* Add RSF-preference option to Movement_Parallel_Functionsmall.cpp
* Add functionality to MakeGrid for creating different landscape structures with landscapemetrics
* Add external input option to FOI_Fast_Matrix.cpp for FOI_cutoff to allow more flexiblity in determining calculation cutoffs-- currently just hardcoded to 5km in code on line 67





