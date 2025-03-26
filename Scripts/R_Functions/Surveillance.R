Surveillance <- function(pop, i, sample.design, grid.list, inc, POSlive, POSdead, POSlive_locs, POSdead_locs, pigs_sampled_timestep) {
  
  # Check if the current fiscal week exists in sample.design
  rows_to_sample <- sample.design[sample.design$fiscal_week == i, ]
  
  if (nrow(rows_to_sample) == 0) {
    POSlive[[i]] <- 0
    POSdead[[i]] <- 0
    POSlive_locs[[i]] <- 0
    POSdead_locs[[i]] <- 0
    pigs_sampled_timestep[[i]] <- 0
    #print("No surveillance this week. Carry on, pigs.")
    #print(paste0("# of pigs sampled this timestep: ", pigs_sampled_timestep[[i]]))
    return(list(pop, POSlive, POSdead, POSlive_locs, POSdead_locs, pigs_sampled_timestep))
  }
  
  #print("Conducting surveillance...")
  
  # Initialize
  pigs_sampled <- 0  # Track pigs sampled so far
  infected_pigs_found <- FALSE  # Flag to track if any infected pigs are found
  
  # Loop through the rows in sample.design for the fiscal week to identify sampling locations
  for (row_index in 1:nrow(rows_to_sample)) {
    sampling_cells <- rows_to_sample$sampling_loc[[row_index]]  # Cells to sample for this location
    current_quantity <- rows_to_sample$Quantity[row_index]  # Quantity to sample for this location
    cells_with_pigs <- character(0)  # Reset cells for this location
    
    # Loop through the cells specified for this location
    for (sampling_cell in sampling_cells) {
      matching_rows <- which(pop[, 3] == sampling_cell)  # Find matching rows in pop
      if (length(matching_rows) > 0) {
        cells_with_pigs <- unique(c(cells_with_pigs, sampling_cell))  # Add unique cells with pigs to the list
      }
    }
    
    # Check if there are enough pigs available for sampling at this location
    pigs_available_to_sample <- sum(pop[, 3] %in% cells_with_pigs)
    
    if (pigs_available_to_sample == 0) {
      # print(paste("No pigs available to sample at location. Skipping location."))
      next  # Skip this iteration and move to the next location
    } else {
      # print(paste("Available pigs to sample at location:", pigs_available_to_sample))
    }
    
    # Sample pigs until the required quantity is met
    while (pigs_sampled < current_quantity && length(cells_with_pigs) > 0) {
      selected_cell <- sample(cells_with_pigs, 1)  # Randomly select a cell
      matching_rows <- which(pop[, 3] == selected_cell)  # Find rows for this cell
      
      pigs_found_in_cell <- FALSE  # Flag to check if pigs are found in the selected cell
      
      # Iterate over the rows in the selected cell
      for (row in matching_rows) {
        # Infectious pigs (columns 9 and 10)
        if (any(pop[row, c(9, 10)] > 0)) {
          pigs_sampled <- pigs_sampled + 1
          pop <- pop[-row, ]  # Remove infected pig
          POSlive[[i]] <- 1
          POSlive_locs[[i]] <- selected_cell
          pigs_found_in_cell <- TRUE
          infected_pigs_found <- TRUE  # Mark that an infected pig was found
        }
        # Non-infectious pigs (columns 8 and 11)
        else if (any(pop[row, c(8, 11)] > 0)) {
          pigs_sampled <- pigs_sampled + 1
          pop <- pop[-row, ]  # Remove non-infected pig
          pigs_found_in_cell <- TRUE
        }
        
        # If we reach the required quantity, break out of the loop
        if (pigs_sampled >= current_quantity) break
      }
      
      # If no pigs were found in the selected cell, remove it from the list
      if (!pigs_found_in_cell) {
        cells_with_pigs <- setdiff(cells_with_pigs, selected_cell)
        # print(paste("No pigs in cell", selected_cell, "- removing from pool."))
      }

      # If we haven't reached the required quantity, continue sampling
      # if (pigs_sampled < current_quantity) {
      #   print(paste("Still need more pigs. Pigs sampled so far:", pigs_sampled))
      # }
    }
    
    # Final message if the sampling doesn't reach the required quantity for this location
    # if (pigs_sampled < current_quantity) {
    #   print(paste("Not enough pigs found to meet sampling requirement at location", row_index))
    # } else {
    #   print(paste("Surveillance complete for location", row_index))
    # }
  }
  
  # Add the final message to indicate if any infected pigs were found
  # if (infected_pigs_found) {
  #   print("At least one infected pig was found during surveillance (columns 9 or 10 > 0).")
  # } else {
  #   print("No infected pigs were found during surveillance.")
  # }
  
  # Save the total number of pigs sampled for this fiscal week
  pigs_sampled_timestep[[i]] <- pigs_sampled
  #print(paste0("# of pigs sampled this timestep: ", pigs_sampled_timestep[[i]]))
  # print(paste0("Length of variable that POSlive is: ", length(POSlive)))
  # print(paste0("Length of variable that pigs_sampled_timestep is: ", length(pigs_sampled_timestep)))
  
  return(list(pop, POSlive, POSdead, POSlive_locs, POSdead_locs, pigs_sampled_timestep))
}