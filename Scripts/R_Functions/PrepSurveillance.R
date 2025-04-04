#######################
####### Purpose #######
#######################

# Processes sampling file that is user inputted indicating
# collection date, longitude, latitude, quantity and acreage of sample effort
  # Currently, only option is to have landscape prediction for county of interest
  # i.e. grid.opt="ras"

#######################
######## Function #####
#######################

# Inputs:
  # sampling- file with sampling scheme
  # inc - grid cell resolution (not sure if this will be true with new MakeGrid implementation)

#Outputs:
#1-sample.design- dataframe containing sampling info to be used in Surveilance.R
  #col1-date
  #col2-lat
  #col3-long
  #col4-fiscal week
  # Madison, finish this.


PrepSurveillance<-function(sampling){

  # Read sampling file
  sample.design <- read.csv("sampling_scheme.csv")

  # Aggregate by week and store weeks to be sampled in variable
  # Convert "collection_date" to Date format
  names(sample.design)[1] <- "dates"  # Rename the first column to 'dates'
  sample.design$dates <- as.Date(sample.design$dates, format = "%m/%d/%Y") # change to standard date format
  
  # Apply mutate on the entire dataframe, not just on a vector
  sample.design <- sample.design %>%
    mutate(# if collection date is in october or later, fiscal year starts in that year
      # can handle out of order dates from spreadsheet
      fiscal_year_start = as.Date(paste0(ifelse(month(dates) >= 10, year(dates), year(dates) - 1), "-10-01")), 
      # Which fiscal week are we in
      fiscal_week = as.integer(difftime(dates, fiscal_year_start, units = "weeks")) + 1,
      # Make week number stays between 1-52 (handle any overflow)
      fiscal_week = ifelse(fiscal_week < 1, fiscal_week + 52, fiscal_week),  # Handle weeks that are before the fiscal year start
      fiscal_week = ifelse(fiscal_week > 52, fiscal_week %% 52, fiscal_week)  # Ensure the week number stays between 1 and 52
    )
  
  # Decide which cells are sampling locations and add "1" to 8th column of grid
  names(sample.design)[2] <- "latitude"
  names(sample.design)[3] <- "longitude"
  
  sample_coords <- st_as_sf(sample.design, coords = c("longitude", "latitude"), crs=4269) # pull out x and y coords from sample.design, convert to sf object
  sample_sf_transformed <- st_transform(sample_coords, crs = 26917) # transform coordinates to meter based CRS
  sample_coords_transformed <- st_coordinates(sample_sf_transformed) # extract the coordinates only
  sample_coords_transformed <- sample_coords_transformed / scale_factor
  
  # Add column to sample.design so the cell # where sampling occurs can be updated
  sample.design$sampling_loc <- 0L
  
  # Loop over each sampling point and check proximity to grid centroids
  for (i in 1:nrow(sample_coords_transformed)) {
    # Extract the x and y coordinates of the current sample point
    sample_x <- sample_coords_transformed[i, 1]
    sample_y <- sample_coords_transformed[i, 2]
    
    # Get the acreage for the current sample point (assuming you have an 'acres' column in the dataframe)
    names(sample.design)[5] <- "acres"  # Rename the first column to 'dates'
    acres <- sample.design$acres[i]
    
    # Convert acres to square kilometers
    area_km2 <- acres * 0.00404686
    
    # Resolution of a grid cell in km2 calculation
    numerator = inc * 1000 * inc * 1000
    denominator = 1000000
    grid_cell_area = numerator/denominator
    
    # Calculate the number of grid cells to sample based on the area (rounding up to ensure entire area is covered)
    num_cells_to_sample <- ceiling(area_km2 / grid_cell_area)
    
    print(paste0("Number of cells to sample: ", num_cells_to_sample))
    
    # Calculate the Euclidean distance (dist between 2 points) from this sample point to each centroid in the grid
    # square root [(xf-xi)^2 + (yf-yi)^2]
    distances <- sqrt((grid[, 6] - sample_x)^2 + (grid[, 7] - sample_y)^2)
    
    # Check if the minimum distance is within the threshold
    # Using threshold of 10 meters
    if (min(distances) <= 10) {
      # Find the index of the closest centroid and set its match to 1
      closest_centroid_index <- which.min(distances)
      grid[, 8][closest_centroid_index] <- 1
      
      # Save the row number which corresponds to the grid cell number where sampling will occur
      sample.design$sampling_loc[i] <- closest_centroid_index
      
      # Now sample additional grid cells based on the area
      sampled_cells <- 1  # We've already sampled the closest one, so start counting from 1
      
      # We want to sample `num_cells_to_sample - 1` more cells
      sampled_cell_indices <- c(closest_centroid_index)  # Initialize with the first sampled cell index
      
      # We want to sample `num_cells_to_sample - 1` more cells
      while (sampled_cells < num_cells_to_sample) {
        # Find the next closest grid cells that haven't been sampled yet
        remaining_distances <- distances
        remaining_distances[grid[, 8] == 1] <- Inf  # Exclude already sampled cells
        
        # Find the index of the next closest centroid
        next_closest_index <- which.min(remaining_distances)
        
        # Mark this grid cell as sampled
        grid[, 8][next_closest_index] <- 1
        sampled_cells <- sampled_cells + 1  # Increment the counter of sampled cells
        
        # Append this grid cell index to the list of sampled cell indices
        sampled_cell_indices <- c(sampled_cell_indices, next_closest_index) 
      }  
      
      # Save the row number (cell number) in the sampling location
      sample.design$sampling_loc[i] <- list(sampled_cell_indices)  # Store as a list
      
    }
  }
  return(sample.design)
}