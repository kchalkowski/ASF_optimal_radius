
#########################
######## Purpose ######## 
#########################

#Initializes grid to run spatially-explicit meta-population disease spread model
#Includes options for either homogeneous, uniform grid, or grid with land class preference probabilities
  #Currently only option is random landscape, alternative clustering options will be available in future updates.

##########################
######## Function ######## 
##########################

#Inputs: 
  #len of each side (in cells);
  #resolution of grid (in km2); 
  #grid.opt="homogeneous","random", or "ras"
  #ras, raster object, optional input

#Outputs: 
#1-grid matrix: 
  #col1- list of cells 1:ncell
  #col2- left-most x coord of cell
  #col3- upper-most y coord of cell
  #col4- right-most x coord of cell
  #col5- lower-most y coord of cell
  #col6- center x coordinate of cell
  #col7- center y coordinate of cell
  #if grid.opt="heterogeneous":
  #col8- preference probability
#2- cells, integer of number of cells in grid
#3- centroids, two col x/y of just centroids of grid of each cell
#Make_Grid(len,inc,grid.opt)


#input options
#required: object
  #raster, or c(len, inc)
#optional: grid.opt, 
#Make_Grid(object)

Make_Grid<-function(object,grid.opt="homogeneous",sample=0,sample.design=NULL){
  require(raster)
  require(NLMR)
  
  if(grid.opt=="homogenous"){
    grid.opt="homogeneous"
  }

  if(class(object)=="numeric"){
    len=object[1]
    inc=object[2]
  }
  
  if(class(object)=="RasterLayer"){
    ras=object
    len=dim(ras)[1]
    inc=res(ras)[1]/1000
    if(dim(ras)[1]!=dim(ras)[2]){stop("raster needs to be square")}

  
      if(dim(ras)[1]!=len){
        stop("dimensions of raster do not match input len")
      }
  
  }
  
  #get number of cells in grid
  cells=len^2
  
  #if grid homogeneous-- will later enter ability to alter LULC
  if("homogeneous"%in%grid.opt & sample != 1){
    #initialize empty grid matrix
    grid=matrix(nrow=round(cells),ncol=7)
  } else {
    #initialize empty grid matrix
    grid=matrix(nrow=round(cells),ncol=8)
  }
  

  #first column is just cell indices
  grid[,1]=1:cells
  
  #Top left X coordinate of each cell
  grid[,2]=rep(seq(0,((inc*len)-inc),inc),times=len)

  #Top left Y coordinate of each cell
  grid[,3]=rep(seq(0,((inc*len)-inc),inc),each=len)
  
  #Top right X coordinate of each cell
  grid[,4]=rep(seq(inc,(inc*len),inc),times=len)

  #Top right Y coordinate of each cell
  grid[,5]=rep(seq(inc,(inc*len),inc),each=len)
  
  #Center X coordinate of each cell
  grid[,6]=rep(seq(((0+inc)/2),(((inc*len)-inc)+(inc*len))/2,inc),times=len)

  #Center Y coordinate of each cell
  grid[,7]=rep(seq(((0+inc)/2),(((inc*len)-inc)+(inc*len))/2,inc),each=len)
  
  #get centroids-only object
  centroids=grid[,c(6,7)]
  
  if(!("homogeneous"%in%grid.opt & sample != 1)){
    
    #simulates a spatially random neutral landscape model with values drawn from a uniform distribution
    #values rescaled to range from 0-1
    if("random"%in%grid.opt){
    r=NLMR::nlm_random(len,len,inc,rescale=TRUE)
    grid[,8]=round(values(r),2)

    centroids=cbind(centroids,grid[,8])
    grid.list=list("cells"=cells,"grid"=grid,"centroids"=centroids,"r"=r)
    
    }
    
    if("ras"%in%grid.opt){
      #need to get values from ras
      grid[,8]=round(values(ras),2)
      #assign to centroids
      centroids=cbind(centroids,grid[,8])
      grid.list=list("cells"=cells,"grid"=grid,"centroids"=centroids)
    }
    
  } else{
    grid.list=list("cells"=cells,"grid"=grid,"centroids"=centroids)
  }
  
  if ("SpatRaster"%in%grid.opt & sample == 1) { # this loop will take in spatraster with land class pref
    # NOTE: This if statement has not been tested. Madison is pausing and reverting back to original
    # surveillance implementation.
    
    # Add a line to resize based on county (remove this in final version, just for testing purposes)
    r <- rast("Input/lands/pland_01.tif") # Not sure where the actual tif gets loaded in this pipeline
    new_extent <- ext(-78.90200, -78.17650, 34.36762, 34.85618)
    ext(r) <- new_extent
    ncell(r)
    plot(r)
    
    # Create grid matrix with info from spatraster
    grid <- matrix(NA, nrow = ncell(r), ncol=8) # total number of cells in grid
    grid[, 1] <- rep(1:ncell(r), each = 1)  # total number of cells in grid
    
    centroids <- xyFromCell(r, 1:ncell(r))
    centroids_df <-as.data.frame(coordinates)
    grid[, 6] <- centroids_df[, 1]  # Centroid x value provided by terra
    grid[, 7] <- centroids_df[, 2]  # Centroid y value provided by terra
    
    # Get the resolution (cell size) of the raster
    cell_size_x <- res(r)[1]  # Resolution in the x-direction
    cell_size_y <- res(r)[2]  # Resolution in the y-direction
    grid[, 2] <- rep(centroids_df[, 1] - (cell_size_x / 2)) # get rid of last one
    grid[, 3] <- rep(centroids_df[, 2] + (cell_size_y / 2)) # get rid of last one
    grid[, 4] <- rep(centroids_df[, 1] + (cell_size_x / 2)) # get rid of first one
    grid[, 5] <- rep(centroids_df[, 2] + (cell_size_y / 2)) # same as grid[, 3]

    # Create centroids matrix
    centroids <- cbind(grid[, 6], grid[, 7])  # Combine center X and Y as a matrix
    
    ######## Build final output needed for code to continue running
    # Kept it the same format as the homogeneous grid like Kayleigh's from above
    # cell index, top left x, top left y, top right x, top right y, centroid x, centroid y, sampling location

    cells <- nrow(grid)

    # Output list
    grid.list <- list(
      "cells" = cells, 
      "grid" = grid,
      "centroids" = centroids
    )
    
    # Transform spat rast from 4326 to 26917 CRS
    r_transformed <- project(r, crs = 26917)  # UTM Zone 33N (example)
    res_x <- res(raster_utm)[1]  # Resolution in x direction (in meters)
    res_y <- res(raster_utm)[2]  # Resolution in y direction (in meters)
    grid_cell_area <- res_x * res_y

    # Add column to sample.design so the cell # where sampling occurs can be updated
    sample.design$sampling_loc <- 0L
    
    # Loop over each sampling point and check proximity to grid centroids
    for (i in 1:nrow(sample_coords_transformed)) {
      # Extract the x and y coordinates of the current sample point
      sample_x <- sample_coords_transformed[i, 1]
      sample_y <- sample_coords_transformed[i, 2]
      
      # Get the acreage for the current sample point (assuming you have an 'acres' column in the dataframe)
      names(sample.design)[5] <- "acres"  # Rename the fifth column to 'acres'
      acres <- sample.design$acres[i]
      
      # Convert acres to square kilometers
      area_km2 <- acres * 0.00404686
      
      # Calculate the number of grid cells to sample based on the area (rounding up to ensure entire area is covered)
      num_cells_to_sample <- ceiling(area_km2 / grid_cell_area)

      
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
          
          # Mark this grid cell as to be sampled
          grid[, 8][next_closest_index] <- 1
          sampled_cells <- sampled_cells + 1  # Increment the counter of sampled cells
          
          # Append this grid cell index to the list of sampled cell indices
          sampled_cell_indices <- c(sampled_cell_indices, next_closest_index) 
        }  
        
        # Save the row number (cell number) in the sampling location
        sample.design$sampling_loc[i] <- list(sampled_cell_indices)  # Store as a list
        
      }
      # Plot the raster and sampling locations
      plot(r_transformed)
      points(sample_coords_transformed$geometry, col = "red", pch = 19)
    }
    
    #  return(list(grid.list = grid.list, sample.design = sample.design))
    grid.list=list(grid.list = grid.list, sample.design = sample.design)
    
    }
  
  if ("homogeneous"%in%grid.opt & sample == 1) { # this loop will address homogeneous landscape and builds its own grid (grid.type = "County")
    
    # county CRS needs to be 26917 so that it can convert to the scale of the grid (meters)
    county_reprojected <- st_transform(county_shapefile, crs = 26917)
    
    # Get the bounding box dimensions
    county_bbox <- st_bbox(county_reprojected)
    
    # want grid to expand past county boundaries
    expansion_percentage <- 0.002  # 5% expansion
    
    # Get boundary around county and expand the bounding box slightly beyond the county boundaries based on expansion %
    xmin_expanded <- county_bbox["xmin"] - (county_bbox["xmin"] * expansion_percentage)
    xmax_expanded <- county_bbox["xmax"] + (county_bbox["xmax"] * expansion_percentage)
    ymin_expanded <- county_bbox["ymin"] - (county_bbox["ymin"] * expansion_percentage)
    ymax_expanded <- county_bbox["ymax"] + (county_bbox["ymax"] * expansion_percentage)
    
    # Calculate the width and height of the expanded bounding box and make sure grid is square
    width <- xmax_expanded - xmin_expanded
    height <- ymax_expanded - ymin_expanded
    grid_size <- max(width, height)
    
    # Calculate the number of cells in x and y direction based on the resolution (in meters)
    # and round to nearest integer [ceiling function]
    len <- ceiling(grid_size / (inc*1000))  # Number of cells along x or y axis
    
    # Create grid coordinates for x and y within the expanded bounding box
    # seq creates a sequence of numbers based on the x and y inputs you give it which are lists of possible points
    grid_x <- seq(xmin_expanded, xmax_expanded, length.out = len + 1)
    grid_y <- seq(ymin_expanded, ymax_expanded, length.out = len + 1)
    
    scale_factor <- 1000
    
    # Create grid matrix
    grid <- matrix(NA, nrow = len * len, ncol=8) # total number of cells in grid
    grid[, 1] <- rep(1:(len * len), each = 1)  # total number of cells in grid
    grid[, 2] <- rep(grid_x[-(len + 1)], times = len)  # create top left x coordinate but get rid of very last one
    grid[, 2] <- (grid[, 2]) / scale_factor  # Convert from meters to kilometers
    grid[, 3] <- rep(grid_y[-(len + 1)], each = len)  # create top left y coordinate but get rid of very last one
    grid[, 3] <- (grid[, 3]) / scale_factor  # Convert from meters to kilometers
    grid[, 4] <- rep(grid_x[-1], times = len)  # create top right x coordinate but get rid of first one
    grid[, 4] <- (grid[, 4]) / scale_factor
    grid[, 5] <- grid[, 3]  # Top-right Y is the same as Top-left Y
    grid[, 6] <- (grid[, 2] + grid[, 4]) / 2  # Center X (average of top-left and top-right X)
    grid[, 7] <- (grid[, 3] + grid[, 5]) / 2  # Center Y (average of top-left and top-right Y)
    
    # Create centroids matrix
    centroids <- cbind(grid[, 6], grid[, 7])  # Combine center X and Y as a matrix
    
    ######## Build final output needed for code to continue running
    # Kept it the same format as the homogeneous grid like Kayleigh's from above
    # cell index, top left x, top left y, top right x, top right y, centroid x, centroid y, sampling location
    
    cells <- nrow(grid)
    
    # Output list
    grid.list <- list(
      "cells" = cells, 
      "grid" = grid,
      "centroids" = centroids
    )
    
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
    
    #  return(list(grid.list = grid.list, sample.design = sample.design))
    grid.list=list(grid.list = grid.list, sample.design = sample.design)
    
  }
  
  return(grid.list)
  
}


