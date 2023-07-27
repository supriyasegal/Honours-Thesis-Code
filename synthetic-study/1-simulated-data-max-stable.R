# This code:
# a) Creates synthetic data on a grid and simulates data from a max-stable model
# b) Fits a max-stable model to the simulated data in (a) and evaluates it for expected extreme rainfall estimation
# c) Fits a GEV to the synthetic data in (a) and evaluates the max-stable model in (b) for return levels estimation

library(SpatialExtremes)
library(ismev)

n.simulations <- 100
average_rmse <- rep(0, n.simulations)
return_level_10_rmse <- rep(0, n.simulations)
return_level_100_rmse <- rep(0, n.simulations)

for (simulation_no in 1:n.simulations) { 
  print(simulation_no)
  
  # CREATE SYNTHETIC DATA
  
  ## 1 - Simulate locations on a grid, and simulate data from a max-stable
  n.site = 35 # number of "seen" locations
  n.unobserved = 15 # number of "unseen" locations
  n.years = 200 # length of data series to simulate
  n.gridsize = 100 
  x <- runif(n.site + n.unobserved, min = 0, max = n.gridsize)
  y <- runif(n.site + n.unobserved, min = 0, max = n.gridsize)
  coord <- cbind(x, y)
  colnames(coord) <- c("lon", "lat")
  coord_observed <- coord[1:n.site, ]
  coord_unobserved <- coord[(n.site + 1):(n.site + n.unobserved), ]
  data <- rmaxstab(n.years, coord, cov.mod = "brown", range = 3, smooth = 0.5)
  
  ## 2 - Transformation to GEV margins
  param.loc <- 10 + 2*coord[,2] 
  param.scale <- 5 + 2 * coord[,1] 
  param.shape <- rep(0.2, (n.site + n.unobserved))
  for (i in 1:(n.site + n.unobserved))
    data[,i] <- frech2gev(data[,i], param.loc[i], param.scale[i], param.shape[i])
  
  ## Save the synthetic data
  file_name <- paste("data/synthetic_study/synthetic_data_", simulation_no, ".csv", sep = "")
  write.csv(data, file = file_name, row.names = FALSE)
  
  data_observed <- data[, 1:n.site]
  data_unobserved <- data[, (n.site + 1):(n.site + n.unobserved)]

  # FIT A MAX-STABLE MODEL TO SYNTHETIC DATA
    
  ## 1 - Fit a max-stable process with the following model for
  ## the GEV parameters
  form.loc <- y ~ lat 
  form.scale <- y ~ lon 
  form.shape <- y ~ 1
  brown_resnick <- fitmaxstab(data_observed, coord_observed, "brown", loc.form = form.loc, scale.form = form.scale, shape.form = form.shape)
  
  ## 2 - Predict the average at the 20 unobserved locations
  predictions <- predict(brown_resnick, coord_unobserved)
  loc.pred <- matrix(predictions[,"loc"], n.unobserved)
  scale.pred <- matrix(predictions[,"scale"], n.unobserved)
  shape.pred <- matrix(predictions[,"shape"], n.unobserved)
  
  n.sim = 10000
  sim <- rmaxstab(n.sim, coord_unobserved, "brown", nugget = brown_resnick$param["nugget"], range = brown_resnick$param["range"], smooth = brown_resnick$param["smooth"])
  for (i in 1:n.unobserved)
    sim[,i] <- frech2gev(sim[,i], loc.pred[i], scale.pred[i], shape.pred[i])
  
  sim_averages <- colMeans(sim)
  true_averages <- colMeans(data_unobserved)
  plot(sim_averages, true_averages)
  abline(a = 0, b = 1, col = "red")
  error <- true_averages - sim_averages
  # Compute MSE
  mse <- mean(error^2)
  rmse <- sqrt(mse)
  average_rmse[simulation_no] <- rmse
  
  ## 3 - Predict return levels
  
  ### 100-year
  max_stable_100_year_return_level_estimation <- numeric(length = ncol(sim))
  for (i in 1:ncol(sim)) {
    sorted_column <- sort(sim[, i], decreasing = TRUE)
    max_stable_100_year_return_level_estimation[i] <- sorted_column[100]
  }
  
  ### 10-year
  max_stable_10_year_return_level_estimation <- numeric(length = ncol(sim))
  for (i in 1:ncol(sim)) {
    sorted_column <- sort(sim[, i], decreasing = TRUE)
    max_stable_10_year_return_level_estimation[i] <- sorted_column[1000]
  }
  
  # FIT A GEV TO EACH STATION
  
  colnames(data) <- paste("col", 1:ncol(data), sep = "")
  
  gev_loc <- rep(0, n.site+n.unobserved)
  gev_scale <- rep(0, n.site+n.unobserved)
  gev_shape <- rep(0, n.site+n.unobserved)
  station_no <- rep(0, n.site+n.unobserved)
  return_level_100 <- rep(0, n.site+n.unobserved)
  return_level_10 <- rep(0, n.site+n.unobserved)
  true_average <- rep(0, n.site+n.unobserved)
  
  ## location, scale, shape in MLE
  for(i in 1:(ncol(data))) {
    fitting <- gev.fit(data[,(i)])
    station_number <- colnames(data)[i]
    loc <- fitting$mle[1]
    scale <- fitting$mle[2]
    shape <- fitting$mle[3]
    gev_loc[i] <- loc
    gev_scale[i] <- scale
    gev_shape[i] <- shape
    station_no[i] <- station_number
    return_level_100[i] <- loc + (scale/shape)*((-log(1 - 1/100))**(-shape) - 1)
    return_level_10[i] <- loc + (scale/shape)*((-log(1 - 1/10))**(-shape) - 1)
    true_average[i] <- mean(data[,(i)])
  }
  
  ## Combine vectors into a data frame
  df <- data.frame(station_no = station_no, 
                   gev_shape = gev_shape, 
                   gev_loc = gev_loc, 
                   gev_scale = gev_scale,
                   return_level_100 = return_level_100,
                   return_level_10 = return_level_10,
                   true_average = true_average)
  
  file_name <- paste("data/synthetic_study/synthetic_gev_", simulation_no, ".csv", sep = "")
  write.csv(df, file = file_name, row.names = FALSE)
  
  rownames(coord) <- paste("col", 1:nrow(coord), sep = "")
  file_name <- paste("data/synthetic_study/synthetic_coord_", simulation_no, ".csv", sep = "")
  write.csv(coord, file = file_name)
  
  
  # CHECK RMSE FOR RETURN LEVELS FOR MAX-STABLE
  
  ## 10-year return level
  true_values <- return_level_10[(length(return_level_10) - (n.unobserved - 1)):length(return_level_10)]
  simulated_values <- max_stable_10_year_return_level_estimation
  # Calculate the error
  error <- true_values - simulated_values
  # Compute MSE
  mse <- mean(error^2)
  # Calculate RMSE
  rmse <- sqrt(mse)
  return_level_10_rmse[simulation_no] <- rmse
  
  ## 100-year return level
  true_values <- return_level_100[(length(return_level_100) - (n.unobserved - 1)):length(return_level_100)]
  simulated_values <- max_stable_100_year_return_level_estimation
  # Calculate the error
  error <- true_values - simulated_values
  # Compute MSE
  mse <- mean(error^2)
  # Calculate RMSE
  rmse <- sqrt(mse)
  return_level_100_rmse[simulation_no] <- rmse
    
}

# See python code for evaluation of kriging (next subsection)
# Copy the output of these into the python file
print(paste(average_rmse, collapse = ", "))
print(paste(return_level_10_rmse, collapse = ", "))
print(paste(return_level_100_rmse, collapse = ", "))
