library(SpatialExtremes)
library(fields)

# Read in data
# Adjust to EAU, CAU, or SAU
rain <- read.table("data/CAU_max_rainfall.csv", sep=",", header=TRUE)
rain <- rain[ -c(1) ]
rain <- data.frame(rain[,-1], row.names = rain[,1])
rain <- data.matrix(rain)

coord_df <- read.table("data/CAU_lon_lat_alt.csv", sep=",", header=TRUE)
coord_df <- coord_df[ -c(1,2) ]
coord <- coord_df[ -c(3) ] # remove altitude
coord <- data.matrix(coord)

loc.form <- y ~ lon.scaled + lat.scaled
scale.form <- y ~ lon.scaled + lat.scaled
shape.form <- y ~ 1

# get true means for each station
true_means <- colMeans(rain)

# initialize vector for simulated averages
sim_averages <- rep(0, ncol(rain))

i <- 1
# perform leave-one-out cross-validation
for (i in 1:ncol(rain)) {
  # remove current station from data
  rain_cv <- rain[,-i]
  coord_cv <- coord[-i,]
  #covariates_cv <- covariates[-i,]
  
  ## Center and scale the covariates
  covariates_cv <- scale(coord_cv)
  colnames(covariates_cv) <- c("lon.scaled", "lat.scaled")
  
  # fit model without current station
  M0 <- fitmaxstab(rain_cv, coord_cv, "twhitmat", nugget = 0, loc.form, scale.form, shape.form, marg.cov = covariates_cv)
  
  ## Since we scale the covariates we need to do the same with our left-out station
  coord.scaled <- coord
  colnames(coord.scaled) <- c("lon.scaled", "lat.scaled")
  coord.scaled <- coord.scaled[i,]
  for (j in 1:2)
    coord.scaled[j] <- (coord.scaled[j] - attributes(covariates_cv)$"scaled:center"[j]) / attributes(covariates_cv)$"scaled:scale"[j]
  
  # get prediction for current station
  current_station_pred <- predict(M0, coord.scaled)
  loc.pred_cv <- current_station_pred[,"loc"]
  scale.pred_cv <- current_station_pred[,"scale"]
  shape.pred_cv <- current_station_pred[,"shape"]
  
  # simulate from fitted model
  n.sim = 500
  sim <- rmaxstab(n.sim, coord[i, , drop = FALSE], "twhitmat", DoF = M0$param["DoF"], nugget = M0$param["nugget"], range = M0$param["range"], smooth = M0$param["smooth"])
  
  # switch to appropriate GEV margins
  for (s in 1:n.sim)
    sim[s] <- frech2gev(sim[s], loc.pred_cv, scale.pred_cv, shape.pred_cv)
  
  # get average at current station
  sim_averages[i] <- mean(sim)
  
}

plot(sim_averages, true_means, xlim = c(0, max(sim_averages, true_means)),
     ylim = c(0, max(sim_averages, true_means)), xlab = "predicted", ylab = "true")
abline(a = 0, b = 1, col = "red")

rmse_list <- c()
# Loop through the arrays pairwise
for (i in 1:length(true_means)) {
  # Compute RMSE and MAE between the ith elements
  rmse <- sqrt(mean((true_means[i] - sim_averages[i])^2))
  mae <- abs(true_means[i] - sim_averages[i])
  # Add RMSE and MAE to their respective lists
  rmse_list <- c(rmse_list, rmse)
}

# Calculate the average RMSE and MAE
avg_rmse <- mean(rmse_list)

# Print the results
print(paste("RMSE:", avg_rmse))

# Create a data frame
data <- data.frame(sim_averages, true_means)
# Save to CSV
write.csv(data, file = "CAU_cross_validation_max_stable.csv", row.names = FALSE)
