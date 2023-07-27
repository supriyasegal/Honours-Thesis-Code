# Instructions:
# - Run this code AFTER running 1-simulated-data-max-stable.R
# - Running this code in a Jupyter notebook is the easiest way to show the box-and-whisker plots
# - Read in the output from 1-simulated-data-max-stable.R into the lists "average_rmse", "return_level_10_rmse", "return_level_100_rmse"
average_rmse = [] # paste in output from R
return_level_10_rmse = [] # paste in output from R
return_level_100_rmse = [] # paste in output from R

# This code:
# a) Uses kriging to interpolate the GEV parameters of the synthetic data from 1-simulated-data-max-stable.R
# b) Evaluates both the max-stable and kriging models on the held-out stations in the synthetic data

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pykrige.kriging_tools as kt
from pykrige.ok import OrdinaryKriging
from pykrige.uk import UniversalKriging
from sklearn.metrics import mean_squared_error
import math

rmse_expected_extreme_rainfall = []
rmse_10_year_return_level = []
rmse_100_year_return_level = []

# Adjust parameters accordingly
n_sim = 100
n_site = 35
n_unobserved = 15

for simulation_no in range(1, n_sim+1):
    coord_df = pd.read_csv(f"data/synthetic_study/synthetic_coord_{simulation_no}.csv")
    coord_df = coord_df.rename(columns={'Unnamed: 0': 'station_no'})
    gev_df = pd.read_csv(f"data/synthetic_study/synthetic_gev_{simulation_no}.csv")
    stations_info_merged = coord_df.merge(gev_df, how="inner", on="station_no")
    
    # Define lambda function
    gev_mean_func = lambda row: row['gev_loc'] + row['gev_scale']*((math.gamma(1-row['gev_shape'])) - 1)/row['gev_shape']
    # Apply the lambda function to each row of the DataFrame
    stations_info_merged['gev_mean'] = stations_info_merged.apply(gev_mean_func, axis=1)
    
    training_data = stations_info_merged.iloc[:n_site]
    prediction_data = stations_info_merged.iloc[n_site:]
    
    # Training data
    x_train = np.array(training_data['lat'])
    y_train = np.array(training_data['lon'])
    
    # Prediction data
    x_pred = np.array(prediction_data['lat'])
    y_pred = np.array(prediction_data['lon'])
    
    # location
    regional_coords = [x_train, y_train]
    phi_train = np.array(training_data['gev_loc'])
    UK1 = UniversalKriging(
        x_train,
        y_train,
        phi_train,
        variogram_model='linear',
        drift_terms=['regional_linear'],
        specified_drift=regional_coords,
        verbose=False,
        enable_plotting=False,
        nlags=70
    )
    # Grid for prediction
    z_pred, ss_pred = UK1.execute("points", x_pred, y_pred)
    predicted_location = z_pred.flatten()
    
    # scale
    phi_train = np.array(training_data['gev_scale'])
    UK1 = UniversalKriging(
        x_train,
        y_train,
        phi_train,
        variogram_model='linear',
        drift_terms=['regional_linear'],
        specified_drift=regional_coords,
        verbose=False,
        enable_plotting=False,
        nlags=70
    )
    # Grid for prediction
    z_pred, ss_pred = UK1.execute("points", x_pred, y_pred)
    predicted_scale = z_pred.flatten()
    
    # shape
    phi_train = np.array(training_data['gev_shape'])
    UK1 = OrdinaryKriging(
        x_train,
        y_train,
        phi_train,
        variogram_model='linear',
        verbose=False,
        enable_plotting=False,
        nlags=70
    )
    # Grid for prediction
    z_pred, ss_pred = UK1.execute("points", x_pred, y_pred)
    predicted_shape = z_pred.flatten()

    # 1) EXPECTED EXTREME RAINFALL
    
    gev_mean_func2 = lambda loc, scale, shape: loc + scale * ((math.gamma(1 - shape)) - 1) / shape
    predicted_mean = list(map(gev_mean_func2, predicted_location, predicted_scale, predicted_shape))
    # Calculate RMSE
    true_gev_loc = np.array(prediction_data['true_average'])
    rmse = mean_squared_error(true_gev_loc, predicted_mean, squared=False)
    rmse_expected_extreme_rainfall.append(rmse)
    
    # 2) RETURN LEVELS 10 YEARS
    return_level_10_func = lambda loc, scale, shape: loc + (scale / shape) * ((-math.log(1 - 1/10)) ** (-shape) - 1)
    predicted_return_level_10 = list(map(return_level_10_func, predicted_location, predicted_scale, predicted_shape))
    true_gev_loc = np.array(prediction_data['return_level_10'])
    rmse = mean_squared_error(true_gev_loc, predicted_return_level_10, squared=False)
    rmse_10_year_return_level.append(rmse)
    
    # 3) RETURN LEVELS 100 YEARS
    return_level_100_func = lambda loc, scale, shape: loc + (scale / shape) * ((-math.log(1 - 1/100)) ** (-shape) - 1)
    predicted_return_level_100 = list(map(return_level_100_func, predicted_location, predicted_scale, predicted_shape))
    true_gev_loc = np.array(prediction_data['return_level_100'])
    rmse = mean_squared_error(true_gev_loc, predicted_return_level_100, squared=False)
    rmse_100_year_return_level.append(rmse)
    
    
# Box-and-whisker plots for expected extreme rainfall
data = [rmse_expected_extreme_rainfall, average_rmse]
labels = ['Kriging', 'Max-stable']
plt.boxplot(data, labels=labels)
plt.xlabel('Method')
plt.ylabel('RMSE')
plt.title(f"n.years = 200\nExpected extreme rainfall")
plt.show()

# Box-and-whisker plots for 10 year return levels
data = [rmse_10_year_return_level, return_level_10_rmse]
labels = ['Kriging', 'Max-stable']
plt.boxplot(data, labels=labels)
plt.xlabel('Method')
plt.ylabel('RMSE')
plt.title('n.years = 200\n10-year return level')
plt.show()

# Box-and-whisker plots for 100 year return levels
data = [rmse_100_year_return_level, return_level_100_rmse]
labels = ['Kriging', 'Max-stable']
plt.boxplot(data, labels=labels)
plt.xlabel('Method')
plt.ylabel('RMSE')
plt.title('n.years = 200\n100-year return level')
plt.show()
