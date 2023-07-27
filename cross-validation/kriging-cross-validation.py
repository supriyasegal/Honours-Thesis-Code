from sklearn.metrics import mean_absolute_error, mean_squared_error
from pykrige.uk import UniversalKriging
from pykrige.ok import OrdinaryKriging

eau_stations = pd.read_csv("data/EAU_lon_lat_alt.csv")
eau_stations['station_no_int'] = eau_stations['station_no'].str[-8:].astype(int)
eau_station_numbers = list(eau_stations["station_no_int"])

sau_stations = pd.read_csv("data/SAU_lon_lat_alt.csv")
sau_stations['station_no_int'] = sau_stations['station_no'].str[-8:].astype(int)
sau_station_numbers = list(sau_stations["station_no_int"])

cau_stations = pd.read_csv("data/CAU_lon_lat_alt.csv")
cau_stations['station_no_int'] = cau_stations['station_no'].str[-8:].astype(int)
cau_station_numbers = list(cau_stations["station_no_int"])

stations_info_merged = pd.read_csv("data/gev_parameters_station_properties.csv")


### CHOOSE REGION:
selected_region = cau_station_numbers.copy()

### LOOCV
mae_list = []
mse_list = []
predicted_rain = []
true_rain = []

for i, row in stations_info_merged.iterrows():
    # Create a new DataFrame with one station removed
    if row['station_no'] not in selected_region:
        continue
    
    print(row['station_no'])
    
    temp_df = stations_info_merged.drop(i)
    
    # Extract the remaining data
    x = np.array(temp_df['Longitude'])
    y = np.array(temp_df['Latitude'])
    
    phi = np.array(temp_df['gev_loc'])
    OK_temp = UniversalKriging(
        x, 
        y, 
        phi, 
        variogram_model='gaussian',
        drift_terms=['regional_linear'],
        verbose=False,
        enable_plotting=False,
        nlags=70
        #coordinates_type='geographic'
    )
    predicted_location, ss_temp = OK_temp.execute('points', np.array([row['Longitude']]), np.array([row['Latitude']]))
    
    phi = np.array(temp_df['gev_scale'])
    OK_temp = UniversalKriging(
        x, 
        y, 
        phi, 
        variogram_model='gaussian',
        drift_terms=['regional_linear'],
        verbose=False,
        enable_plotting=False,
        nlags=70
        #coordinates_type='geographic'
    )
    predicted_scale, ss_temp = OK_temp.execute('points', np.array([row['Longitude']]), np.array([row['Latitude']]))
    
    phi = np.array(temp_df['gev_shape'])
    OK_temp = OrdinaryKriging(
        x,
        y,
        phi,
        variogram_model='linear',
        verbose=False,
        enable_plotting=False,
        nlags=70
    )
    predicted_shape, ss_temp = OK_temp.execute('points', np.array([row['Longitude']]), np.array([row['Latitude']]))
    
    
    gev_mean_func2 = lambda loc, scale, shape: loc + scale * ((math.gamma(1 - shape)) - 1) / shape
    predicted_mean = list(map(gev_mean_func2, predicted_location, predicted_scale, predicted_shape))
    
    # Compute the absolute error between the predicted and actual GEV mean
    mae = mean_absolute_error(np.array([row['avg_rainfall']]), predicted_mean)
    
    # Compute the ROOT mean squared error between the predicted and actual GEV mean
    mse = mean_squared_error(np.array([row['avg_rainfall']]), predicted_mean, squared=False)
    
    
    # store the prediction
    predicted_rain.append(predicted_mean)
    true_rain.append(np.array([row['avg_rainfall']]))
    
    # Add the MAE and MSE to the lists
    mae_list.append(mae)
    mse_list.append(mse)

# Compute the average MAE and RMSE over all stations
mean_mae = np.mean(mae_list)
mean_rmse = np.mean(mse_list)

print(mean_mae)
print(mean_rmse)

## Plot of predictions vs true

plt.plot(predicted_rain, true_rain, 'bo')
plt.plot([0, 140], [0, 140], 'r--')
# Set plot limits
plt.xlim(0, 110)
plt.ylim(0, 110)
# Set axis labels
plt.xlabel('predicted')
plt.ylabel('true')
plt.gcf().set_size_inches(7, 6)  # Adjust the size 
# Show the plot
plt.show()
