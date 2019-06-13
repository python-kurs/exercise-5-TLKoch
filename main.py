# Exercise 5
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

input_dir  = Path("data")
output_dir = Path("solution")

# 1. Go to http://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php#datafiles
#    and download the 0.25 deg. file for daily mean temperature.
#    Save the file into the data directory but don't commit it to github!!! [2P]


# 2. Read the file using xarray. Get to know your data. What's in the file?
#    Calculate monthly means for the reference periode 1981-2010 for Europe (Extent: Lon_min:-13, Lon_max: 25, Lat_min: 30, Lat_max: 72). [2P]

tg_ds = xr.open_dataset(input_dir / "tg_ens_mean_0.25deg_reg_v19.0e.nc")

tg_mean_monthly = tg_ds.sel(latitude = slice(30, 72), longitude = slice(-13, 25), time = slice("1981", "2010")).groupby("time.month").mean("time")


# 3. Calculate monthly anomalies from the reference period for the year 2018 (use the same extent as in #2).
#    Make a quick plot of the anomalies for the region. [2P]

tg_2018 = tg_ds.sel(latitude = slice(30, 72), longitude = slice(-13, 25), time = slice("2018", "2018")).groupby("time.month").mean("time")

tg_anomalies_2018 = tg_2018 - tg_mean_monthly 

tg_anomalies_2018["tg"].plot()


# 4. Calculate the mean anomaly for the year 2018 for Europe (over all pixels of the extent from #2) 
#    Compare this overall mean anomaly to the anomaly of the pixel which contains Marburg. 
#    Is the anomaly of Marburg lower or higher than the one for Europe? [2P] 

tg_mean_all = tg_ds.sel(latitude = slice(30, 72), longitude = slice(-13, 25), time = slice("1981", "2010")).tg.mean("time")

tg_mean_2018 = tg_ds.sel(latitude = slice(30, 72), longitude = slice(-13, 25), time = slice("2018", "2018")).tg.mean("time")

tg_mean_anomaly_2018 = tg_mean_2018 - tg_mean_all

tg_mean_anomaly_2018_allpixels = tg_mean_anomaly_2018.mean().values


# Marburg: lat: 50.81 lon: 8. 77
tg_mean_anomalies_2018_mr = tg_mean_anomaly_2018.sel(latitude=50.81, longitude = 8.77, method="nearest").values


if tg_mean_anomaly_2018_allpixels > tg_mean_anomalies_2018_mr:
    print("The anomaly of the mean temperature in 2018 is higher in whole Europa than in Marburg")
else:
    print("The anomaly of the mean temperature in 2018 is higher in Marburg than in whole Europe")


# 5. Write the monthly anomalies from task 3 to a netcdf file with name "europe_anom_2018.nc" to the solution directory.
#    Write the monthly anomalies for Marburg to a csv file with name "marburg_anom_2018.csv" to the solution directory. [2P]


tg_anomalies_2018.to_netcdf(output_dir / "europe_anom_2018.nc")


tg_monthly_mean_anomalies_2018_mr = tg_anomalies_2018.sel(latitude=50.81, longitude = 8.77, method="nearest").to_dataframe()

tg_monthly_mean_anomalies_2018_mr = tg_monthly_mean_anomalies_2018_mr.rename({'tg': 'anomaly'}, axis='columns')

tg_monthly_mean_anomalies_2018_mr.to_csv(output_dir / "marburg_anom_2018.csv", header = True)