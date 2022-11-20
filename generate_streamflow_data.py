import numpy as np
import pandas as pd

import xarray as xr
import geopandas as gpd

# From the PyNHD library, import data acuistion tools
import pydaymet as daymet
import pynhd as pynhd
from pynhd import NHD
from pynhd import NLDI, NHDPlusHR, WaterData
import pygeohydro as gh
from pygeohydro import NWIS, plot


#%%#############################################################################
### Step 0) Data and specifications
################################################################################
dates = ("1999-01-01", "2020-12-31")
dates = ("1999-01-01", "2010-12-31")
bbox = (-76.5, 37.5, -74.0, 44.0)
bbox = (-77.8, 37.5, -74.0, 44.0)  # Expanded on 11/2/22
nav_dist = 100

# Geo conversion
crs = 4386
crs_code = 'epsg:4386'

_ = xr.set_options(display_expand_attrs=False)

#%%#############################################################################
### 1. Get flow data
################################################################################
print("Getting flow data from USGS.")

# Initialize waterdata tool
wd = WaterData("nhdflowline_network")

# Use the national water info system (NWIS)
nwis = NWIS()

# Send a query for all gage info in the bbox
query = {"bBox": ",".join(f"{b:.06f}" for b in bbox),
        "hasDataTypeCd": "dv",
        "outputDataTypeCd": "dv"}

info_box = nwis.get_info(query)

# Filter stations to have only those with proper dates
stations = info_box[(info_box.begin_date <= dates[0]) & (info_box.end_date >= dates[1])].site_no.tolist()

## Source flows
print("Requesting data...")
qobs = nwis.get_streamflow(stations, dates, mmd=False)

# Remove gages with nans
qobs = qobs.dropna(axis = 1, how = 'any')
stations = qobs.columns.to_numpy()

print(f'All gage data sourced. You have {qobs.shape[1]} gages after cleaning.')


#%%#############################################################################
### 2. Get Location data (long,lat)
################################################################################
# Initialize storage
loc_data = np.zeros((len(qobs.columns), 2))
loc_data_names = np.array(['long', 'lat'])

# Loop through gage info queried earlier
for i,st in enumerate(qobs.columns):
    # Pull the lat-long
    long = info_box.set_index('site_no').loc[st.split('-')[1]]['dec_long_va']
    lat = info_box.set_index('site_no').loc[st.split('-')[1]]['dec_lat_va']
    # Store
    loc_data[i,0] = long if len(long.shape) == 0 else long[0]
    loc_data[i,1] = lat if len(lat.shape) == 0 else lat[0]

print("Long and lat data found.")
