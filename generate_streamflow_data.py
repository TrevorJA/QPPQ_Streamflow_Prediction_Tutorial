"""
Trevor Amestoy
Cornell University
Fall 2022

Collect USGS streamflow data at all locations within a specified region and 
data range. 

Exports streamflowdata (cms) and location data (long, lat).

"""


import numpy as np
import pandas as pd

# From the PyNHD library, import data acuistion tools
from pygeohydro import NWIS


#%%#############################################################################
### Step 0) Data and specifications
################################################################################
# Specify time-range and region of interest
dates = ("2000-01-01", "2010-12-31")
region = (-108.0, 38.0, -105.0, 40.0)


#%%#############################################################################
### 1. Get flow data
################################################################################
print("Getting flow data from USGS.")

# Use the national water info system (NWIS)
nwis = NWIS()

# Send a query for all gage info in the region
query = {"bBox": ",".join(f"{b:.06f}" for b in region),
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
### 2. Get Location data (lat, long)
################################################################################
# Initialize storage
loc_data = np.zeros((len(qobs.columns), 2))

# Loop through gage info queried earlier
for i,st in enumerate(qobs.columns):
    # Pull the lat-long
    long = info_box.set_index('site_no').loc[st.split('-')[1]]['dec_long_va']
    lat = info_box.set_index('site_no').loc[st.split('-')[1]]['dec_lat_va']
    # Store
    loc_data[i,0] = lat if len(lat.shape) == 0 else lat[0]
    loc_data[i,1] = long if len(long.shape) == 0 else long[0]
    

print("Long and lat data found.")

#%%#############################################################################
### 3. Export data
################################################################################

np.savetxt('./data/observed_gage_locations.csv', loc_data, delimiter = ',')
np.savetxt('./data/observed_streamflow.csv', qobs.T, delimiter = ',')