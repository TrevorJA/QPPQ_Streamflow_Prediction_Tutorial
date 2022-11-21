"""
Trevor Amestoy
Cornell University
Fall 2022


Run the streamflow prediction at ungauged locations.

"""

import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
from QPPQ import StreamflowGenerator

from plot_functions import plot_predicted_and_observed

### Load Data
gage_locations = np.loadtxt('./data/observed_gage_locations.csv', delimiter = ',')
observed_flows = np.loadtxt('./data/observed_streamflow.csv', delimiter = ',')


### Calculate FDCs
fdc_quantiles = np.linspace(0,1,200)

observed_fdc = np.quantile(observed_flows, fdc_quantiles, axis =1).T


### Choose one location as a test site
test_site = np.random.randint(0, gage_locations.shape[0])

test_location = gage_locations[test_site,:].reshape(1,-1)
test_flow = observed_flows[test_site, :]
test_site_fdc = observed_fdc[test_site, :]

gage_locations = np.delete(gage_locations, test_site, axis = 0)
observed_flows = np.delete(observed_flows, test_site, axis = 0)
observed_fdc = np.delete(observed_fdc, test_site, axis = 0)


# Specify model prediction_inputs
QPPQ_args = {'observation_locations' : gage_locations,
            'historic_streamflows' : observed_flows,
            'K' : 20}

# Initialize the model
QPPQ_model = StreamflowGenerator(QPPQ_args)


# Specify the prediction arguments
prediction_args = {'prediction_location': test_location,
                    'prediction_fdc' : test_site_fdc,
                    'fdc_quantiles' : fdc_quantiles}

# Run the prediction
predicted_flow = QPPQ_model.predict_streamflow(prediction_args)


### VISUALIZE

p1 = plot_predicted_and_observed(predicted_flow.flatten(), test_flow.flatten(), yscale = 'linear')

"""
### Visualize gage locations on a map

import folium
latlon = gage_locations
mapit = folium.Map( location=[39, -107], zoom_start=6 )
for coord in latlon:
    folium.CircleMarker( location=[ coord[0], coord[1] ], fill_color='#43d9de', radius=12 ).add_to( mapit )
mapit.save( 'map.html')
"""


