"""
Trevor Amestoy
Cornell University
Fall 2022


Run the streamflow prediction at ungauged locations.

"""

import numpy as np

from QPPQ import StreamflowGenerator

### Load Data

gage_locations = np.loadtxt()
observed_flows = np.loadtxt()

### Choose one location as a test site
test_site = np.random.randint(0, gage_locations.shape[0])

test_location = gage_locations[test_site,:]
test_flow = observed_flows[test_site, :]

gage_locations = np.delete(gage_locations, test_site, axis = 0)
observed_flows = np.delete(observed_flows, test_site, axis = 0)


# Specify model prediction_inputs
QPPQ_args = {'observed_locations' : gage_locations,
            'historic_flows' : observed_flows,
            'K' : 20}

# Initialize the model
QPPQ_model = StreamflowGenerator(QPPQ_args)


# Specify the prediction arguments
prediction_args = {'prediction_location': test_location,
                    'prediction_fdc' : test_site_fdc,
                    'fdc_quantiles' : fdc_quantiles}

# Run the prediction
predicted_flow = QPPQ_model.predict_streamflow(prediction_args)
