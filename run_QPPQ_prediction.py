"""
Trevor Amestoy
Cornell University
Fall 2022


Run the streamflow prediction at ungauged locations.

"""

import numpy as np

from QPPQ import StreamflowGenerator

QPPQ_args = {'historic_flows' : observed_flows,
            'observed_locations' : observed_locations}
