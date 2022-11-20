"""
Trevor Amestoy
Cornell University
Fall 2022

The function IDW_generator takes:
- A set of observed flow timeseries
- FDC at a target site
- The locations (lat, long) of the observed sites
- The location (lat, long) of the target site


And performs an inverse-distance weighted prediction of flow timeseries at the
target site based on flow at `K` surrounding sites.
"""


import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler
from geopy.distance import geodesic

from utils import interpolate_FDC


class StreamflowGenerator():
    def __init__(self, args, **kwargs):
        # Store variables
        self.observation_locations = args['observation_locations']
        self.historic_flows = args['historic_streamflows']
        self.K = args['K']

        # Define constants of interest
        self.n_observations = self.observation_locations.shape[0]
        self.T = self.historic_flows.shape[1]


    def get_knn(self):
        """Find the distances and indices of the K nearest gage locations."""
        distances = np.zeros((self.n_observations))
        for i in range(self.n_observations):
            distances[i] = geodesic(self.prediction_location, self.observation_locations[i,:]).kilometers
        self.knn_distances = np.argsort(distances, axis = 0)[0:self.K].flatten()
        self.knn_distances = np.sort(distances, axis = 0)[0:self.K].flatten()
        return


    def calculate_nep(self, KNN_Q, Q_t):
        "Finds the NEP at time t based on historic observatons."
        # Get historic FDC
        quants = np.linspace(0,1,200)
        fdc = np.quantile(KNN_Q, quants, axis = 1).T
        # Find nearest FDC value
        nearest_quantile = np.argsort(abs(fdc - Q_t), axis = 1)[:,0]
        nep = quants[nearest_quantile]
        return nep


    def predict_streamflow(self, *args):
        """
        Run the QPPQ prediction method for a single locations.

        Parameters:
        ----------
        prediction_inputs : ndarray
            The matrix of feature values for the location of interest.
        predicted_fdc : ndarray
            An array of discrete FDC values.

        Returns:
        --------
        predicted_flow : ndarray
            A timeseries of predicted streamflow at the location.
        """
        self.prediction_location, self.prediction_fdc, self.fdc_quantiles = args
        self.n_predictions = self.prediction_location.shape[0]

        ### Find nearest K observations
        self.get_knn()
        knn_flows = self.historic_flows[self.knn_distances, :]

        ### Calculate weights as inverse square distance
        self.wts = 1/self.knn_distances**2

        # Normalize weights
        self.norm_wts = self.wts/np.sum(self.wts)

        ### Create timeseries of NEP at observation locations
        self.observed_neps = np.zeros_like(knn_flows)
        for t in range(self.T):
            self.observed_neps[:,t] = self.calculate_nep(knn_flows, knn_flows[:,t:t+1])

        ### Calculate predicted NEP timeseries using weights
        self.predicted_nep = np.zeros((self.n_predictions, self.T))
        for t in range(self.T):
            self.predicted_nep[:,t] = np.sum(self.observed_neps[:,t:t+1].T * self.norm_wts)

        ### Convert NEP timeseries to flow timeseries
        self.predicted_flow = np.zeros_like(self.predicted_nep)
        for t in range(self.T):
            nep_t = self.predicted_nep[0,:][t]
            self.predicted_flow[0,t] = interpolate_FDC(nep_t, self.prediction_fdc, self.fdc_quantiles)

        return self.predicted_flow