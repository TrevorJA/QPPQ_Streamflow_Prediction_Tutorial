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


class StreamflowGenerator():
    def __init__(self, args, **kwargs):
        # Store variables
        self.observation_locations = args['observation_locations']
        self.historic_flows = args['historic_streamflows']
        self.K = args['K']

        # Define constants of interest
        self.n_observations = self.observation_locations.shape[0]
        self.T = self.historic_flows.shape[1]
        self.fdc_neps = np.linspace(0,1,200)


    def get_knn(self):
        """Find the distances and indices of the K nearest gage locations."""
        distances = np.zeros((self.n_observations))
        for i in range(self.n_observations):
            distances[i] = geodesic(self.prediction_location, self.observation_locations[i,:]).kilometers
        self.knn_indices = np.argsort(distances, axis = 0)[0:self.K].flatten()
        self.knn_distances = np.sort(distances, axis = 0)[0:self.K].flatten()
        return


    def calculate_nep(self, KNN_Q, Q_t):
        "Finds the NEP at time t based on historic observatons."
        # Get historic FDC
        fdc = np.quantile(KNN_Q, self.fdc_neps, axis = 1).T
        # Find nearest FDC value
        nearest_quantile = np.argsort(abs(fdc - Q_t), axis = 1)[:,0]
        nep = self.fdc_neps[nearest_quantile]
        return nep

    def interpolate_fdc(self, nep, fdc):
        """
        Performs linear interpolation of discrete FDC values to find flow at a NEP.

        Parameters
        ----------
        nep :: float
            Non-exceedance probability at a specific time.
        fdc :: array
            Array of discrete FDC points
        self.fdc_neps :: array
            Array of quantiles for discrete FDC.

        Returns
        -------
        flow :: float
            A single flow value given the NEP and FDC points.
        """
        tol = 0.0000001
        assert(len(fdc) == len(self.fdc_neps)), f'FDC and self.fdc_neps should be same length, but are {len(fdc)} and {len(self.fdc_neps)}.'
        if nep == 0:
            nep = np.array(tol)
        sq_diff = (self.fdc_neps - nep)**2

        # Index of nearest discrete NEP
        ind = np.argmin(sq_diff)

        # Handle edge-cases
        if nep <= self.fdc_neps[0]:
            return fdc[0]
        elif nep >= self.fdc_neps[-1]:
            return fdc[-1]

        if self.fdc_neps[ind] <= nep:
            flow_range = fdc[ind:ind+2]
            nep_range = self.fdc_neps[ind:ind+2]
        else:
            flow_range = fdc[ind-1:ind+1]
            nep_range = self.fdc_neps[ind-1:ind+1]

        slope = (flow_range[1] - flow_range[0])/(nep_range[1] - nep_range[0])
        flow = flow_range[0] + slope*(nep_range[1] - nep)
        return flow

    def predict_streamflow(self, args):
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
        self.prediction_location = args['prediction_location']
        self.prediction_fdc = args['prediction_fdc']
        self.fdc_quantiles = args['fdc_quantiles']
        self.n_predictions = self.prediction_location.shape[0]

        ### Find nearest K observations
        self.get_knn()
        knn_flows = self.historic_flows[self.knn_indices, :]

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
            self.predicted_flow[0,t] = self.interpolate_fdc(nep_t, self.prediction_fdc)

        return self.predicted_flow
