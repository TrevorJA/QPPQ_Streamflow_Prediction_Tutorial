## Streamflow Predictions in Ungauged Basins

Predicting streamflow at ungauged locations is a classic problem in hydrology which has motivated significant research over the last several decades ([Hrachowitz et al., 2013](Hrachowitz, Markus, et al. "A decade of Predictions in Ungauged Basins (PUB)—a review." _Hydrological sciences journal_ 58.6 (2013): 1198-1255.)).  

There are numerous different methods for performing predictions in ungauged basins, but here I focus on the common *QPPQ method*.

Below, I describe the method and further down I provide a walkthrough demonstration of QPPQ streamflow prediction in Python.  

The supporting code can be found on my GitHub here: [QPPQ_Streamflow_Prediction_Tutorial](https://github.com/TrevorJA/QPPQ_Streamflow_Prediction_Tutorial).

### QPPQ-Method for Streamflow Prediction

[Fennessey (1994)](https://scholar.google.com/citations?view_op=view_citation&hl=en&user=efFhgZ8AAAAJ&citation_for_view=efFhgZ8AAAAJ:zYLM7Y9cAGgC) introduced the *QPPQ method* for streamflow estimation at ungauged locations. 

The QPPQ method is commonly used and encouraged by the USGS, and is described at length in their publication [*Estimation of Daily Mean Streamflow for Ungaged Stream locations...* (2016)](https://pubs.usgs.gov/sir/2015/5157/sir20155157.pdf).

**QPPQ consists of four key steps:**
1. Estimating an FDC for the target catchment of interest, $\hat{FDC}_{pred}$.
2. Identify $K$ donor locations, nearest to the target point.
3. Transferring the timeseries of nonexceedance probabilities ($\mathbf{P}$) from the donor site(s) to the target.
4. Using estimated FDC for the target to map the donated nonexceedance timeseries, $\mathbf{P}$ back to streamflow.

To limit the scope of this tutorial, let's assume that an estimate of the FDC at the target site, $\hat{FDC}_{pred}$, has already been determined through some other statistical or observational study. 

Then the QPPQ method can be described more formally. Given an ungauged location with an estimated FDC, $\hat{FDC}_{pred}$, and set of observed streamflow timeseries $\mathbf{q_i}$ at $K$ neighboring sites, such that:

$Q_{obs} = \set{\mathbf{q_1}, \mathbf{q_2}, ..., \mathbf{q_k}}$

With corresponding $K$ FDCs at the observation locations:

$FDC_{obs} = \set{FDC_1, FDC_2, ... , FDC_k}$

The FDCs are used to convert the observed streamflow timeseries, $\mathbf{q_{obs, i}}$, to non-exceedance probability timeseries, $\mathbf{p_{obs, i}}$.

$FDC_i : \mathbf{q_{i}} \to \mathbf{p_i}$

We can then perform a weighted-aggregation of the non-exceedance probability timeseries to estimate the non-exceedance timeseries at the ungauged location. It is most common to apply an inverse-squared-distance weight to each observed timeseries such that:

$\mathbf{p_{pred}} = \sum^k (\mathbf{p_i}w_i)$

Where $w_i = 1 / d_i^2$ where $d_i$ is the distance from the observation $i$ to the ungauged location, and $\sum^k w_i = 1$. 

Finally, the estimated FDC at the ungauged location, $\hat{FDC}_{pred}$, is used to convert the non-exceedance timeseries to streamflow timeseries:

$\hat{FDC_{pred}} : \mathbf{p_{pred}} \to \mathbf{q_{pred}}$

Looking at this formulation, and the sequence of transformations that take place, I hope it is clear why the method is rightfully called the *QPPQ method*.

This method is summarized well by the taken from the [USGS Report on the topic](https://pubs.usgs.gov/sir/2015/5157/sir20155157.pdf):

![QPPQ method graphic](https://github.com/TrevorJA/QPPQ_Streamflow_Prediction_Tutorial/blob/main/QPPQ_Method_Graphic.PNG)

In the following section, I step through an implementation of this method in Python.

## Tutorial

All Python scripts used in this tutorial can be found in my GitHub repository: [QPPQ_Streamflow_Prediction_Tutorial](https://github.com/TrevorJA/QPPQ_Streamflow_Prediction_Tutorial).

In order run the scripts in this tutorial yourself, you will need to have installed the Python libraries listed in `requirements.txt`, including:
```
numpy
pandas
matplotlib
sklearn
geopy
geopandas
pygeohydro   # Only if generating your own data
```
Running `pip install -r requirements.txt` from the command line, while inside a local copy of the directory will install all of these packages.

#### Data retrieval

I collected USGS streamflow data from $N$ gages using the [[HyRiver]]suite for Python. 

If you would like to learn more about hydro-environmental data acquisition in Python, check out my old post on [*Efficient hydroclimatic data accessing with HyRiver for Python*](https://waterprogramming.wordpress.com/2022/09/20/efficient-hydroclimatic-data-accessing-with-hyriver-for-python/).

The script used to retrieve the data is available [here](https://github.com/TrevorJA/QPPQ_Streamflow_Prediction_Tutorial/blob/main/generate_streamflow_data.py). If you would like to experiment with this method in other regions, you can change the `region` variable on line 21, which specifies the corners of a bounding-box within which gage data will be retrieved:

```python
# Specify time-range and region of interest
dates = ("2000-01-01", "2010-12-31")
region = (-108.0, 38.0, -105.0, 40.0)
```
Above, I specify a region West of the Rocky Mountains in Colorado. Running the `generate_streamflow_data.py`, I found 73 USGS gage locations (blue circles). 


![Gage locations plotted on a map.](https://github.com/TrevorJA/QPPQ_Streamflow_Prediction_Tutorial/blob/main/gage_location_map.PNG)

#### QPPQ Model

The file `QPPQ.py` contains the method outlined above, defined as the `StreamflowGenerator` class object. 

The `StreamflowGenerator` has four key methods (or functions):

```python
class StreamflowGenerator():
    def __init__(self, args):  
	    ...
	def get_knn(self):
		...
	def calculate_nep(self):
		...
	def interpolate_fdc(self):
		...
	def predict_streamflow(self):
		...
		return predicted_flow
```

The method `get_knn` finds the $k$ observation, gage locations nearest to the prediction location, and stores the distances to these observation locations (`self.knn_distances`) and the indices associated with these locations (`self.knn_indices`).

```python
    def get_knn(self):
        """Find distances and indices of the K nearest gage locations."""
        distances = np.zeros((self.n_observations))
        for i in range(self.n_observations):
            distances[i] = geodesic(self.prediction_location, self.observation_locations[i,:]).kilometers
        self.knn_indices = np.argsort(distances, axis = 0)[0:self.K].flatten()
        self.knn_distances = np.sort(distances, axis = 0)[0:self.K].flatten()
        return
```

The next method, `calculate_nep`, calculates the NEP of a flow at an observation location at time $t$, or $P(Q \leq q_t)_{i}$.

```python 
    def calculate_nep(self, KNN_Q, Q_t):
        "Finds the NEP at time t based on historic observatons."
        # Get historic FDC
        fdc = np.quantile(KNN_Q, self.fdc_neps, axis = 1).T
        # Find nearest FDC value
        nearest_quantile = np.argsort(abs(fdc - Q_t), axis = 1)[:,0]
        nep = self.fdc_neps[nearest_quantile]
        return nep	
```

The `interpolate_fdc` performs a linear interpolate of the discrete FDC, and estimates flow for some given NEP. 

```python
    def interpolate_fdc(self, nep, fdc):
        "Performs linear interpolation of discrete FDC at a NEP."
        tol = 0.0000001
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
```

Finally, `predict_streamflow(self, *args)` combines these other methods and performs the full QPPQ prediction. 

```python
    def predict_streamflow(self, args):
        "Run the QPPQ prediction method for a single locations."
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
```

The `predict_streamflow` method is the only function called directly by the user.  While `get_knn`, `calculate_nep`, and `interpolate_fdc` are all used by `predict_streamflow`.

#### Generate streamflow predictions

The script [`run_QPPQ_predictions.py`](https://github.com/TrevorJA/QPPQ_Streamflow_Prediction_Tutorial/blob/main/run_QPPQ_prediction.py) runs the model and produces predictions at a test site. First, the data generated by [`generate_streamflow_data.py`](https://github.com/TrevorJA/QPPQ_Streamflow_Prediction_Tutorial/blob/main/generate_streamflow_data.py) is loaded:

```python
import numpy as np
from QPPQ import StreamflowGenerator

### Load Data
gage_locations = np.loadtxt('./data/observed_gage_locations.csv', delimiter = ',')
observed_flows = np.loadtxt('./data/observed_streamflow.csv', delimiter = ',')
```

The FDCs at each site are estimated at 200 discrete quantiles:
```python
fdc_quantiles = np.linspace(0,1,200)
observed_fdc = np.quantile(observed_flows, fdc_quantiles, axis =1).T
```
A random test site is selected, and removed from the observation data:
```python
# Select a test site and remove it from observations
test_site = np.random.randint(0, gage_locations.shape[0])

# Store test site
test_location = gage_locations[test_site,:]
test_flow = observed_flows[test_site, :]
test_site_fdc = observed_fdc[test_site, :]

# Remove test site from observations
gage_locations = np.delete(gage_locations, test_site, axis = 0)
observed_flows = np.delete(observed_flows, test_site, axis = 0)
```

When initializing the `StreamflowGenerator`, we must provide an array of gage location data (longitude, latitude), historic streamflow data at each gage, and the $K$ number of nearest neighbors to include in the timeseries aggregation. 

I've set-up the `StreamflowGenerator` class to receive these inputs as a dictionary, such as:

```python
# Specify model prediction_inputs
QPPQ_args = {'observed_locations' : gage_locations,
            'historic_flows' : observed_flows,
            'K' : 20}

# Intialize the model
QPPQ_model = StreamflowGenerator(QPPQ_args)
```

Similarly, the prediction arguments are provided as a dictionary to the `predict_streamflow` function:

```python 
# Specify the prediction arguments
prediction_args = {'prediction_location': test_location,
                    'prediction_fdc' : test_site_fdc,
                    'fdc_quantiles' : fdc_quantiles}
                    
# Run the prediction
predicted_flow = QPPQ_model.predict_streamflow(prediction_args)
```
I made a function, `plot_predicted_and_observed`, which allows for a quick visual check of the predicted timeseries compared to the observed timeseries:

```python
from plot_functions import plot_predicted_and_observed
plot_predicted_and_observed(predicted_flow, test_flow)
```
Which shows some very-nice quality predictions!

![Comparison of prediction and actual streamflow.](https://github.com/TrevorJA/QPPQ_Streamflow_Prediction_Tutorial/blob/main/output/streamflow_prediction_fig.png)

One benefit of working with the `StreamflowGenerator` as a Python `class` object is that we can retrieve the internal variables for further inspection.  

For example, I can call `QPPQ_model.knn_distances` to retrieve the distances to the $k$ nearest neighbors used to predict the flow at the ungauged location.  In this case, the gages used to make the prediction above were located $4.44, 13.23,. 18.38,...$ kilometers away.  

## Caveat and Conclusion

It is worth highlighting one major caveat to this example, which is that the FDC used for the prediction site was perfectly known from the historic record. In most cases, the FDC will not be known when making predictions in ungauged basins.  Rather, estimations of the FDC will need to be used, and thus the prediction quality shown above is somewhat of a *ideal-case* when performing a QPPQ in ungauged basins.  

There are numerous methods for estimating FDCs at the ungauged site, including the [Generalized Pareto distribution approximation proposed by Fennessey (1994)](https://scholar.google.com/citations?view_op=view_citation&hl=en&user=efFhgZ8AAAAJ&citation_for_view=efFhgZ8AAAAJ:zYLM7Y9cAGgC) or, more recently, through the use of Neural Networks, as highlighted in [Worland,  et al. (2019)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2018WR024463).

Hopefully this tutorial helped to get you familiar with a foundational streamflow prediction method.  

$\hat{FDC}_{pred} : \mathbf{p_{pred}} \to \mathbf{q_{pred}}, \, \mathbf{P} = \set{\mathbf{p}_1, …, \mathbf{p}_k}$


## References

Fennessey, Neil Merrick. "A Hydro-Climatological Model of Daily Stream Flow for the Northeast United States." Order No. 9510802 Tufts University, 1994. Ann Arbor: _ProQuest._ Web. 21 Nov. 2022.

Hrachowitz, Markus, et al. "A decade of Predictions in Ungauged Basins (PUB)—a review." _Hydrological sciences journal_ 58.6 (2013): 1198-1255.

Razavi, Tara, and Paulin Coulibaly. "Streamflow prediction in ungauged basins: review of regionalization methods." _Journal of hydrologic engineering_ 18.8 (2013): 958-975.

Stuckey, M.H., 2016, Estimation of daily mean streamflow for ungaged stream locations in the Delaware River Basin, water years 1960–2010: U.S. Geological Survey Scientific Investigations Report 2015–5157, 42 p., http://dx.doi.org/10.3133/sir20155157.

Worland, Scott C., et al. "Prediction and inference of flow duration curves using multioutput neural networks." _Water Resources Research_ 55.8 (2019): 6850-6868.