"""
Trevor Amestoy
Cornell University 
Fall 2022

Simple plot functions.

"""

import matplotlib.pyplot as plt
import numpy as np

import geopandas as gpd

def plot_predicted_and_observed(P, Q, yscale = 'linear'):
	assert(len(P) == len(Q)), 'Data are not the same size.'
	xs = np.arange(len(P))
	plt.figure(figsize = (8,3), dpi =150)
	plt.plot(xs, P, label = 'Predicted Flow', color = 'blue', alpha = 0.8)
	plt.plot(xs, Q, label = 'Observed Flow', color = 'grey', alpha = 0.7, linewidth = 2)
	if yscale == 'log':
		plt.ylabel('Log-Flow (cms)')
	else:
		plt.ylabel('Flow (cms)')
	plt.xlabel('Day')
	plt.yscale(yscale)
	plt.legend()
	plt.show()
	return plt