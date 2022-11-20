import numpy as np
import pandas as pd

def find_continuous_FDC(data, quants):
    """
    Calculates FDC values as specified `quants` or non-exceedance probabilites.
    """
    flows = np.quantile(data, quants, axis = 1).T
    return flows


################################################################################

def find_NEPs(hist_data, flow):
    """
    Finds the non-exceedance probability (NEP) for `flow` given hist_data
    timeseries.
    """

    quants = np.linspace(0,1,200)
    fdc_data = find_continuous_FDC(hist_data, quants)

    # Find nearest FDC value
    diff = fdc_data - flow

    small_diff = np.argsort(abs(diff), axis = 1)[:,0]
    nep = quants[small_diff]
    return nep

################################################################################

def interpolate_FDC(nep, fdc, quants):
    """
    Performs linear interpolation of discrete FDC values to find flow at a NEP.

    Parameters
    ----------
    nep :: float
        Non-exceedance probability at a specific time.
    fdc :: array
        Array of discrete FDC points
    quants :: array
        Array of quantiles for discrete FDC.

    Returns
    -------
    flow :: float
        A single flow value given the NEP and FDC points.
    """
    tol = 0.0000001
    assert(len(fdc) == len(quants)), f'FDC and quants should be same length, but are {len(fdc)} and {len(quants)}.'
    if nep == 0:
        nep = np.array(tol)
    sq_diff = (quants - nep)**2

    # Index of nearest discrete NEP
    ind = np.argmin(sq_diff)

    # Handle edge-cases
    if nep <= quants[0]:
        return fdc[0]
    elif nep >= quants[-1]:
        return fdc[-1]

    if quants[ind] <= nep:
        flow_range = fdc[ind:ind+2]
        nep_range = quants[ind:ind+2]
    else:
        flow_range = fdc[ind-1:ind+1]
        nep_range = quants[ind-1:ind+1]


    #print(f'flow_range: {flow_range}')
    #print(f'nep_range: {nep_range}')
    slope = (flow_range[1] - flow_range[0])/(nep_range[1] - nep_range[0])
    flow = flow_range[0] + slope*(nep_range[1] - nep)

    return flow
