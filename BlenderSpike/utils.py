import numpy as np
from scipy.interpolate import interp1d

def linear_interpolation(source_data, n_points):
    '''
        Linearly resamples an array with n_points
    '''

    if(len(source_data)==1):
        return np.ones(n_points)*source_data[0]
    source_prop = np.linspace(0,1, len(source_data))
    output_prop = np.linspace(0,1, n_points)
    interF = interp1d(source_prop, source_data)
    output= interF(output_prop)
    
    return output