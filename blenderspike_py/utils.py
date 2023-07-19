import numpy as np
from scipy.interpolate import interp1d
import pickle

def type_from_section_name(section_name):
    if "apic" in section_name:
        return "apic"
    if "dend" in section_name:
        return "dend"
    if "axon" in section_name or "hillock" in section_name or "initial" in section_name:
        return "axon"
    if "soma" in section_name:
        return "soma"
    return section_name

def linear_interpolation(source_data, n_points):
    '''
        Linearly resamples an array with n_points
    '''
    source_prop = np.linspace(0,1, len(source_data))
    output_prop = np.linspace(0,1, n_points)
    interF = interp1d(source_prop, source_data)
    output= interF(output_prop)
    
    return output

def matrix_from_monitor(monitor):
    '''
        Returns an NsegxT array of section's voltage in each segment (Nseg) at each time point
    '''
    return np.vstack([np.array(vec) for vec in monitor.Vectors])

