from neuron import h
import numpy as np
import pickle
from .utils import *

def blender_dict_from_section(
                            section, ID,
                            voltage_data, 
                            frame_num):
    '''
    Creates an output dictionary for a single branch, containing 3D coordinates of Section points (N=n3d)
    and activity frame-by-frame data at each segment (N=nseg)

    Parameters:
        - section (nrn.Section) - NEURON section object (to retrieve 3D geometry)
        - ID (Int) - ID of a Section used to label branches in Blender (usually equal to an index in h.allseg)
        - voltage_data (Nseg x t) array_like segment-wise voltage data
        - frame_num (Int) - number of frames used in Blender animation (linear interpolation is performed to match that to voltage_data time)
    '''
    output = dict()
    output["ID"] = ID
    output["type"] = type_from_section_name(section.name())
    output["X"] = np.array([np.round(section.x3d(k),3) for k in range(section.n3d())])
    output["Y"] = np.array([np.round(section.y3d(k),3) for k in range(section.n3d())])
    output["Z"] = np.array([np.round(section.z3d(k),3) for k in range(section.n3d())])
    output["DIAM"] = np.array([np.round(section.diam3d(k),3) for k in range(section.n3d())])
    output["Voltage"] = dict()

    if voltage_data.shape[0] != section.nseg:
        raise Exception("Dimension of voltage_data and sec.nseg mismatch! ({} and {})".format(voltage_data.shape[0], section.nseg))

    interp_matrix = np.vstack([linear_interpolation(voltage_data[k,:],frame_num) for k in range(section.nseg)])
    for frame in range(frame_num):
        output["Voltage"][frame] = interp_matrix[:,frame]
    return output



class SectionMonitor():
    def __init__(self, section,dt=0.025):
        '''
            Class for recording activity of all segments in a NEURON Section

            dt – fixed sampling rate of a monitor
        '''
        self.nseg = section.nseg
        self.section = section

        # Creating NEURON Vectors
        self.Vectors = []
        for seg in list(section.allseg())[1:-1]:
            self.Vectors.append(h.Vector())
            self.Vectors[-1].record(seg._ref_v, dt)

class CellRecorder():
    '''
        Top-level class to record voltage at each segment and export the cell's geometry and animation data into Blender
    '''
    def __init__(self, sections,dt=0.025):
        '''
            Initialize a cell recorder

            sections (List) - list of NEURON sections in a cell. 
            dt (Float) fixed sampling rate of a monitor
        '''
        self.sections = sections
        self.dt=dt

        self.construct_monitors()

    def construct_monitors(self):
        self.monitors = []
        for sec in self.sections:
            self.monitors.append(SectionMonitor(
                sec, dt=self.dt
            ))
    
    def save_pickle(self, filename, FRAME_NUM=100, voltage_data=None):
        '''
            Save geometry and voltage into .pickle to bring it to Blender. 

            Arguments:
            - filename - path to the pickle file to be saved
            - FRAME_NUM - number of Blender frames to use for resulting animation
            - voltage_data (List of arrays. Default - None). 
                Section-wise voltage data. Specify only if you wish to use external simulation results.
                Otherwise, leave unspecified to use the recorded data from monitors
        '''

        if voltage_data is None:
            voltage_data = [matrix_from_monitor(mon) for mon in self.monitors]
        
        output = [
            blender_dict_from_section(
                self.sections[k],k,voltage_data[k],FRAME_NUM
            ) for k in range(len(self.sections))
        ]

        with open(filename, "wb") as f:
            pickle.dump(output,f)
