import sys
import os

from goniometer import GenericStateMatrixProvider

# Add the parent directory to the system path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from model import Mantex
import numpy as np
from copy import deepcopy
from mantex_run import GetAllPoleFigures, GetSinglePoleFigure
from helper_funcs import *

class NyrtexMantex():
    def __init__(self, source, detectors, sample,
                 exp_data, sample_view_axes, detector_colors, q_probe = 1, probe_window = 0.05, run_index = 0):
        self.source = source
        self.detectors = detectors
        self.sample = sample
        self.reference_sample = deepcopy(sample)
        self.reference_sample.smp = GenericStateMatrixProvider(np.eye(3), np.zeros(3))
        self.raw_exp_data = deepcopy(exp_data)
        self.sample_view_axes = sample_view_axes
        self.detector_colors = detector_colors
        self.q_probe = q_probe
        self.probe_window = probe_window
        self.run_index = run_index
        self.q_range = self.raw_exp_data.detector_readouts[0,0,:,0].copy()

        self.r_dict = {'recip_sample': 1,
                                   'lab_sample': 3,
                                   'k_vecs': 1,
                                   'gonio_r': 1,
                                   'gonio_v': 1}

        self.Alg = GetAllPoleFigures(deepcopy(self.raw_exp_data), self.source, deepcopy(self.sample),
                                     deepcopy(self.sample_view_axes), self.q_probe, self.probe_window)
        self.Alg.execute()
        # set to current run
        self.update_mantex()

    def update_mantex(self):
        new_sample = deepcopy(self.sample)
        new_sample.smp = self.raw_exp_data.smps[self.run_index]
        self.mantex = Mantex(self.source, self.detectors, new_sample, self.raw_exp_data.smps[self.run_index],
                             self.sample_view_axes, self.q_probe, self.probe_window)
        self.mantex.calc_pole()

    def get_current_readouts(self):
        dat = deepcopy(self.raw_exp_data.detector_readouts[self.run_index])
        #return dat
        return [self.mantex.correct_spectrum_for_attenuation(dat[idet], det) for idet, det in enumerate(self.detectors)]

    def get_detector_vectors(self):
        return [det.position - self.sample.position for det in self.detectors]

    def get_norm_detector_vectors(self):
        return [x/np.linalg.norm(x) for x in [det.position - self.sample.position for det in self.detectors]]

    def get_cube(self, size = (0.1, 0.1, 0.1),
                  offset = (0,0,0), use_gonio = True):
        if use_gonio == True:
            cube_data = cuboid_data(offset, size, self.sample.lab_frame_axes)
        else:
            cube_data = cuboid_data(offset, size, np.eye(3,3))
        return cube_data

    def get_lab_frame(self):
        source_repr = self.get_cube(offset = self.source.position, size=[2,1,0.5], use_gonio= False)
        source_col = ['grey']*3
        det_reprs = []
        det_cols = []
        for idet, det in enumerate(self.detectors):
            det_reprs.append(self.get_cube(offset=det.position, size=[0.1, 0.1, 0.1], use_gonio=False))
            det_cols.append([self.detector_colors[idet]] * 3)
        return source_repr, source_col, det_reprs, det_cols
    def get_mesh_size(self):
        return np.linalg.norm(self.sample.get_mesh_vectors()[:,0,:], axis = 1).max()





