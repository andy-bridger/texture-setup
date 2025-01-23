import sys
import os
# Add the parent directory to the system path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from mantex_sim import MantexSim
import numpy as np

class NyrtexMantex(MantexSim):
    def __init__(self, pf_ax, source, detectors, sample,
                 exp_data, sample_view_axes, detector_colors, q_probe = 1, probe_window = 0.05, run_index = 0):
        self.exp_data = exp_data
        self.pf_ax = pf_ax
        self.run = run_index
        self.q_range = exp_data.detector_readouts[0, 0, :, 0]
        ewald_radii = (np.min(self.q_range), np.max(self.q_range))
        super().__init__(source, detectors, sample, exp_data.smps[run_index], ewald_radii,
                 sample_view_axes, detector_colors, ewald_steps = 2, q_probe = q_probe, probe_window = probe_window, exp_data=exp_data)
        self.q_range = exp_data.detector_readouts[0, 0, :, 0]
    def update_mantex(self, run_index):
        super().__init__(self.source, self.detectors, self.sample,
                                  self.exp_data.smps[run_index], self.sample_view_axes, self.detector_colors, 2,
                                  self.q_probe, self.probe_window, exp_data=self.exp_data)

    def setup_qrange_rl(self):
        pass
    def get_inplane_pole_figure(self):
        pass
    def get_detector_signal(self, thresh = 0.1):
        q_space_readouts = self.exp_data.detector_readouts[self.run]
        self.q_range = self.exp_data.detector_readouts[0, 0, :, 0]
        self.detector_readout = []
        for idet, det in enumerate(self.detectors):
            ipl, entry_point, exit_point = self.sample.get_internal_path_length(self.ki, det.position)
            scaled_det_r = (q_space_readouts[idet,:, 1]/ipl)
            self.detector_readout.append(scaled_det_r)
            det.spectrum = np.concatenate((self.q_range[:,None], scaled_det_r[:,None]),
                                          axis =1)
    def update(self):
        self.calc_pole()
        self.update_pole_positions()
        self.update_proj_Ks()
        self.get_detector_signal()
        self.get_inplane_pole_figure()
