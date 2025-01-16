import os
import numpy as np

class ExperimentalData:
    def __init__(self, detectors):
        self.goniometer_positions = []
        self.detectors = detectors
        self.detector_readouts = []

    def save_exp(self, exp_name='test'):
        np.save(f"./experiments/{exp_name}_goniometer_angles", np.asarray(self.goniometer_positions))
        np.save(f"./experiments/{exp_name}_detector_readouts", np.asarray(self.detector_readouts))
        np.save(f"./experiments/{exp_name}_detector_positions",
                np.asarray([det.position for det in self.detectors]))