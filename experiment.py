import os
import numpy as np

class ExperimentalData:
    def __init__(self, detectors, name = 'test'):
        self.goniometer_positions = []
        self.detectors = detectors
        self.detector_readouts = []
        self.name = name

    def save_exp(self):
        np.save(f"./experiments/{self.name}_goniometer_angles", np.asarray(self.goniometer_positions))
        np.save(f"./experiments/{self.name}_detector_readouts", np.asarray(self.detector_readouts))
        np.save(f"./experiments/{self.name}_detector_positions",
                np.asarray([det.position for det in self.detectors]))