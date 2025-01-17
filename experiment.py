import os
import numpy as np
from detector import Detector

class ExperimentalData:
    def __init__(self, detectors, name = 'test', from_data = False):
        if from_data == False:
            self.goniometer_positions = []
            self.detectors = detectors
            self.detector_readouts = []
            self.name = name
        else:
            pass

    def save_exp(self):
        np.save(f"./experiments/{self.name}_goniometer_angles", np.asarray(self.goniometer_positions))
        print("Readout Shape: ",  np.asarray(self.detector_readouts).shape)
        np.save(f"./experiments/{self.name}_detector_readouts", np.asarray(self.detector_readouts))
        np.save(f"./experiments/{self.name}_detector_positions",
                np.asarray([det.position for det in self.detectors]))
    def from_data(self, name, dir = "./experiments"):
        dets = [Detector(pos, 'default_detector',) for pos in np.load(f"{dir}/{name}_detector_positions.npy")]
        self.detectors = dets
        self.detector_readouts = np.load(f"{dir}/{name}_detector_readouts.npy")
        self.goniometer_positions = np.load(f"{dir}/{name}_goniometer_angles.npy")
