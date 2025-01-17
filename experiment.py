import os
import numpy as np
from detector import Detector
from goniometer import GenericStateMatrixProvider

class ExperimentalData:
    def __init__(self, detectors, name = 'test', from_data = False):
        if from_data == False:
            self.state_data = []
            self.detectors = detectors
            self.detector_readouts = []
            self.name = name
            self.smps = ()
        else:
            pass

    def save_exp(self):
        np.save(f"./experiments/{self.name}_state_data", np.asarray(self.state_data))
        print("Readout Shape: ",  np.asarray(self.detector_readouts).shape)
        np.save(f"./experiments/{self.name}_detector_readouts", np.asarray(self.detector_readouts))
        np.save(f"./experiments/{self.name}_detector_positions",
                np.asarray([det.position for det in self.detectors]))
    def from_data(self, name, dir = "./experiments"):
        dets = [Detector(pos, 'default_detector',) for pos in np.load(f"{dir}/{name}_detector_positions.npy")]
        self.detectors = dets
        self.detector_readouts = np.load(f"{dir}/{name}_detector_readouts.npy")
        self.state_data = np.load(f"{dir}/{name}_state_data.npy")
        self.smps = [GenericStateMatrixProvider(None, None) for i in range(self.state_data.shape[0])]
        [self.smps[i].from_state_data(sd) for i, sd in enumerate(self.state_data)]

