import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation

from helper_funcs import sphere_line_intercept


class Detector:
    def __init__(self, position, name):
        self.name = f'detector {name}'
        self.position = np.asarray(position)
        self.lab_frame_axes = np.array(((1,0,0),(0,1,0),(0,0,1)))
        self.K = None
        self.ki = None
        self.readout = np.zeros(2)
        self.readout_probe = 0
    def __repr__(self):
        return self.name
    def get_polychromatic_ks(self, rs, ki):
        pKs = []
        for r in rs:
            roots = sphere_line_intercept(np.zeros(3), self.K, -ki*r, r)
            pKs.append(roots[roots != 0] * self.K)
        self.pKs = np.asarray(pKs)
        return self.pKs
    def get_readout(self):
        return np.log(self.readout[self.readout_probe] + 1)