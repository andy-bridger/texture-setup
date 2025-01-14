import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation

class Detector:
    def __init__(self, position, name):
        self.name = f'detector {name}'
        self.position = np.asarray(position)
        self.lab_frame_axes = np.array(((1,0,0),(0,1,0),(0,0,1)))
    def __repr__(self):
        return self.name