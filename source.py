import numpy as np

class Source:
    def __init__(self, position):
        self.name = 'source'
        self.position = np.asarray(position)
    def __repr__(self):
        return self.name