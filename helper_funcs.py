import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation

def sphere_line_intercept(p1, p2, p3, r):
    """
    p1: one point on line (ideally vector representation of line so this can be (0,0,0)
    p2: another point on line (if vector repr, p2 = vec)
    p3: centre of sphere
    r:  radius of sphere
    """

    a = np.sum((p2 - p1) ** 2)
    b = 2 * (np.dot(p2 - p1, p1 - p3))
    c = np.sum(p3 ** 2) + np.sum(p1 ** 2) - 2 * np.dot(p3, p1) - r ** 2

    roots = np.array((-b + np.sqrt((b ** 2) - 4 * a * c)) / (2 * a), (-b - np.sqrt((b ** 2) - 4 * a * c)) / (2 * a))

    return roots