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

def get_sphere_cart_array(r = 1, res = 100, offset=(0,0,0)):
    u = np.linspace(0, 2 * np.pi, res)
    v = np.linspace(0, np.pi, res)
    x = r*np.outer(np.cos(u), np.sin(v)) + offset[0]
    y = r*np.outer(np.sin(u), np.sin(v))+ offset[1]
    z = r*np.outer(np.ones(np.size(u)), np.cos(v))+ offset[2]
    return np.concatenate((x[None,:],y[None,:],z[None,:]), axis = 0)

def cuboid_data(center, size, rotation_matrix):
    """
    Generate the X, Y, Z coordinates for a cuboid.

    Parameters:
        center: tuple of 3 floats
            Coordinates of the cuboid's center (Px, Py, Pz).
        size: tuple of 3 floats
            Dimensions of the cuboid (a, b, c).
        rotation_matrix: 3x3 numpy array
            Rotation matrix to orient the cuboid.

    Returns:
        X, Y, Z: 2D numpy arrays
            Coordinates for the cuboid's surface.
    """
    # Extract cuboid size
    a, b, c = size
    # Generate corner points of the cuboid in local coordinates
    x = np.array([-0.5, 0.5]) * a
    y = np.array([-0.5, 0.5]) * b
    z = np.array([-0.5, 0.5]) * c
    # Create a meshgrid for the cuboid
    x, y, z = np.meshgrid(x, y, z, indexing="ij")

    # Flatten and combine into an array of points
    points = np.array([x.flatten(), y.flatten(), z.flatten()]).T

    # Apply the rotation
    rotated_points = points @ rotation_matrix.T

    # Translate to the center
    translated_points = rotated_points + np.array(center)

    # Reshape into 2D arrays for plotting
    X = translated_points[:, 0].reshape(2, 2, 2)
    Y = translated_points[:, 1].reshape(2, 2, 2)
    Z = translated_points[:, 2].reshape(2, 2, 2)

    return X, Y, Z

def equator(r = 1, res = 100, offset=(0,0,0)):
    u = np.linspace(0, 2 * np.pi, res)
    x = r*np.cos(u) + offset[0]
    y = r*np.sin(u) + offset[1]
    z = np.zeros_like(x) + offset[2]
    return np.concatenate((x[None,:],y[None,:],z[None,:]), axis = 0)

def get_rot_to_orient_pole_to_z(view_axis):
    bisect = view_axis + ((np.array((0,0,1)) - view_axis)/2)
    bisect = (bisect/np.linalg.norm(bisect))*-np.pi
    return Rotation.from_rotvec(bisect)

def orient_to_pole(eq, view_axis):
    rot = get_rot_to_orient_pole_to_z(view_axis)
    return rot.apply(eq.T).T