import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation

from helper_funcs import equator


class Goniometer:
    def __init__(self, scale = 1, exp_runs = ()):
        self.z_norm_init = np.array((0,0,1))
        self.x_prime_norm_init = np.array((1,0,0))
        self.z_prime_norm_init = np.array((0,0,1))
        self.z_norm = np.array((0,0,1))
        self.x_prime_norm = np.array((1,0,0))
        self.z_prime_norm = np.array((0,0,1))
        self.scale = scale
        self.update_equators()
        self.phi_frac = 0
        self.theta_frac = 0
        self.psi_frac = 0
        self.exp_runs= exp_runs

    def equator(self, r = 1.0, res = 100, offset=(0,0,0)):
        return equator(r*self.scale, res, offset)

    def orient_to_pole(self, eq, view_axis):
        bisect = view_axis + ((np.array((0,0,1)) - view_axis)/2)
        bisect = (bisect/np.linalg.norm(bisect))*-np.pi
        rot = Rotation.from_rotvec(bisect)
        return rot.apply(eq.T).T

    def update_norms(self, phi, theta):
        rot1 = Rotation.from_euler('Z', phi, degrees=True)
        rot2 = Rotation.from_euler('ZX', [phi, theta], degrees=True)
        self.x_prime_norm = rot1.apply(self.x_prime_norm_init)
        self.z_prime_norm = rot2.apply(self.z_prime_norm_init)

    def update_equators(self):
        self.z_eq = self.orient_to_pole(self.equator(r = 3), self.z_norm)
        self.x_prime_eq = self.orient_to_pole(self.equator(r = 2.95), self.x_prime_norm)
        self.z_prime_eq = self.orient_to_pole(self.equator(r = 2.9), self.z_prime_norm)

    def update_fracs(self, phi, theta, psi):
        self.phi_frac = int((phi/360)*self.z_eq.shape[1])
        self.theta_frac = int((theta/360)*self.x_prime_eq.shape[1])
        self.psi_frac = int((psi/360)*self.z_prime_eq.shape[1])