from xml.sax.saxutils import prepare_input_source

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation

from helper_funcs import *


class GenericStateMatrixProvider:
    def __init__(self, rot_mat, translation):
        self.rot_mat = rot_mat
        self.translation = translation
    def get_rot(self):
        return self.rot_mat
    def get_translation(self):
        return self.translation
    def to_state_data(self):
        return np.concatenate((self.rot_mat.reshape((-1)), self.translation))
    def from_state_data(self, state_data):
        self.translation = state_data[-3:]
        self.rot_mat = state_data[:-3].reshape((3,3))


class Goniometer(GenericStateMatrixProvider):
    def __init__(self, scale = 1, phi =0,theta=0, psi=0, exp_runs = (), offset = np.zeros(3), scheme = 'euler'):
        super().__init__(None, None)
        self.scheme = scheme
        self.phi = phi
        self.theta = theta
        self.psi = psi
        self.update_rot()
        self.z_norm_init = np.array((0,0,1))
        self.x_prime_norm_init = np.array((1,0,0))
        self.z_prime_norm_init = np.array((0,0,1))
        self.z_norm = np.array((0,0,1))
        self.x_prime_norm = np.array((1,0,0))
        self.z_prime_norm = np.array((0,0,1))
        self.scale = scale
        self.update_equators()
        self.update_fracs()
        self.exp_runs= exp_runs
        self.translation = offset

    def equator(self, r = 1.0, res = 100, offset=(0,0,0)):
        return equator(r*self.scale, res, offset)

    def update_norms(self):
        rot1 = Rotation.from_euler('Z', self.phi, degrees=True)
        rot2 = Rotation.from_euler('ZX', [self.phi, self.theta], degrees=True)
        self.x_prime_norm = rot1.apply(self.x_prime_norm_init)
        self.z_prime_norm = rot2.apply(self.z_prime_norm_init)

    def update_equators(self):
        self.z_eq = orient_to_pole(self.equator(r = 3), self.z_norm)
        self.x_prime_eq = orient_to_pole(self.equator(r = 2.95), self.x_prime_norm)
        self.z_prime_eq = orient_to_pole(self.equator(r = 2.9), self.z_prime_norm)

    def update_fracs(self):
        self.phi_frac = int((self.phi/360)*self.z_eq.shape[1])
        self.theta_frac = int((self.theta/360)*self.x_prime_eq.shape[1])
        self.psi_frac = int((self.psi/360)*self.z_prime_eq.shape[1])

    def update_rot(self):
        if self.scheme == 'euler':
            self.rot = Rotation.from_euler('ZXZ', (self.phi, self.theta, self.psi), degrees=True)
            self.rot_mat = self.rot.as_matrix()
        else:
            init_rot = Rotation.from_euler('yz', (self.phi, self.theta), degrees=True)
            final_rotvec = init_rot.apply(np.array((0,0,1)))*(self.psi/360)
            second_rot = Rotation.from_rotvec(final_rotvec, degrees=False)
            self.rot_mat = second_rot.as_matrix()@init_rot.as_matrix()
            self.rot = Rotation.from_matrix(self.rot_mat)

    def update_angles(self, phi, theta, psi):
        self.phi = phi
        self.theta = theta
        self.psi = psi
        self.update_rot()
    def get_rot(self):
        return self.rot_mat
    def get_translation(self):
        return self.translation

