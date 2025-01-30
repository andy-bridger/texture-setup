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
        self.scheme = scheme.lower()
        self.phi = phi
        self.theta = theta
        self.psi = psi
        self.init_rot_ax()
        if self.scheme == 'euler':
            self.update_norms = self.update_norms_euler
            self.update_rot = self.update_rot_euler
            self.labels = [r'$\phi_{1}$', r'$\Phi$', r'$\phi_{2}$']
        else:
            self.update_norms = self.update_norms_rot1
            self.update_rot = self.update_rot_rot1
            self.labels = [r'$\chi$', r'$\omega$', r'$\phi$']
        self.update_rot()
        self.scale = scale
        self.update_equators()
        self.update_fracs()
        self.exp_runs= exp_runs
        self.translation = offset

    def equator(self, r = 1.0, res = 100, offset=(0,0,0)):
        return equator(r*self.scale, res, offset)

    def init_rot_ax(self):
        if self.scheme == 'euler':
            self.ax1_init = np.array((0, 0, 1))
            self.ax2_init = np.array((1, 0, 0))
            self.ax3_init = np.array((0, 0, 1))
            self.ax1 = np.array((0, 0, 1))
            self.ax2 = np.array((1, 0, 0))
            self.ax3 = np.array((0, 0, 1))
        else:
            self.ax1_init = np.array((0, 1, 0))
            self.ax2_init = np.array((0, 0, 1))
            self.ax3_init = np.array((0, 0, 1))
            self.ax1 = np.array((0, 1, 0))
            self.ax2 = np.array((0, 0, 1))
            self.ax3 = np.array((0, 0, 1))

    def update_norms_euler(self):
        rot1 = Rotation.from_euler('Z', self.phi, degrees=True)
        rot2 = Rotation.from_euler('ZX', [self.phi, self.theta], degrees=True)
        self.ax2 = rot1.apply(self.ax2_init)
        self.ax3 = rot2.apply(self.ax3_init)

    def update_norms_rot1(self):
        rot1 = Rotation.from_euler('yz', [self.phi, self.theta], degrees=True)
        self.ax3 = rot1.apply(self.ax3_init)

    def update_equators(self):
        self.ax1_eq = orient_to_pole(self.equator(r = 3), self.ax1)
        self.ax2_eq = orient_to_pole(self.equator(r = 2.9), self.ax2)
        self.ax3_eq = orient_to_pole(self.equator(r = 2.8), self.ax3)

    def update_fracs(self):
        self.phi_frac = int((self.phi/360)*self.ax1_eq.shape[1])
        self.theta_frac = int((self.theta/360)*self.ax2_eq.shape[1])
        self.psi_frac = int((self.psi/360)*self.ax3_eq.shape[1])


    def update_rot_euler(self):
        self.rot = Rotation.from_euler('ZXZ', (self.phi, self.theta, self.psi), degrees=True)
        self.rot_mat = self.rot.as_matrix()

    def update_rot_rot1(self):
        init_rot = Rotation.from_euler('yz', (self.phi, self.theta), degrees=True)
        final_rotvec = init_rot.apply(np.array((0,0,1)))*(self.psi*np.pi/180)
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
    def update(self):
        self.update_rot()
        self.update_norms()
        self.update_equators()
        self.update_fracs()

