import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
from goniometer import Goniometer
from detector import Detector
from source import Source
from sample import Sample
from helper_funcs import *
from experiment import ExperimentalData


class Mantex():
    def __init__(self, source, detectors, sample, goniometer, ewald_radii, pole_radius,
                 sample_view_axis, detector_colors = None, ewald_steps = 2, q_probe = 1):
        self.source = source
        self.detectors = detectors
        self.sample = sample
        self.goniometer = goniometer
        self.sample_view_axis = np.asarray(sample_view_axis)
        self.no_pfi = True
        self.pole_radius = pole_radius
        self.ki_raw = self.sample.position - self.source.position
        self.ki_raw_scale = np.linalg.norm(self.ki_raw)
        self.ki = self.ki_raw/self.ki_raw_scale
        self.Ks = self.get_Ks()


    def get_Ks(self):
        norm = lambda x : x/np.linalg.norm(x)
        ki = norm(self.sample.position - self.source.position)
        Ks = []
        for det in self.detectors:
            K = norm(norm(det.position - self.sample.position) - ki)
            Ks.append(K)
            det.K = K
        return Ks

    def calc_pole(self):
        self.eq_cart = orient_to_pole(equator(r=self.pole_radius), self.pole_view_axis)
        self.pole_Ks = self.get_pole_K_intercepts()
        self.pole_K_projs = np.asarray([self.get_orientated_stereographic_projection(pK) for pK in self.pole_Ks])
        self.pole_K_projs_vec = np.asarray([self.pole_K_projs[i] - K for i, K in enumerate(self.pole_Ks)])


    def get_orientated_stereographic_projection(self, point):
        spoint = orient_to_pole(point, self.north_pole)
        oppole = np.array((0,0,-np.sign(spoint[2])*self.pole_radius))
        svec = spoint - oppole
        u = -oppole[2]/svec[2]
        proj_point = (u*svec) + oppole
        return orient_to_pole(proj_point, self.north_pole)


    def get_pole_K_intercepts(self):
        Gs = []
        for K in self.Ks:
            roots = sphere_line_intercept(np.zeros(3), K, self.sample.position, self.pole_radius)
            Gs.append(roots[roots!=0]*K)
        return np.asarray(Gs)

    def get_pole_figure_intensities(self, readouts):
        self.pole_view_axis = self.sample.get_view_in_lab_frame(self.sample_view_axis)
        self.north_pole = self.pole_view_axis / np.linalg.norm(self.pole_view_axis)
        self.south_pole = -self.north_pole

        self.calc_pole()
        self.pole_figure_points = self.sample.rot.apply(self.pole_K_projs, inverse=True)
        pfi = []
        for di in range(len(self.detectors)):#just read the probe channel for now
            dr = readouts[di]
            npfi = np.concatenate((self.pole_figure_points[di], np.array((dr,))))
            pfi.append(npfi)
        pfi = np.asarray(pfi)
        return pfi[np.argsort(pfi[:,3])]

