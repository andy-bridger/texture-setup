import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
from goniometer import Goniometer
from detector import Detector
from source import Source
from helper_funcs import *
from experiment import ExperimentalData
from scipy.stats import norm
from scipy.optimize import curve_fit


class Mantex():
    def __init__(self, source, detectors, sample, smp,
                 sample_view_axes, q_probe, probe_window = 0.05):
        self.source = source
        self.detectors = detectors
        self.sample_position = sample.position
        self.sample = sample
        self.smp = smp
        self.sample_view_axes = np.asarray(sample_view_axes)
        self.sample_view_axis = np.cross(sample_view_axes[0], sample_view_axes[1])
        self.sample_view_mat = np.concatenate((self.sample_view_axes, self.sample_view_axis[None, :]), axis = 0).T
        self.ki_raw = self.sample_position - self.source.position
        self.ki_raw_scale = np.linalg.norm(self.ki_raw)
        self.ki = self.ki_raw/self.ki_raw_scale
        self.Ks = self.get_Ks()
        self.q_probe = q_probe
        self.probe_window = probe_window


    def get_Ks(self):
        norm = lambda x : x/np.linalg.norm(x)
        ki = norm(self.sample_position - self.source.position)
        Ks = []
        for det in self.detectors:
            K = norm(norm(det.position - self.sample_position) - ki)
            Ks.append(K)
            det.K = K
        return Ks

    def calc_pole(self):
        self.pole_Ks = self.get_pole_K_intercepts()
        self.update_pole_positions()

    def get_pole_K_intercepts(self):
        Gs = []
        for K in self.Ks:
            roots = sphere_line_intercept(np.zeros(3), K, self.sample_position, 1)
            Gs.append(roots[roots!=0]*K)
        return np.asarray(Gs)

    def update_pole_positions(self):
        # for calculating the pole figure and doing the stereographic projection, we want to map the desired inplane view axes
        # U and V, onto (1,0,0) and (0,1,0). This transformation can then be applied to all the K vectors and the projection
        # becomes straightforward as the north/south pole is (0,0 +/-1). The in plane projection can just be read off when the line
        # between the reprojected K and the pole crosses z = 0.

        # To find the mapping we need to know what U and V are, within the lab frame (this way we can use our instrument coordinates straight away)
        # so we apply the State Matrix Providers rotation to our sample view matrix (vectors U, V and W) to get U', V', and W' I guess.
        # to find the transformation of these onto (1,0,0) and (0,1,0) is then conveniently just the inverse of this matrix.
        # This is because matrix multiplication is really just mapping points to other points and the inverse maps to the identity matrix.


        self.pole_view_mat = (self.smp.get_rot()@(self.sample_view_mat))
        self.to_pole_view = np.linalg.inv(self.pole_view_mat)


        # Project the K vectors through the equator of the new orientation
        pK_projs = []
        pK_sters = []
        for pK in self.pole_Ks:
            pK_proj, pK_ster =  self.get_orientated_stereographic_projection(pK)
            pK_projs.append(pK_proj)
            pK_sters.append(pK_ster)
        self.pole_K_projs = np.asarray(pK_projs)
        self.pole_figure_points = np.asarray(pK_sters)

    def get_orientated_stereographic_projection(self, point):
        # as mentioned before, we want to use our to_pole_view to reproject the points
        spoint = (self.to_pole_view@point)
        # we want the vector to whichever pole is on the opposite hemisphere
        oppo_pole = np.array((0,0,-np.sign(spoint[2])))
        svec = spoint - oppo_pole

        if spoint[2] == 0: # point is within the projection plane
            return point, spoint
        else:
            u = -oppo_pole[2]/svec[2]
            proj_point = (u*svec) + oppo_pole
            return (self.pole_view_mat@proj_point), proj_point


    def get_pole_figure_intensities(self, spectra):
        self.calc_pole() # calculate the new projections of the detector Ks

        # append the readout intensity to the sample space position information
        pfi = []
        for di, det in enumerate(self.detectors):
            dr = self.readout_spectrum(spectra[di], det)
            npfi = np.concatenate((self.pole_figure_points[di], np.array((dr,))))
            pfi.append(npfi)
        pfi = np.asarray(pfi)
        return pfi[np.argsort(pfi[:,3])]

    def readout_spectrum(self, spectrum, det):
        return self.BadFitSpectrum(self.correct_spectrum_for_attenuation(spectrum, det))

    def correct_spectrum_for_attenuation(self, spectrum, det):
        ipl, entry_point, exit_point = self.sample.get_internal_path_length(self.ki, det.position)
        spectrum[:,1] *= ipl
        return spectrum

    def FitSpectrum(self, spectrum):
        # can imagine how in mantid implementation this will be a bit jazzier
        popt, pcov = curve_fit(gaussian, spectrum[:,0], spectrum[:,1],
                               bounds = ((0, self.q_probe-self.probe_window, 0),
                                         (1e6, self.q_probe+self.probe_window, 2)))
        return popt[0] # return the amplitude of the fit gaussian
    def BadFitSpectrum(self, spectrum):
        lb = np.searchsorted(spectrum[:,0], self.q_probe-self.probe_window)
        ub = np.searchsorted(spectrum[:, 0], self.q_probe + self.probe_window)
        spec_win = spectrum[lb:ub, 1]
        return spec_win.max()-spec_win.min() #Gaussian fit is kind of rubbish anyway
                                             #might as well have a simple fit that can be interpreted
                                             #on visual inspection

def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-((x - mean) / 4 / stddev)**2)
