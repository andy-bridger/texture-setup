import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
from goniometer import Goniometer
from detector import Detector
from source import Source
from sample import Sample
from helper_funcs import *
from experiment import ExperimentalData
from model import Mantex

class MantexSim(Mantex):
    def __init__(self, source, detectors, sample, goniometer, ewald_radii,
                 sample_view_axis, detector_colors = None, ewald_steps = 2, q_probe = 1):

        super().__init__(source, detectors, sample, goniometer, sample_view_axis)

        self.exp_data = []
        self.r_dict = {'recip_sample': 1,
                                   'lab_sample': 1,
                                   'k_vecs': 1,
                                   'gonio_r': 1,
                                   'gonio_v': 1,
                                    'q_probe': q_probe}
        self.ewald_radii = ewald_radii
        self.ewald_radius = ewald_radii[0]
        self.q_range = np.linspace(ewald_radii[0], ewald_radii[1], ewald_steps)
        self.pole_figure_intensities = np.array(())
        if detector_colors == None:
            self.detector_colors = ['black']*len(detectors)
        else:
            self.detector_colors = detector_colors
        for det in self.detectors:
            det.ki = self.ki
            det.radii = ewald_radii
            det.r_steps = ewald_steps

        self.Ks = self.get_Ks()
        self.get_detector_Ks()
        self.setup_qrange_rl()
        self.setup_pole()

        self.get_probe_ind()

        self.show_lab_K_vecs = True
        self.lab_K_vecs = []
        self.get_inplane_pole_figure()

        self.update()
        self.create_experiment()

    def setup_pole(self):
        self.pole_cart_array = get_sphere_cart_array(r=1, offset=(0, 0, 0))

    def get_probe_ind(self):
        self.probe_ind = np.searchsorted(self.q_range, self.r_dict['q_probe'], 'left')
        for det in self.detectors:
            det.readout_probe = self.probe_ind

    def update(self):
        self.update_pole_positions()
        self.update_proj_Ks()
        self.get_detector_signal()
        self.get_probe_ind()
        self.get_inplane_pole_figure()

    def setup_qrange_rl(self):
        self.qs_of_interest = []
        for iq, q_probe in enumerate(self.q_range):
            qs_of_interest = []
            for det in self.detectors:
                qs_of_interest.append(self.sample.get_q_shell(np.linalg.norm(det.pKs[iq]),0.2))
            self.qs_of_interest.append(np.asarray(qs_of_interest).sum(axis = 0)>0)

    def get_detector_Ks(self):
        pKs = []
        for det in self.detectors:
            det_pKs = det.get_polychromatic_ks(self.q_range, self.ki)

            pKs.append(det_pKs)
        self.pKs = np.asarray(pKs)

    def update_proj_Ks(self):
        self.pole_K_projs_vec = np.asarray([self.pole_K_projs[i] - K for i, K in enumerate(self.pole_Ks)])


    def get_detector_signal(self, thresh = 0.1):
        q_space_readouts = []
        for iq, q_probe in enumerate(self.q_range):

            q_of_interest = self.qs_of_interest[iq]
            # detector, grain, coord, reflection  det qrange coord
            pKs = self.pKs[:,iq, :]
            disp_vecs = self.sample.lab_space_rlatts[None, :,:,q_of_interest] - pKs[:,None,:,None]
            dists = np.linalg.norm(disp_vecs, axis = 2)
            relevant_dists = np.where(dists < thresh, np.exp(-500*dists**2), 0)
            q_space_readouts.append(relevant_dists.sum(axis = (1,2)))
        q_space_readouts = np.asarray(q_space_readouts)

        # placeholder to keep other functions working
        self.detector_readout = q_space_readouts.T
        for idet, det in enumerate(self.detectors):
            det.readout = self.detector_readout[idet]

    def get_inplane_pole_figure(self):
        self.get_detector_signal()

        readouts = [det.get_readout() for det in self.detectors]
        if len(self.pole_figure_intensities) == 0:
            self.pole_figure_intensities = self.get_pole_figure_intensities(readouts)
        else:
            new_pfi = self.get_pole_figure_intensities(readouts)
            insertion_inds = np.searchsorted(self.pole_figure_intensities[:,-1], new_pfi[:,-1])
            self.pole_figure_intensities = np.insert(self.pole_figure_intensities, insertion_inds, new_pfi, axis = 0)


    def calc_ewald(self, er):
        return get_sphere_cart_array(r = er, offset = -self.ki*er)


    def get_cube(self, size = (0.1, 0.1, 0.1),
                  offset = (0,0,0), use_gonio = True):
        if use_gonio == True:
            cube_data = cuboid_data(offset, size, self.sample.lab_frame_axes)
        else:
            cube_data = cuboid_data(offset, size, np.eye(3,3))
        return cube_data

    def get_lab_frame(self):
        source_repr = self.get_cube(offset = self.source.position, size=[2,1,1], use_gonio= False)
        source_col = ['grey']*3
        det_reprs = []
        det_cols = []
        for idet, det in enumerate(self.detectors):
            det_reprs.append(self.get_cube(offset=det.position, size=[1, 1, 2], use_gonio=False))
            det_cols.append([self.detector_colors[idet]] * 3)
        return source_repr, source_col, det_reprs, det_cols

    def get_sample(self, ratio = 0.1):
        cube = self.get_cube(offset = (0,0, 0), size = self.ewald_radius*np.ones(3)*ratio)
        return cube

    def create_experiment(self):
        self.exp_data = ExperimentalData(self.detectors)

    def get_detector_vectors(self):
        return [det.position - self.sample.position for det in self.detectors]