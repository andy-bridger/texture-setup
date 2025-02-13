import sys
import os

import numpy as np

# Add the parent directory to the system path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from presenter import Presenter
from helper_funcs import *

class NyrtexPresenter(Presenter):
    def __init__(self, Model, View):
        self.Model = Model
        self.View = View
        self.init_view()

    def init_view(self):
        self.View.setup_view()
        self.View.update_view_axes()
        #self.View.add_lab_k_widgets(self.lab_k_update)
        #self.View.add_sample_scale_widget(self.sample_scale_update, self.Model.r_dict['lab_sample'])
        #self.View.add_gonio_ring_scale_widget(self.gr_scale_update)
        #self.View.add_gonio_axis_scale_widget(self.ga_scale_update)
        self.View.add_probe_pos_widget(self.update_probe_pos,
                                       self.Model.q_range[0],
                                       self.Model.q_range[-1],
                                       self.Model.q_probe)
        self.View.add_run_num_widget(self.update_run, 0, self.Model.raw_exp_data.detector_readouts.shape[0], 0)

    def plot_readouts(self):
        self.det_readout_plot_artist = self.View.show_detector_readout(self.View.det_ax,
                                                                       self.Model.get_current_readouts(),
                                                                       self.Model.detector_colors)

    def plot_lab_sample(self, ax):
        self.lab_sample_artist = self.View.plot_mesh(ax, self.Model.mantex.sample.get_mesh_vectors(self.Model.r_dict['lab_sample']))

    def plot_ref_sample(self, ax):
        self.lab_sample_artist = self.View.plot_mesh(ax, self.Model.reference_sample.get_mesh_vectors(self.Model.r_dict['lab_sample']))

    def plot_lab_beam_paths(self):
        self.det_path_artists = self.View.plot_lab_detector_paths(self.View.lab_ax, self.Model.sample.position,
                                                                  self.Model.get_detector_vectors())
        self.source_path_artists = self.View.plot_lab_source_path(self.View.lab_ax, self.Model.source.position,
                                                                  self.Model.mantex.ki_raw)

    def plot_sample_beam_paths(self):
        self.det_path_artists = self.View.plot_lab_detector_paths(self.View.sample_ax, self.Model.sample.position,
                                                                  self.Model.get_norm_detector_vectors())
        self.source_path_artists = self.View.plot_lab_source_path(self.View.sample_ax, -self.Model.mantex.ki,
                                                                  self.Model.mantex.ki)

    def plot_lab_Ks(self):
        self.lab_K_vec_artists = self.View.plot_lab_K_vecs(self.View.lab_ax,
                                                           self.Model.sample.position,
                                                           self.Model.mantex.Ks,
                                                           self.Model.detector_colors,
                                                           scale = self.Model.mantex.ki_raw_scale*self.Model.r_dict['k_vecs']/4)

    def plot_sample_lab_Ks(self):
        self.lab_K_vec_artists = self.View.plot_lab_K_vecs(self.View.sample_ax,
                                                           self.Model.sample.position,
                                                           self.Model.mantex.Ks,
                                                           self.Model.detector_colors,
                                                           scale = self.Model.r_dict['k_vecs'])

    def plot_pole_figure_detectors(self):
        self.pole_figure_detector_artists = self.View.plot_pole_figure_position(self.View.pole_proj_ax,
                                                                                self.Model.mantex.pole_figure_points,
                                                                                self.Model.detector_colors,
                                                                                equator())

    def plot_goniometer(self, ax):
        self.Model.mantex.smp.update()
        self.z_artist = [self.View.plot_line(ax, self.Model.mantex.smp.ax1_eq, self.Model.r_dict['gonio_r']),
                         self.View.plot_line(ax, self.Model.mantex.smp.ax1_eq[:,:self.Model.mantex.smp.phi_frac],
                                             self.Model.r_dict['gonio_r'], 'red')]
        self.xp_artist = [self.View.plot_line(ax, self.Model.mantex.smp.ax2_eq, self.Model.r_dict['gonio_r']),
                          self.View.plot_line(ax,
                                              self.Model.mantex.smp.ax2_eq[:, :self.Model.mantex.smp.theta_frac],
                                              self.Model.r_dict['gonio_r'], 'blue')]
        self.zp_artist = [self.View.plot_line(ax, self.Model.mantex.smp.ax3_eq, self.Model.r_dict['gonio_r']),
                         self.View.plot_line(ax, self.Model.mantex.smp.ax3_eq[:, :self.Model.mantex.smp.psi_frac],
                                self.Model.r_dict['gonio_r'], 'green')]
    def plot_goniometer_axes(self, ax):
        self.View.plot_goniometer_axis(ax, 'Z',
                                       self.Model.mantex.sample.position,
                                       self.Model.mantex.smp.ax1,
                                       self.Model.r_dict['gonio_v'])
        self.View.plot_goniometer_axis(ax, "X'",
                                       self.Model.mantex.sample.position,
                                       self.Model.mantex.smp.ax2,
                                       self.Model.r_dict['gonio_v'])
        self.View.plot_goniometer_axis(ax, "Z'",
                                       self.Model.mantex.sample.position,
                                       self.Model.mantex.smp.ax3,
                                       self.Model.r_dict['gonio_v'])

    def plot_readouts(self):
        self.det_readout_plot_artist = self.View.show_detector_readout(self.View.det_ax,
                                                                       self.Model.get_current_readouts(),
                                                                       self.Model.detector_colors,)

    def plot_q_probe(self):
        self.det_probe_artist = self.View.show_detector_probe(self.View.det_ax,
                                                              self.Model.Alg.q_probe,
                                                              self.Model.probe_window)
    def add_ax_to_plot(self, ax, mat, pos, length = 1, dim='3d'):
        if dim == '2d':
            return self.View.plot_cart_ax_2d(ax, mat, pos, length)
        else:
            return self.View.plot_cart_ax(ax, mat, pos, length)

    def plot_all(self):
        # Lab Frame
        self.plot_lab_sample(self.View.lab_ax)
        self.plot_lab_beam_paths()
        self.plot_lab_Ks()
        self.plot_lab_components()
        self.plot_goniometer(self.View.lab_ax)
        #self.add_ax_to_plot(self.View.lab_ax, (3, 3, 3))

        # Sample Frame
        self.plot_lab_sample(self.View.sample_ax)
        self.plot_goniometer(self.View.sample_ax)
        self.plot_sample_beam_paths()
        self.plot_sample_lab_Ks()
        s = self.Model.get_mesh_size()
        self.add_ax_to_plot(self.View.sample_ax, self.Model.mantex.sample.smp.get_rot(), np.array((3,0,0))*s, 1)

        # Object Reference
        self.plot_ref_sample(self.View.object_ax)
        s = self.Model.get_mesh_size()
        self.add_ax_to_plot(self.View.object_ax, np.eye(3,3), np.array((3,0,0))*s, s)
        self.add_angles_titles()

        # Detector
        self.plot_readouts()

        # Pole Detector Frame
        self.plot_pole_figure_detectors()
        self.plot_q_probe()
        self.add_ax_to_plot(self.View.pole_proj_ax, np.eye(2,2), (-1,-1), dim ='2d', length=0.2)

        # Pole Figure Frame
        self.plot_pole_figure_intensities()
        self.add_ax_to_plot(self.View.calc_pf_ax, np.eye(2,2),(-1, -1), dim='2d', length=0.2)

        self.View.fix_aspect()

    def add_angles_titles(self):
        smp = self.Model.mantex.sample.smp
        self.View.sample_ax.set_title(f"{smp.labels[0]}: {smp.phi}, {smp.labels[1]}: {smp.theta}, {smp.labels[2]}: {smp.psi}")

    def sample_scale_update(self, val):
        self.remove_artist_set(self.lab_sample_artist)
        self.Model.r_dict['lab_sample'] = val
        #if not self.Model.sample.primitive:
        #    self.Model.sample.scale = val
        self.plot_lab_sample(self.View.lab_ax)
        self.plot_lab_sample(self.View.sample_ax)
        self.View.fix_aspect()

    def plot_pole_figure_intensities(self):
        self.pole_figure_artists = []
        self.pole_figure_artists += self.View.plot_pole_figure_intensities(self.View.calc_pf_ax,
                                                                                self.Model.Alg.pole_figure_intensities,self.Model.mantex.pole_figure_points,
                                                                                self.Model.detector_colors, equator())

    def update_probe_pos(self, val):
        self.remove_artist_set(self.det_probe_artist)
        self.remove_artist_set(self.pole_figure_artists)
        self.Model.Alg.q_probe = val
        self.Model.Alg.execute()
        self.plot_pole_figure_intensities()
        self.plot_q_probe()

    def update_run(self, val):
        self.Model.run_index = val
        self.Model.update_mantex()
        self.View.lab_ax.clear()
        self.View.sample_ax.clear()
        self.View.pole_proj_ax.clear()
        self.View.det_ax.clear()
        self.View.calc_pf_ax.clear()
        self.View.update_view_axes()
        self.plot_all()
