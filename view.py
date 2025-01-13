import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
from model import *

from matplotlib.widgets import Slider, CheckButtons

class View():
    def __init__(self):
        pass
    def setup_view(self, pole_radius):
        self.fig = plt.figure()
        gs0 = self.fig.add_gridspec(7, 10)
        self.lab_ax = self.fig.add_subplot(gs0[0:3, :8], projection='3d')
        self.recip_ax = self.fig.add_subplot(gs0[3:6, :8], projection='3d')
        self.pole_proj_ax = self.fig.add_subplot(gs0[:3, 8:])
        self.calc_pf_ax = self.fig.add_subplot(gs0[3:6, 8:])
        self.pole_proj_ax.set_xlim([-pole_radius-0.1, pole_radius+0.1])
        self.pole_proj_ax.set_ylim([-pole_radius-0.1, pole_radius+0.1])
        self.recip_ax.set_axis_off()
        self.pole_proj_ax.set_axis_off()
        self.lab_ax.set_axis_off()
    def update_view_axes(self):
        self.recip_ax.clear()
        self.pole_proj_ax.clear()
        self.recip_ax.set_axis_off()
        self.pole_proj_ax.set_axis_off()
        self.calc_pf_ax.set_axis_off()
        self.lab_ax.set_axis_off()
    def fix_aspect(self):
        self.lab_ax.set_aspect('equal')
        self.pole_proj_ax.set_aspect('equal')
        self.recip_ax.set_aspect('equal')
        self.calc_pf_ax.set_aspect('equal')

    def show_detector_Ks(self, ax, Ks, detector_colors):
        artists = []
        for i, K in enumerate(Ks):
            artists.append(ax.quiver(0,0,0, K[0], K[1], K[2], arrow_length_ratio = 0.05, color = detector_colors[i]))
        return artists

    def show_reciprocal_lattices(self, ax, lab_space_rlatts, cell_colors, rl_alphas):
        artists = []
        for i,rl in enumerate(lab_space_rlatts):
            artists.append(ax.scatter(xs=rl[0], ys= rl[1], zs= rl[2], c=cell_colors[i], s = rl_alphas))
        return artists

    def plot_sphere(self, ax = None, cart_array = None):
        # Plot the surface
        x,y,z = cart_array
        return ax.plot_surface(x, y, z, color='grey', alpha = 0.25)

    def plot_line(self, ax = None, cart_array = None, scale = 1, color = 'grey'):
        # Plot the surface
        x,y,z = cart_array*scale
        return ax.plot(x, y, z, color=color, alpha = 1)

    def plot_cube(self, ax, cube_data, col):
        X, Y, Z = cube_data
        cube = []
        for i in range(2):  # Bottom/Top faces
            cube.append(ax.plot_surface(X[i, :, :], Y[i, :, :], Z[i, :, :], alpha=0.5, color=col[0]))
        for i in range(2):  # Front/Back faces
            cube.append(ax.plot_surface(X[:, i, :], Y[:, i, :], Z[:, i, :], alpha=0.5, color=col[1]))
        for i in range(2):  # Left/Right faces
            cube.append(ax.plot_surface(X[:, :, i], Y[:, :, i], Z[:, :, i], alpha=0.5, color=col[2]))
        return cube

    def plot_pole_figure_vectors(self, ax, pole_Ks, detector_colors, pole_K_projs_vec):
        artists = []
        for i, K in enumerate(pole_Ks):
            artists.append(ax.quiver(0,0,0, K[0], K[1], K[2], arrow_length_ratio = 0.05, color = detector_colors[i]))
            artists.append(ax.quiver(K[0], K[1], K[2],
                      pole_K_projs_vec[i][0], pole_K_projs_vec[i][1], pole_K_projs_vec[i][2],
                      arrow_length_ratio = 0.05, color = detector_colors[i]))
        return artists

    def plot_pole_figure_position(self, ax, sample_view_axis, pole_figure_points, detector_colors, eq):
        artists = []
        if np.where(sample_view_axis ==1)[0] == 0:
            artists.append(ax.scatter(pole_figure_points[:,1], pole_figure_points[:,2], c = detector_colors))
        if np.where(sample_view_axis ==1)[0] == 1:
            artists.append(ax.scatter(pole_figure_points[:,0], pole_figure_points[:,2], c = detector_colors))
        if np.where(sample_view_axis ==1)[0] == 2:
            artists.append(ax.scatter(pole_figure_points[:,0], pole_figure_points[:,1], c = detector_colors))
        artists.append(ax.plot(eq[0], eq[1], c='grey'))
        return artists

    def plot_pole_figure_intensities(self, ax, sample_view_axis, pole_figure_intensities, eq):
        artists = []
        if np.where(sample_view_axis ==1)[0] == 0:
            artists.append(ax.scatter(pole_figure_intensities[:,1], pole_figure_intensities[:,2],
                                c = pole_figure_intensities[:,3]))
        if np.where(sample_view_axis ==1)[0] == 1:
            artists.append(ax.scatter(pole_figure_intensities[:,0], pole_figure_intensities[:,2],
                                c = pole_figure_intensities[:,3]))
        if np.where(sample_view_axis ==1)[0] == 2:
            artists.append(ax.scatter(pole_figure_intensities[:,0], pole_figure_intensities[:,1],
                                c = pole_figure_intensities[:,3]))
        #eq = self.equator()
        artists.append(ax.plot(eq[0], eq[1], c = 'grey'))
        return artists

    def plot_lab_detector_paths(self,ax, sample_position, det_vec):
        artists = []
        for dv in det_vec:
            artists.append(ax.quiver(X=sample_position[0], Y=sample_position[1], Z=sample_position[2],
                               U = dv[0], V = dv[1], W = dv[2],
                               color= 'grey', arrow_length_ratio= 0.1, alpha= 0.25))
        return artists

    def plot_lab_source_path(self, ax, source_position, beam_vec):
        return ax.quiver(X = source_position[0], Y = source_position[1], Z = source_position[2],
                         U = beam_vec[0], V = beam_vec[1], W = beam_vec[2],
                         color= 'grey', arrow_length_ratio= 0.1, alpha= 0.25)

    def plot_lab_K_vecs(self, ax, sample_position, Ks, detector_colors, scale):
        artists = []
        for i, K in enumerate(Ks):
            artists.append(ax.quiver(X=sample_position[0], Y=sample_position[1], Z=sample_position[2],
                               U = K[0]*scale, V = K[1]*scale, W = K[2]*scale,
                               color= detector_colors[i], arrow_length_ratio= 0.1, alpha= 0.25))
        return artists

    def plot_goniometer_axis(self, ax, label, sample_position, vector, scale, color = 'grey'):
        artists = []
        artists.append(ax.quiver(X=sample_position[0], Y=sample_position[1], Z=sample_position[2],
                               U = vector[0]*scale, V = vector[1]*scale, W = vector[2]*scale,
                               color= color, arrow_length_ratio= 0.1, alpha= 0.25))
        return artists

    def add_goniometer_widgets(self, update_plot):
        ax_phi = plt.axes([0.2, 0.09, 0.3, 0.03], facecolor='lightgoldenrodyellow')  # Slider for φ
        ax_theta = plt.axes([0.2, 0.05, 0.3, 0.03], facecolor='lightgoldenrodyellow')  # Slider for θ
        ax_psi = plt.axes([0.2, 0.01, 0.3, 0.03], facecolor='lightgoldenrodyellow')  # Slider for ψ

        self.slider_phi = Slider(ax_phi, 'Phi (°)', 0, 360, valinit=0)
        self.slider_theta = Slider(ax_theta, 'Theta (°)', 0, 360, valinit=0)
        self.slider_psi = Slider(ax_psi, 'Psi (°)', 0, 360, valinit=0)

        # Connect sliders to the update function
        self.slider_phi.on_changed(update_plot)
        self.slider_theta.on_changed(update_plot)
        self.slider_psi.on_changed(update_plot)

    def add_lab_k_widgets(self, update_lab_k):
        pos = plt.axes([0.95, 0.09, 0.025, 0.03])
        self.lab_k_vec_slider = Slider(pos, 'K scale', 0.05, 2, valinit=1.0, valstep = 0.05)
        self.lab_k_vec_slider.on_changed(update_lab_k)
    def add_sample_scale_widget(self, update_sample_scale):
        pos = plt.axes([0.95, 0.05, 0.025, 0.03])
        self.sample_slider = Slider(pos, 'Sample scale', 0.1, 3, valinit=1.0, valstep = 0.1)
        self.sample_slider.on_changed(update_sample_scale)
    def add_gonio_ring_scale_widget(self, update_gr_scale):
        pos = plt.axes([0.95, 0.01, 0.025, 0.03])
        self.ring_slider = Slider(pos, 'Ring scale', 0.1, 3, valinit=1.0, valstep = 0.1)
        self.ring_slider.on_changed(update_gr_scale)
    def add_gonio_axis_scale_widget(self, update_ga_scale):
        pos = plt.axes([0.85, 0.01, 0.025, 0.03])
        self.ga_slider = Slider(pos, 'GA scale', 0.1, 5, valinit=3.0, valstep = 0.1)
        self.ga_slider.on_changed(update_ga_scale)