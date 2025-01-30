import sys
import os

# Add the parent directory to the system path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from view import View
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
from model import *
from mpl_toolkits import mplot3d

from matplotlib.widgets import Slider, Button, TextBox

class NyrtexView(View):
    def __init__(self):
        pass
    def setup_view(self):
        self.fig = plt.figure()
        gs0 = self.fig.add_gridspec(10, 10, wspace = 0.1, hspace= 0.0)
        self.lab_ax = self.fig.add_subplot(gs0[:8, :4], projection='3d')
        #self.lab_ax.view_init(azim=90, elev=-120, roll=180)
        self.sample_ax = self.fig.add_subplot(gs0[2:6, 3:8], projection='3d')
        #self.sample_ax.view_init(azim=90, elev=-120, roll=180)
        self.object_ax = self.fig.add_subplot(gs0[2:6, 7:10], projection='3d')
        #self.object_ax.view_init(azim=90, elev=-120, roll=180)
        self.pole_proj_ax = self.fig.add_subplot(gs0[6:8, 0:5])
        self.det_ax = self.fig.add_subplot(gs0[6:8, 4:7])
        self.calc_pf_ax = self.fig.add_subplot(gs0[6:8, 6:10])
        self.pole_proj_ax.set_xlim([-1.1, 1.1])
        self.pole_proj_ax.set_ylim([-1.1, 1.1])
        self.pole_proj_ax.set_axis_off()
        self.lab_ax.set_axis_off()
    def update_view_axes(self):
        self.pole_proj_ax.clear()
        self.pole_proj_ax.set_axis_off()
        self.calc_pf_ax.set_axis_off()
        self.lab_ax.set_axis_off()
        self.sample_ax.set_axis_off()
        self.object_ax.set_axis_off()
    def fix_aspect(self):
        self.lab_ax.set_aspect('equal')
        self.pole_proj_ax.set_aspect('equal')
        self.calc_pf_ax.set_aspect('equal')
        self.sample_ax.set_aspect('equal')
        self.object_ax.set_aspect('equal')
    def show_detector_readout(self, ax, readouts, cols):
        artists = []
        alpha = np.clip(2/len(readouts), 0, 1)
        for i, r in enumerate(readouts):
            artists.append(ax.plot(r[:,0], (r[:,1])+(i*0.002), c = cols[i], alpha = alpha))
        return artists
    def plot_pole_figure_position(self, ax, pole_figure_points, detector_colors, eq):
        artists = []
        size = 100 / len(detector_colors)
        artists.append(ax.scatter(pole_figure_points[:,0], pole_figure_points[:,1], c = detector_colors, s = size))
        artists.append(ax.plot(eq[0], eq[1], c='grey'))
        return artists

    def plot_pole_figure_intensities(self, ax, pole_figure_intensities, det_pKs, det_cols, eq):
        artists = []
        size = 100/len(det_pKs)
        artists.append(ax.scatter(pole_figure_intensities[:,0], pole_figure_intensities[:,1],
                            c = pole_figure_intensities[:,3], s = size))
        for idet, pK in enumerate(det_pKs):
            artists.append(ax.scatter(pK[0], pK[1], facecolors = 'none',
                                edgecolors=det_cols[idet],s = 1.25*size))
        #eq = self.equator()
        artists.append(ax.plot(eq[0], eq[1], c = 'grey'))
        return artists
    def plot_cart_ax(self, ax, mat, pos, length = 1):
        artists = []
        artists.append(ax.quiver(*pos, *mat[:,0]*length, color = 'cadetblue'))
        artists.append(ax.quiver(*pos, *mat[:,1]*length, color = 'magenta'))
        artists.append(ax.quiver(*pos, *mat[:,2]*length, color = 'dodgerblue'))
        return artists
    def plot_cart_ax_2d(self, ax, mat, pos, length = 1):
        artists = []
        artists.append(ax.quiver(*pos, *mat[:,1]*length, color = 'magenta', scale = 1, width = 0.01))
        artists.append(ax.quiver(*pos, *mat[:,0]*length, color = 'cadetblue', scale = 1, width = 0.01))
        return artists
    def add_probe_pos_widget(self, update_probe_p_scale, min, max, init):
        pos = plt.axes([0.44, 0.2, 0.22, 0.03])
        self.probe_p_slider = Slider(pos, 'probe Q', min, max, valinit=init, valstep = 0.02)
        self.probe_p_slider.on_changed(update_probe_p_scale)
    def add_run_num_widget(self, update_run, min, max, init):
        pos = plt.axes([0.44, 0.1, 0.22, 0.03])
        self.run_slider = Slider(pos, 'run_num', min, max, valinit=init, valstep = 1)
        self.run_slider.on_changed(update_run)
    def add_sample_scale_widget(self, update_sample_scale, init_scale):
        pos = plt.axes([0.85, 0.05, 0.025, 0.03])
        self.sample_slider = Slider(pos, 'Sample scale', 0.00, 5, valinit=init_scale, valstep = 0.1)
        self.sample_slider.on_changed(update_sample_scale)
