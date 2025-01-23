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
        self.lab_ax = self.fig.add_subplot(gs0[:10, :5], projection='3d')
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
    def fix_aspect(self):
        self.lab_ax.set_aspect('equal')
        self.pole_proj_ax.set_aspect('equal')
        self.calc_pf_ax.set_aspect('equal')