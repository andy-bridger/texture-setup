from concurrent.futures.process import ProcessPoolExecutor

import numpy as np

from model import Mantex
from goniometer import Goniometer
from experiment import ExperimentalData
import matplotlib.pyplot as plt
from helper_funcs import *
from concurrent.futures import ThreadPoolExecutor
import time


class GetPoleFigure:
    def __init__(self, exp_name, source, sample_pos, sample_view_axis = (1,0,0), exp_fp = "./experiments",
               q_probe = 1, probe_window = 0.05):
        self.exp_data = ExperimentalData(None, None, from_data= True)
        self.exp_data.from_data(exp_name, exp_fp)
        self.sample_view_axis = np.asarray(sample_view_axis)
        self.goniometer_arr = [Goniometer(phi = init_gonio[0], theta = init_gonio[1], psi = init_gonio[2],
                                     exp_runs=()) for init_gonio in self.exp_data.goniometer_positions]

        self.mantex_arr = [Mantex(source, self.exp_data.detectors, np.asarray(sample_pos),
                                  goniometer, sample_view_axis,
                                  q_probe, probe_window) for goniometer in self.goniometer_arr]
    def execute(self):
        results = map(self.dispatch_pfi_calc, self.mantex_arr, self.exp_data.detector_readouts)
        self.pole_figure_intensities = np.concatenate([r for r in results], axis = 0)
        self.pole_figure_intensities = self.pole_figure_intensities[np.argsort(self.pole_figure_intensities[:,-1])]

    def dispatch_pfi_calc(self, mantex, spectra):
        return mantex.get_pole_figure_intensities(spectra)

class PoleFigurePlot:
    def __init__(self):
        pass

    def setup_view(self):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(1,1,1)

    def fix_aspect(self, ax):
        ax.set_aspect('equal')

    def plot_pole_figure(self, ax, sample_view_axis, pole_figure_intensities, eq):
        artists = []
        if np.where(sample_view_axis == 1)[0] == 0:
            artists.append(ax.scatter(pole_figure_intensities[:, 1], pole_figure_intensities[:, 2],
                                      c=pole_figure_intensities[:, 3]))
        if np.where(sample_view_axis == 1)[0] == 1:
            artists.append(ax.scatter(pole_figure_intensities[:, 0], pole_figure_intensities[:, 2],
                                      c=pole_figure_intensities[:, 3]))
        if np.where(sample_view_axis == 1)[0] == 2:
            artists.append(ax.scatter(pole_figure_intensities[:, 0], pole_figure_intensities[:, 1],
                                      c=pole_figure_intensities[:, 3]))
        artists.append(ax.plot(eq[0], eq[1], c='grey'))
        return artists

class PoleFigurePresenter:
    def __init__(self, model, view):
        self.Model = model
        self.View = view
    def init_view(self):
        self.View.setup_view()
    def plot_pole_figure(self):
        self.init_view()
        self.Model.execute()
        self.View.plot_pole_figure(self.View.ax, self.Model.sample_view_axis, self.Model.pole_figure_intensities, equator())
        self.View.fix_aspect(self.View.ax)



