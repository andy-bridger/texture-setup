from concurrent.futures.process import ProcessPoolExecutor

import numpy as np

from model import Mantex
from goniometer import Goniometer, GenericStateMatrixProvider
from experiment import ExperimentalData
from sample import SampleObject
import matplotlib.pyplot as plt
from helper_funcs import *
from concurrent.futures import ThreadPoolExecutor
from copy import deepcopy


class GetSinglePoleFigure:
    def __init__(self, exp_data, source, sample, sample_view_axes = ((1,0,0), (0,1,0)),
               q_probe = 1, probe_window = 0.05, ind = 0):
        self.exp_data = exp_data
        self.sample = sample

        self.sample_view_axes = np.asarray(sample_view_axes)
        self.smp = exp_data.smps[ind]
        self.sample.set_smp(self.smp)

        self.mantex = Mantex(source, self.exp_data.detectors, self.sample,
                                  self.smp, sample_view_axes,
                                  q_probe, probe_window)
        self.spectra = self.exp_data.detector_readouts[ind]

    def execute(self):
        return self.mantex.get_pole_figure_intensities(self.spectra)

class GetAllPoleFigures:
    def __init__(self, exp_data, source, sample, sample_view_axes = ((1,0,0),(0,1,0)),  q_probe = 1, probe_window = 0.05):
        self.exp_data = exp_data
        self.source = source
        self.sample = sample
        self.sample_view_axes = np.asarray(sample_view_axes)
        self.q_probe = q_probe
        self.probe_window = probe_window
    def gen_spfs(self):
        spfs = []
        for i, smp in enumerate(self.exp_data.smps):
            new_sample = deepcopy(self.sample)
            self.sample.smp = smp
            spfs.append(GetSinglePoleFigure(deepcopy(self.exp_data), self.source, new_sample, self.sample_view_axes,
                                q_probe=self.q_probe, probe_window=self.probe_window, ind=i))
        self.spfs = spfs
    def execute(self):
        self.gen_spfs()
        pfis = []
        for i, spf in enumerate(self.spfs):
            pfis.append(spf.execute())
        pfis = np.concatenate(pfis, axis = 0)
        self.pole_figure_intensities = pfis[np.argsort(pfis[:,-1])]


class PoleFigurePlot:
    def __init__(self):
        pass

    def setup_view(self):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(1,1,1)

    def fix_aspect(self, ax):
        ax.set_aspect('equal')

    def plot_pole_figure(self, ax, pole_figure_intensities, eq):
        artists = []
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
    def calc_pole_figure(self):
        self.Model.execute()
    def plot_pole_figure(self):
        self.init_view()
        self.calc_pole_figure()
        self.View.plot_pole_figure(self.View.ax, self.Model.pole_figure_intensities, equator())
        self.View.fix_aspect(self.View.ax)



