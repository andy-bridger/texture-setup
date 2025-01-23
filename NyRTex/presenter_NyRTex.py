import sys
import os

# Add the parent directory to the system path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from presenter import Presenter
from helper_funcs import *

class NyrtexPresenter(Presenter):
    def __init__(self, Model, View):
        super().__init__(Model, View)

    def init_view(self):
        self.View.setup_view()
        self.View.update_view_axes()
        self.View.add_lab_k_widgets(self.lab_k_update)
        self.View.add_sample_scale_widget(self.sample_scale_update)
        self.View.add_gonio_ring_scale_widget(self.gr_scale_update)
        self.View.add_gonio_axis_scale_widget(self.ga_scale_update)
        self.View.add_probe_pos_widget(self.update_probe_pos,
                                       self.Model.ewald_radii[0],
                                       self.Model.ewald_radii[1],
                                       self.Model.q_probe)

    def plot_all(self):
        # Lab Frame
        self.plot_lab_sample()
        self.plot_lab_beam_paths()
        self.plot_lab_Ks()
        self.plot_lab_components()
        self.plot_goniometer()
        self.plot_readouts()

        # Pole Detector Frame
        self.plot_pole_figure_detectors()
        self.plot_q_probe()

        # Pole Figure Frame

        self.View.fix_aspect()