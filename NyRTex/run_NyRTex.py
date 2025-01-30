import sys
import os

# Add the parent directory to the system path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


from mantex_run import GetAllPoleFigures, PoleFigurePlot, PoleFigurePresenter
from model import Mantex
from source import Source
from sample import Sample, SampleObject
from goniometer import GenericStateMatrixProvider, Goniometer
import numpy as np
import matplotlib.pyplot as plt
from experiment import ExperimentalData
from detector import Detector
from stl import mesh
from view_NyRTex import NyrtexView
from presenter_NyRTex import NyrtexPresenter
from model_NyRTex import NyrtexMantex
from NyRTex_helpers import *

info_dir = r"C:\Users\kcd17618\Documents\NyRTex\Cu_bolt_Forbes\Cu_bolt_info"
det_r_dir = r"C:\Users\kcd17618\Documents\NyRTex\Cu_bolt_Forbes\Cu_bolt_Forbes_npy"

north_detectors_pos = np.load(f"{info_dir}/detector_positions/north_pos.npy")
south_detectors_pos = np.load(f"{info_dir}/detector_positions/south_pos.npy")

ndp = np.mean(north_detectors_pos, axis = 0)
sdp = np.mean(south_detectors_pos, axis = 0)

detectors = [Detector(ndp, 'north'), Detector(sdp, 'South')]

with open(f"{info_dir}/runs_pA.txt", 'r') as f:
    runs = [x.split('\t')[0] for x in f.readlines()]
with open(f"{info_dir}/rotation_pA.txt", 'r') as f:
    angs = [x.rstrip('\n').split('\t') for x in f.readlines()]


sig_ax = np.load(f"{det_r_dir}/dspacing.npy")


drs = np.asarray([np.concatenate((load_det_red(f"{det_r_dir}/ENOR00{run}.npy", sig_ax),
                                   load_det_red(f"{det_r_dir}/ESUR00{run}.npy", sig_ax)), axis = 0) for run in runs])
# run,det, Q, Qval/I

exp_data = ExperimentalData(detectors, 'Cu_bolt', from_data= False)
exp_data.detector_readouts = drs #(run, det, reading_num,Q/sig)
exp_data.smps = [Goniometer(1, int(ang[0]), int(ang[1]) ,int(ang[2]), scheme = 'rot1' ) for ang in angs]

source = Source(np.load(f"{info_dir}/source_pos.npy"))
sample = SampleObject((0.,0.,0.), GenericStateMatrixProvider(np.eye(3), np.zeros(3)),
                      mesh = mesh.Mesh.from_file(f"{info_dir}/mesh.stl"), mesh_scale=1)
sample_view_axes = ((1,0,0), (0,0,1))

q_probe = 2.08

nview = NyrtexView()
mantex = NyrtexMantex(source, detectors, sample,
                                  exp_data, sample_view_axes, ['red', 'blue'], q_probe=q_probe, probe_window=0.1)
presenter = NyrtexPresenter(mantex, nview)

presenter.plot_all()
plt.show()