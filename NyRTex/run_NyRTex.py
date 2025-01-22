import sys
import os

# Add the parent directory to the system path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


from mantex_run import GetAllPoleFigures, PoleFigurePlot, PoleFigurePresenter
from source import Source
from sample import Sample, SampleObject
from goniometer import GenericStateMatrixProvider, Goniometer
import numpy as np
import matplotlib.pyplot as plt
from experiment import ExperimentalData
from detector import Detector
from stl import mesh

info_dir = r"C:\Users\kcd17618\Documents\NyRTex\Cu_bolt_Forbes\Cu_bolt_info"
det_r_dir = r"C:\Users\kcd17618\Documents\NyRTex\Cu_bolt_Forbes\Cu_bolt_Forbes_npy"

north_detectors_pos = np.load(f"{info_dir}/detector_positions/north_pos.npy")
south_detectors_pos = np.load(f"{info_dir}/detector_positions/south_pos.npy")

ndp = np.mean(north_detectors_pos, axis = 0)
sdp = np.mean(south_detectors_pos, axis = 0)
print(ndp, sdp)

detectors = [Detector(ndp, 'north'), Detector(sdp, 'South')]

with open(f"{info_dir}/runs_pA.txt", 'r') as f:
    runs = [x.split('\t')[0] for x in f.readlines()]
with open(f"{info_dir}/rotation_pA.txt", 'r') as f:
    angs = [x.rstrip('\n').split('\t') for x in f.readlines()]

drs = [np.moveaxis(np.concatenate((np.load(f"{det_r_dir}/ENOR00{run}.npy").sum(axis = 1)[:,None,:],
                                   np.load(f"{det_r_dir}/ESUR00{run}.npy").sum(axis = 1)[:,None,:]), axis = 1), 1,0) for run in runs]

exp_data = ExperimentalData(detectors, 'Cu_bolt', from_data= False)
exp_data.detector_readouts = drs
exp_data.smps = [Goniometer(1, int(ang[1]), int(ang[1]) ,int(ang[2]), scheme = 'rot1' ) for ang in angs]

source = Source(np.load(f"{info_dir}/source_pos.npy"))
sample = SampleObject((0.,0.,0.), GenericStateMatrixProvider(np.eye(3), np.zeros(3)), mesh = mesh.Mesh.from_file(f"{info_dir}/mesh.stl"))
sample_view_axes = ((1,0,0), (0,1,0))


alg = GetAllPoleFigures(exp_data, source, sample, sample_view_axes)

view = PoleFigurePlot()
presenter = PoleFigurePresenter(alg, view)

presenter.plot_pole_figure()
plt.show()