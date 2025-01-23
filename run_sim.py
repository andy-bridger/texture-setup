import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
from model import Mantex
from mantex_sim import MantexSim
from goniometer import Goniometer
from detector import Detector
from source import Source
from sample import Sample
from presenter import Presenter
from view import View
from stl import mesh


from matplotlib.widgets import Slider, CheckButtons

# Create instances of Source, Sample, and SphereConstructions
ndet = Detector([0, 10, 0], 'north')
ndet2 = Detector([5*np.sqrt(2), 5*np.sqrt(2), 0], 'north2')
sdet = Detector([0, -10, 0], 'south')
source = Source([20, 0 , 0])

# default runs

#exp_runs = [(phi, theta, 0) for theta in np.arange(0, 360, 45) for phi in np.arange(0, 360, 45) ]

#test loading some real runs

with open(f"C:/Users/kcd17618/Documents/NyRTex/DDSteel/DDSteel_info/rotations.txt", 'r') as f:
    angs = [x.rstrip('\n').split('\t') for x in f.readlines()]

exp_runs = [(int(ang[0]), int(ang[1]) ,int(ang[2])) for ang in angs]

goniometer = Goniometer(exp_runs = exp_runs, offset=np.array((0,0,-7)), scheme='rot1')

sample = Sample([0,0,0], goniometer,
                (1,1,1), 2,
                [(0,0,0)], q_probe=1,
                cell_colors=['purple'], mesh=mesh.Mesh.from_file("./shapes/utah_teapot.stl"))


model = MantexSim(source,
              [ndet, sdet],
              sample,
              goniometer,
              (0.7, 1.4),
                  ((0,1,0), (1,0,0)),
              detector_colors=('red','blue'), ewald_steps=100)

view = View()

presenter = Presenter(model, view)
presenter.plot_all()
plt.show()

