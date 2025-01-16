import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
from model import Model
from goniometer import Goniometer
from detector import Detector
from source import Source
from sample import Sample
from presenter import Presenter
from view import View

from matplotlib.widgets import Slider, CheckButtons

# Create instances of Source, Sample, and SphereConstructions
ndet = Detector([0, 10, 0], 'north')
ndet2 = Detector([5*np.sqrt(2), 5*np.sqrt(2), 0], 'north2')
sdet = Detector([0, -10, 0], 'south')
source = Source([20, 0 , 0])
sample = Sample([0,0,0], [0,0,0],
                (1,1,1), 2,
                [(0,0,0)], q_probe=1,
                cell_colors=['purple'])

exp_runs = [(phi, theta, 0) for theta in np.arange(0, 360, 45) for phi in np.arange(0, 360, 45) ]

goniometer = Goniometer(exp_runs = exp_runs)

model = Model(source,
              [ndet,ndet2, sdet],
              sample,
              goniometer,
              (0.7, 1.4),
              1,
              (1,0,0),
              detector_colors=('red','green', 'blue'), ewald_steps=100)

view = View()

presenter = Presenter(model, view)
presenter.plot_all()
plt.show()

