import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
from funcs import *

from matplotlib.widgets import Slider

# Create instances of Source, Sample, and SphereConstructions
ndet = Detector([0, 10, 0], 'north')
ndet2 = Detector([1, 11, 0], 'north2')
sdet = Detector([0, -10, 0], 'south')
source = Source([10, 0 , 0])
sample = Sample([0,0,0], [0,0,0],
                (1,1,1), 2,
                [(0,0,0),(45,0,0)], q_probe=1.41,
                cell_colors=['pink', 'purple'])

sphere_construction = SphereConstructions(source, [ndet,ndet2, sdet],
                                          sample, (1, 2),1,
                                          (1,0,0),
                                          detector_colors=('red','green', 'blue'))

# Function to update the visualization when sliders change
def update_plot(val):
    # Update the sample's orientation
    sample.orient_array = [slider_phi.val, slider_theta.val, slider_psi.val]
    sample.update()
    sphere_construction.update()
    sphere_construction.plot_all()

# Setup the main figure and axes for visualization
sphere_construction.plot_all()

# Create sliders below the visualization
fig = sphere_construction.fig
ax_phi = plt.axes([0.2, 0.01, 0.65, 0.03], facecolor='lightgoldenrodyellow')  # Slider for φ
ax_theta = plt.axes([0.2, 0.05, 0.65, 0.03], facecolor='lightgoldenrodyellow')  # Slider for θ
ax_psi = plt.axes([0.2, 0.09, 0.65, 0.03], facecolor='lightgoldenrodyellow')  # Slider for ψ

slider_phi = Slider(ax_phi, 'Phi (°)', 0, 360, valinit=0)
slider_theta = Slider(ax_theta, 'Theta (°)', 0, 360, valinit=0)
slider_psi = Slider(ax_psi, 'Psi (°)', 0, 360, valinit=0)

# Connect sliders to the update function
slider_phi.on_changed(update_plot)
slider_theta.on_changed(update_plot)
slider_psi.on_changed(update_plot)

plt.show()