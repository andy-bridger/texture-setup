from mantex_run import GetAllPoleFigures, PoleFigurePlot, PoleFigurePresenter
from source import Source
from sample import Sample, SampleObject
from goniometer import GenericStateMatrixProvider
import numpy as np
import matplotlib.pyplot as plt
from experiment import ExperimentalData
from stl import mesh

exp_data = ExperimentalData(None, None, from_data= True)
exp_data.from_data('test')

source = Source([20, 0 , 0])
sample = SampleObject((0.,0.,0.), GenericStateMatrixProvider(np.eye(3), np.zeros(3)), mesh = mesh.Mesh.from_file("./shapes/utah_teapot.stl"))
sample_view_axes = ((1,0,0), (0,1,0))

alg = GetAllPoleFigures(exp_data, source, sample, sample_view_axes)

## TODO:
# Fix mesh scale and movement system. The mesh position should be described using a SMP relative to the centre of the gauge volume at sample_position.
# add a proof of concept attenuation system so the intensity is somewhat related to the mesh
# this will probably just be some sort of path length thing

view = PoleFigurePlot()
presenter = PoleFigurePresenter(alg, view)

presenter.plot_pole_figure()
plt.show()