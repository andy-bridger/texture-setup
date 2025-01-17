from mantex_run import GetPoleFigure, PoleFigurePlot, PoleFigurePresenter
from source import Source
from sample import Sample, SampleObject
import numpy as np
import matplotlib.pyplot as plt

source = Source([20, 0 , 0])
exp_runs = [(phi, theta, 0) for theta in np.arange(0, 360, 45) for phi in np.arange(0, 360, 45) ]
sample = SampleObject((0.,0.,0.))
sample_view_axis = (1,0,0)

alg = GetPoleFigure("test3", source, sample, sample_view_axis)
view = PoleFigurePlot()
presenter = PoleFigurePresenter(alg, view)

presenter.plot_pole_figure()
plt.show()