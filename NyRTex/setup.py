import scipy.io
import numpy as np
import matplotlib.pyplot as plt

# Load the .mat file
sample_mat_file = scipy.io.loadmat(r"C:\Users\kcd17618\Documents\NyRTex\DDSteel\DDSteel_info\mysample_Fe.mat")

savdir = "C:/Users/kcd17618/Documents/NyRTex/DDSteel/DDSteel_info"

# Access the variables in the .mat file
sample_info = dict(zip((0,1,2,3,4,5,6), sample_mat_file['sam'][-1,0]))

crystallography = dict(zip(('name',1,2,3,4,5,6,'cell_parameters','space_group',9), sample_info[0][0][0]))

shape, shape_parameters = sample_info[3][0][0]

print(sample_info, crystallography)

sample_axes = sample_info[4]

np.save(f"{savdir}/shape.npy", shape)
np.save(f"{savdir}/shape_params.npy", shape_parameters)
np.save(f"{savdir}/sample_ax.npy", sample_axes)
np.save(f"{savdir}/cell_params.npy", crystallography['cell_parameters'])
