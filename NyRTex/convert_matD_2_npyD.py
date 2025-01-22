import scipy.io
import numpy as np
from pathlib import Path

dir = Path(r'C:\Users\kcd17618\Documents\NyRTex\DDSteel_d\DDSteel_d')
npdir = str(dir) + "_npy"

if not Path(npdir).exists():
    Path(npdir).mkdir()

files = list(dir.walk())[0][-1]

for f in files:

    d_space_mat_file = scipy.io.loadmat(f"{str(dir)}/{f}")

    data = np.array(d_space_mat_file['w2d'][0][0][2])
    data_ax1 = np.array(d_space_mat_file['w2d'][0][0][5][0])
    data_ax2 = np.array(d_space_mat_file['w2d'][0][0][8][0])

    np_data = np.concatenate((data_ax1[:-1], data.sum(axis = 1))) # I think we have summed across the time axis and have this in dspace vs counts-ish

    np.save(npdir+"/"+f.replace(".mat", ".npy"), np_data)