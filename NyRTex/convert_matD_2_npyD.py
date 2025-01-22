import scipy.io
import numpy as np
from pathlib import Path

dir = Path(r"C:\Users\kcd17618\Documents\NyRTex\Cu_bolt_Forbes\Cu_bolt_Forbes")
npdir = str(dir) + "_npy"

if not Path(npdir).exists():
    Path(npdir).mkdir()

files = list(dir.walk())[0][-1]

for f in files:

    d_space_mat_file = scipy.io.loadmat(f"{str(dir)}/{f}")

    data = np.array(d_space_mat_file['w2d'][0][0][2])
    data_ax1 = np.array(d_space_mat_file['w2d'][0][0][5][0])[:-1]
    data_ax2 = np.array(d_space_mat_file['w2d'][0][0][8][0])

    dspace = np.repeat(data_ax1[:,None], data.shape[1], axis = 1)

    np.save(npdir+"/"+f.replace(".mat", ".npy"), np.concatenate((dspace[:,:,None], data[:,:,None]), axis = 2))
#np.save(npdir+"/"+"dspacing.npy", data_ax1[:-1])