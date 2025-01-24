import scipy.io
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

dir = Path(r"C:\Users\kcd17618\Documents\NyRTex\Cu_bolt_Forbes\Cu_bolt_Forbes")
npdir = str(dir) + "_npy"

if not Path(npdir).exists():
    Path(npdir).mkdir()

files = list(dir.walk())[0][-1]

plt.figure()
for f in files:

    d_space_mat_file = scipy.io.loadmat(f"{str(dir)}/{f}")

    data = np.array(d_space_mat_file['w2d'][0][0][2])
    data_ax1 = np.array(d_space_mat_file['w2d'][0][0][5][0])[:-1]
    data_ax2 = np.array(d_space_mat_file['w2d'][0][0][8][0])

    #plt.plot(data.sum(axis = 1))

    np.save(npdir+"/"+f.replace(".mat", ".npy"), data)

#plt.show()
np.save(npdir+"/"+"dspacing.npy", data_ax1)