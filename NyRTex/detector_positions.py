import numpy as np
import matplotlib.pyplot as plt

det_pos = np.load(r"C:\Users\kcd17618\Documents\NyRTex\DDSteel_info\detector_positions\det_pos.npy")

north_pos = det_pos[det_pos[:,0]>1]
south_pos = det_pos[det_pos[:,0]<-1]

print(north_pos.shape)

fig = plt.figure()
ax=fig.add_subplot(1,1,1,projection = '3d')
ax.scatter(**dict(zip(('xs','ys','zs'), north_pos.T)))
ax.scatter(**dict(zip(('xs','ys','zs'), south_pos.T)))
plt.show()