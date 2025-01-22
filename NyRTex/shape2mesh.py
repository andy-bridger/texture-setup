import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from mesh_utils import *

infodir = "C:/Users/kcd17618/Documents/NyRTex/DDSteel_info"

shape = np.load(f"{infodir}/shape.npy")
shape_params = np.load(f"{infodir}/shape_params.npy")[0]/1000 # as shape_params are in mm

if shape[0] == 'Cilindrical':

    def gen_func(radius, height, n_rad_step = 20, n_vert_steps = 10, offset = np.zeros(3)):

        mesh_top = mesh_circ(radius, offset + np.array((0, 0, height/2)), n_rad_step)
        mesh_bottom = mesh_circ(radius, offset + np.array((0, 0, -height/2)), n_rad_step)
        mesh_sides = mesh_cylindrical_side(radius, height, n_rad_step, n_vert_steps, offset+ np.array((0, 0, -height/2)))
        return np.concatenate((mesh_top, mesh_sides, mesh_bottom), axis =0)

test = gen_func(*shape_params)


fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection='3d')
ax.add_collection3d(mplot3d.art3d.Poly3DCollection(test, edgecolor = 'grey', facecolor = 'white', lw = 0.1))
plt.show()

