import numpy as np
import matplotlib.pyplot as plt
from mesh_utils import *
from stl import mesh
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

infodir = "C:/Users/kcd17618/Documents/NyRTex/Cu_bolt_Forbes/Cu_bolt_info"

shape = np.load(f"{infodir}/shape.npy")
shape_params = np.load(f"{infodir}/shape_params.npy")[0]/1000 # as shape_params are in mm

print(shape, shape_params)

if shape[0] == 'Cilindrical':

    def gen_func(height, radius, n_rad_step = 20, n_vert_steps = 10, offset = np.zeros(3)):

        mesh_top = mesh_circ(radius, offset + np.array((0, 0, height/2)), n_rad_step)
        mesh_bottom = mesh_circ(radius, offset + np.array((0, 0, -height/2)), n_rad_step)
        mesh_sides = mesh_cylindrical_side(radius, height, n_rad_step, n_vert_steps, offset+ np.array((0, 0, -height/2)))
        return np.concatenate((mesh_top, mesh_sides, mesh_bottom), axis =0)



npy_mesh = gen_func(*shape_params)
np.save(f"{infodir}/mesh.npy", npy_mesh)
poly3d = Poly3DCollection(npy_mesh,  edgecolor = 'grey', facecolor = 'white', lw = 0.1)

num_triangles=len(npy_mesh)
data = np.zeros(num_triangles, dtype=mesh.Mesh.dtype)
for i in range(num_triangles):
    #I did not know how to use numpy-arrays in this case. This was the major roadblock
    # assign vertex co-ordinates to variables to write into mesh
    data["vectors"][i] = npy_mesh[i]
m=mesh.Mesh(data)
m.save(f"{infodir}/mesh.stl")

fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection = '3d')
ax.add_collection3d(poly3d)
ax.set_aspect('equal')
plt.show()