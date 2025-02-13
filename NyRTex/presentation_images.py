from stl import mesh
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import math


def cart2sph(x, y, z):
    xy = np.sqrt(x ** 2 + y ** 2)  # sqrt(x² + y²)

    x_2 = x ** 2
    y_2 = y ** 2
    z_2 = z ** 2

    r = np.sqrt(x_2 + y_2 + z_2)  # r = sqrt(x² + y² + z²)

    theta = np.arctan2(y, x)

    phi = np.arctan2(xy, z)

    return r, theta, phi

def point2col(p):
    norm_p = p/np.linalg.norm(p)
    x,y,z = (norm_p/2)+0.5
    return np.array((x,y,z, 0.9999))


def stereographic_projection(rot_mat, point):
    # as mentioned before, we want to use our to_pole_view to reproject the points
    spoint = (rot_mat @ point)
    # we want the vector to whichever pole is on the opposite hemisphere
    oppo_pole = np.array((0, 0, -np.sign(spoint[2])))
    svec = spoint - oppo_pole

    if spoint[2] == 0:  # point is within the projection plane
        return point, spoint, spoint
    else:
        u = -oppo_pole[2] / svec[2]
        proj_point = (u * svec) + oppo_pole
        return (np.linalg.inv(rot_mat) @ proj_point), proj_point, spoint

def fibonacci_sphere(samples=1000, r=1):

    points = []
    phi = math.pi * (math.sqrt(5.) - 1.)  # golden angle in radians

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = math.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = math.cos(theta) * radius
        z = math.sin(theta) * radius

        points.append((x*r, y*r, z*r))

    return points

def fig1():
    info_dir = r"C:\Users\kcd17618\Documents\NyRTex\DDSteel\DDSteel_info"
    mesh_vecs = mesh.Mesh.from_file(f"{info_dir}/mesh.stl").vectors



    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection = "3d")
    ax.add_collection3d(mplot3d.art3d.Poly3DCollection(mesh_vecs, edgecolor = 'grey', facecolor = 'white', lw = 0.1))
    mat = np.eye(3,3)
    length = np.linalg.norm(mesh_vecs[:,0,:], axis = 1).max() * 1
    pos = np.ones(3)*length*3
    ax.quiver(*pos, *mat[:,0]*length, color = 'cadetblue')
    ax.quiver(*pos, *mat[:,1]*length, color = 'magenta')
    ax.quiver(*pos, *mat[:,2]*length, color = 'dodgerblue')
    for p in fibonacci_sphere(500, 0.02):
        ax.quiver(*np.zeros(3), *p, alpha = 0.1, arrow_length_ratio = 0.1, color = point2col(p))
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_xlim([-4*length, 4*length])
    ax.set_ylim([-4*length, 4*length])
    ax.set_zlim([-4*length, 4*length])
    plt.show()

def fig2():
    info_dir = r"C:\Users\kcd17618\Documents\NyRTex\DDSteel\DDSteel_info"
    mesh_vecs = mesh.Mesh.from_file(f"{info_dir}/mesh.stl").vectors

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection="3d")
    mat = np.eye(3, 3)
    length = np.linalg.norm(mesh_vecs[:, 0, :], axis=1).max() * 1
    pos = np.ones(3) * length * 3
    ax.quiver(*pos, *mat[:, 0] * length, color='cadetblue')
    ax.quiver(*pos, *mat[:, 1] * length, color='magenta')
    ax.quiver(*pos, *mat[:, 2] * length, color='dodgerblue')
    for p in fibonacci_sphere(500, 0.02):
        sp = stereographic_projection(np.eye(3,3), p/np.linalg.norm(p))[0]*0.02
        ax.scatter(*sp, alpha=0.5, color=point2col(p))
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_xlim([-4 * length, 4 * length])
    ax.set_ylim([-4 * length, 4 * length])
    ax.set_zlim([-4 * length, 4 * length])
    plt.show()

fig1()