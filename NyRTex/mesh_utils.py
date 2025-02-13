import numpy as np

def mesh_circ(radius, offset, n_rad_step):
    us = np.linspace(0, 2 * np.pi, n_rad_step)
    xs = (radius * np.cos(us) + offset[0])[:, None]
    ys = (radius * np.sin(us) + offset[1])[:, None]
    zs = np.ones_like(xs) * offset[2]
    cent_xs = np.ones_like(xs)* offset[0]
    cent_ys = np.ones_like(xs) * offset[1]
    v0 = np.concatenate((xs, ys, zs), axis=1)
    v1 = np.concatenate((cent_xs, cent_ys, zs), axis=1)
    v2 = np.concatenate((np.roll(xs, 1), np.roll(ys, 1), zs), axis=1)
    return np.concatenate([v0[:,None], v1[:,None], v2[:,None]], axis = 1)

def MakeTriangulatedGridIndex(Nr,Nc):

    out = np.empty((Nr-1,Nc-1,2,3),dtype=int)

    r = np.arange(Nr*Nc).reshape(Nr,Nc)

    out[:,:, 0,0] = r[:-1,:-1]
    out[:,:, 1,0] = r[:-1,1:]
    out[:,:, 0,1] = r[:-1,1:]

    out[:,:, 1,1] = r[1:,1:]
    out[:,:, :,2] = r[1:,:-1,None]

    out.shape =(-1,3)
    return out

def mesh_cylindrical_side(radius, height, n_rad_step, n_vert_step, offset = np.zeros(3)):
    face_tris = MakeTriangulatedGridIndex(n_vert_step, n_rad_step)

    zs = np.linspace(0, height, n_vert_step)

    us = np.linspace(0, 2*np.pi, n_rad_step)
    xs = (radius*np.cos(us))[face_tris%n_rad_step]
    ys = (radius*np.sin(us))[face_tris%n_rad_step]
    zs = zs[face_tris//n_rad_step]

    return np.concatenate([xs[:,:,None], ys[:,:,None], zs[:,:,None]], axis = -1) + offset[None, None,:]

def mesh_cuboid(a, b, c, offset = np.zeros(3)):
    vec = np.array((a/2, b/2, c/2))
    p1 = (np.array((-1, -1,-1))*vec)[None,:]
    p2 = (np.array((1,-1,-1))*vec)[None,:]
    p3 = (np.array((1,1,-1))*vec)[None,:]
    p4 = (np.array((-1,1,-1))*vec)[None,:]
    p5 = -p3
    p6 = -p4
    p7 = -p1
    p8 = -p2
    tris = np.asarray((np.concatenate((p1, p2, p3), axis = 0),
            np.concatenate((p1, p4, p3), axis=0),
            np.concatenate((p1, p2, p6), axis=0),
            np.concatenate((p1, p5, p6), axis=0),
            np.concatenate((p2, p3, p7), axis=0),
            np.concatenate((p2, p6, p7), axis=0),
            np.concatenate((p3, p4, p8), axis=0),
            np.concatenate((p3, p7, p8), axis=0),
            np.concatenate((p1, p4, p8), axis=0),
            np.concatenate((p1, p5, p8), axis=0),
            np.concatenate((p5, p6, p7), axis=0),
            np.concatenate((p5, p8, p7), axis=0),))
    return tris + offset[None,None,:]
