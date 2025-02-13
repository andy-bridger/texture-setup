import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation


def get_mat_from_angles(ang, scheme, prog):
    if scheme == 'euler':
        if prog == 'nyrtex':
            seq = 'ZXZ' #perpendicular to beam and dets, beam, perp
        else:
            seq = 'YZY' #perpendicular to beam and dets, beam, perp
        return Rotation.from_euler(seq, (ang[0], ang[1], ang[2]), degrees=True).as_matrix()

    else:
        if prog == 'nyrtex':
            seq = 'yz'
            ax = np.array((0,0,1))
        else:
            seq = 'xy'
            ax = np.array((0,1,0))

        init_rot = Rotation.from_euler(seq, (ang[0], ang[1]), degrees=True)
        final_rotvec = init_rot.apply(ax * (ang[2] * np.pi / 180))
        second_rot = Rotation.from_rotvec(final_rotvec, degrees=False)
        return second_rot.as_matrix() @ init_rot.as_matrix()

prog = 'mantid' # nyrtex or mantid, they have different reference frames

fp = r"C:\Users\kcd17618\Documents\NyRTex\DDSteel\DDSteel_info\rotations.txt"
scheme = 'euler'
save_path = fr"C:\Users\kcd17618\Documents\NyRTex\DDSteel\DDSteel_info\state_transforms_{prog}_{scheme}.npy"
translations = None

with open(fp, "r") as f:
    angles = [[float(y) for y in x.rstrip('\n').split('\t')] for x in f.readlines()]

if type(translations) != type(None):
    translations = np.load(translations) # load from a file
else:
    translations = np.zeros((len(angles), 3)) # add null translation


exp_mats = np.asarray([get_mat_from_angles(ang, scheme, prog).reshape(-1) for ang in angles])

print(exp_mats.shape)

np.save(save_path, np.concatenate([exp_mats, translations], axis = 1))

