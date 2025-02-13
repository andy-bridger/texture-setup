import sys
import os

# Add the parent directory to the system path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


from mantex_run import GetAllPoleFigures, PoleFigurePlot, PoleFigurePresenter
from model import Mantex
from source import Source
from sample import Sample, SampleObject
from goniometer import GenericStateMatrixProvider, Goniometer
import numpy as np
import matplotlib.pyplot as plt
from experiment import ExperimentalData
from detector import Detector
from stl import mesh
from view_NyRTex import NyrtexView
from presenter_NyRTex import NyrtexPresenter
from model_NyRTex import NyrtexMantex
from NyRTex_helpers import *


def project_onto_cylinder(arr):
    xyarr = arr[:,:2]
    nxyarr = xyarr/np.linalg.norm(xyarr, axis = 1)[:,None]
    xs, ys = nxyarr.T
    theta = np.ones_like(xs)*np.pi*0.5
    nperp = xs!=0
    theta[nperp] = np.arctan(ys[nperp]/xs[nperp])#+ 0.25*np.pi*(np.sign(ys[nperp])+1)
    theta = np.where(theta < 1, theta+np.pi, theta)
    return np.concatenate((theta[:,None], arr[:,2][:,None]), axis = 1)


def get_det_grid_masks(det_pos, n_tbins, n_zbins):
    proj = project_onto_cylinder(det_pos.copy())
    all_zbins = np.unique(np.round(proj[:,1], 1))
    z_bin_masks = []
    for zbin in np.array_split(all_zbins, n_zbins):
        mask = np.zeros_like(proj[:,1], dtype = bool)
        for z in zbin:
            mask += np.where(np.abs(proj[:,1] - z) < 0.05,True,False)
        z_bin_masks.append(mask)

    bin_masks = []
    for zm in z_bin_masks:
        t_argsort = np.argsort(proj[zm][:,0])
        for t_bin in np.array_split(t_argsort, n_tbins):
            t_mask = np.zeros_like(proj[:,0], dtype=bool)
            bin_mask = np.zeros_like(t_argsort, dtype=bool)
            bin_mask[t_bin] = True
            t_mask[zm] = bin_mask
            bin_masks.append(t_mask)
    return np.asarray(bin_masks)

n_zbins = 5
n_tbins = 1

info_dir = r"C:\Users\kcd17618\Documents\NyRTex\DDSteel\DDSteel_info"

north_detectors_pos = np.roll(np.load(f"{info_dir}/detector_positions/north_pos.npy"), 1)
south_detectors_pos = np.roll(np.load(f"{info_dir}/detector_positions/south_pos.npy"),1)

n_proj = project_onto_cylinder(north_detectors_pos.copy())
s_proj = project_onto_cylinder(south_detectors_pos.copy())

ngrid_masks = get_det_grid_masks(north_detectors_pos, n_tbins, n_zbins)
sgrid_masks = get_det_grid_masks(south_detectors_pos, n_tbins, n_zbins)

np.save(f"{info_dir}/detector_positions/north_bins_{n_tbins}x{n_zbins}.npy", ngrid_masks)
np.save(f"{info_dir}/detector_positions/south_bins_{n_tbins}x{n_zbins}.npy", sgrid_masks)

