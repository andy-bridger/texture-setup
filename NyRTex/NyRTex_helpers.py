import numpy as np

def bank_dets(dat, bank_masks):
    if type(bank_masks) == type(None):
        return dat.sum(axis = 1)[None,:]
    else:
        banked = np.asarray([dat.T[bm].sum(axis = 0) for bm in bank_masks])
        return banked

def load_det_red(fp, Q, bank_masks=None):
    dat = np.load(fp) #(Q_pos, det_pix)
    dat = bank_dets(dat, bank_masks) #(det_banks, Q)
    fd =  combine_with_Q(dat, Q) #(det_banks, Q, Qval/I)
    return fd


def combine_with_Q(dat, Q):
    rQ = np.repeat(Q[None, :], dat.shape[0], axis = 0)
    return np.concatenate((rQ[:,:,None], dat[:,:,None]), axis = -1)

def convert_pfi_to_polar(pfi):
    x, y = pfi[:,0], pfi[:,1]
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    out_pfi = pfi.copy()
    out_pfi[:,0] = rho
    out_pfi[:,1] = phi
    return out_pfi