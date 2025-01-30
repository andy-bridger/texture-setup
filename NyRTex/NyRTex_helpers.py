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