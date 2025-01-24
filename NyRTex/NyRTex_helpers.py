import numpy as np

def bank_dets(dat):
    return dat.sum(axis = 1)[None,:]

def load_det_red(fp, Q):
    dat = np.load(fp) #(Q_pos, det_pix)
    dat = bank_dets(dat) #(det_banks, Q)
    fd =  combine_with_Q(dat, Q) #(det_banks, Q, Qval/I)
    return fd


def combine_with_Q(dat, Q):
    return np.concatenate((Q[None,:,None], dat[:,:,None]), axis = -1)