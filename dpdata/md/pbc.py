import numpy as np


def posi_diff(box, r0, r1) :
    rbox = np.linalg.inv(box)
    rbox = rbox.T
    p0 = (np.dot(rbox, r0))
    p1 = (np.dot(rbox, r1))
    dp = p0 - p1
    shift = np.zeros(3)
    for dd in range(3) :
        if dp[dd] >= 0.5 : 
            dp[dd] -= 1
        elif dp[dd] < -0.5 :
            dp[dd] += 1
    dr = np.dot(box.T, dp)    
    return dr


def posi_shift(box, r0, r1) :
    rbox = np.linalg.inv(box)
    rbox = rbox.T
    p0 = (np.dot(rbox, r0))
    p1 = (np.dot(rbox, r1))
    dp = p0 - p1
    shift = np.zeros(3)
    for dd in range(3) :
        if dp[dd] >= 0.5 : 
            shift[dd] -= 1
        elif dp[dd] < -0.5 :
            shift[dd] += 1
    return shift


def dir_coord(coord, box) :
    rbox = np.linalg.inv(box)
    return np.matmul(coord, rbox)


def system_pbc_shift(system) :
    f_idx = 0
    prev_ncoord = dir_coord(system['coords'][f_idx], 
                            system['cells' ][f_idx])
    shifts = np.zeros([system.get_nframes(), system.get_natoms(), 3], dtype = int)
    curr_shift = np.zeros([system.get_natoms(), 3], dtype = int)
    half = 0.5 * np.ones([system.get_natoms(), 3])
    for ii in range(system.get_nframes()) :
        curr_ncoord = dir_coord(system['coords'][ii], 
                                system['cells' ][ii])
        diff_ncoord = curr_ncoord - prev_ncoord
        curr_shift -= (diff_ncoord > half)
        curr_shift += (diff_ncoord <-half)
        shifts[ii] = np.copy(curr_shift)
        prev_ncoord = curr_ncoord
    return np.array(shifts, dtype = int)
