import numpy as np
from .pbc import system_pbc_shift

def _msd(coords, cells, pbc_shift, begin):
    nframes = cells.shape[0]
    natoms = coords.shape[1]
    ff = begin
    prev_coord = coords[ff] + np.matmul(pbc_shift[ff], cells[ff])
    msds = [0.]
    for ff in range(begin+1,nframes) :
        curr_coord = coords[ff] + np.matmul(pbc_shift[ff], cells[ff])
        diff_coord = curr_coord - prev_coord
        msds.append(np.sum(diff_coord * diff_coord) / natoms)
    return np.array(msds)

def _msd_win(coords, cells, pbc_shift, begin, window):
    nframes = cells.shape[0]
    natoms = coords.shape[1]
    ncoords = np.zeros(coords.shape)
    msd = np.zeros([window])
    for ff in range(nframes) :
        ncoords[ff] = coords[ff] + np.matmul(pbc_shift[ff], cells[ff]) 
    cc = 0
    for ii in range(begin,nframes-window+1) :
        start = np.tile(ncoords[ii], (window, 1, 1))
        diff_coord = ncoords[ii:ii+window] - start
        diff_coord = np.reshape(diff_coord, [-1, natoms * 3])
        msd += np.sum(diff_coord * diff_coord, axis = 1) / natoms
        cc += 1
    return np.array(msd) / cc

def msd(system, sel = None, begin = 0, window = 0) :
    natoms = system.get_natoms()    
    if sel is None :
        sel_idx = np.arange(natoms)
    else :
        sel_idx = []
        for ii in range(natoms) :
            if sel[ii] :
                sel_idx.append(ii)
        sel_idx = np.array(sel_idx, dtype = int)
    nsel = sel_idx.size
    nframes = system.get_nframes()
    pbc_shift = system_pbc_shift(system)
    coords = system['coords']
    cells = system['cells']
    pbc_shift = pbc_shift[:,sel_idx,:]
    coords = coords[:,sel_idx,:]
    if window <= 0 :
        return _msd(coords, cells, pbc_shift, begin)
    else :
        return _msd_win(coords, cells, pbc_shift, begin, window)
        
