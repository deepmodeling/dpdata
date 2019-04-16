import numpy as np
from .pbc import posi_diff

def compute_bonds(box,
                  posis, 
                  atype, 
                  oh_sel = [0,1],
                  max_roh = 1.3, 
                  uniq_hbond = True):
    natoms = len(posis)
    bonds = []
    o_type = oh_sel[0]
    h_type = oh_sel[1]
    for ii in range(natoms) :
        bonds.append([])
    for ii in range(natoms) :
        if atype[ii] == o_type :
            for jj in range(natoms) :
                if atype[jj] == h_type :
                    dr = posi_diff(box, posis[ii], posis[jj])
                    if np.linalg.norm(dr) < max_roh :
                        bonds[ii].append(jj)
                        bonds[jj].append(ii)
    if uniq_hbond :
        for jj in range(natoms) :
            if atype[jj] == h_type :
                if len(bonds[jj]) > 1 :
                    orig_bonds = bonds[jj]
                    min_bd = 1e10
                    min_idx = -1
                    for ii in bonds[jj] :
                        dr = posi_diff(box, posis[ii], posis[jj])
                        drr = np.linalg.norm(dr)
                        # print(ii,jj, posis[ii], posis[jj], drr)
                        if drr < min_bd :
                            min_idx = ii
                            min_bd = drr
                    bonds[jj] = [min_idx]
                    orig_bonds.remove(min_idx)
                    # print(min_idx, orig_bonds, bonds[jj])
                    for ii in orig_bonds :
                        bonds[ii].remove(jj)
    return bonds


# def ions_count (atype, 
#                 bonds, 
#                 oh_sel = [0, 1]) :
#     no = 0
#     noh = 0
#     noh2 = 0
#     noh3 = 0
#     nh = 0
#     natoms = len(atype)
#     o_type = oh_sel[0]
#     h_type = oh_sel[1]
#     for ii in range(natoms) :
#         if atype[ii] == o_type :
#             if len(bonds[ii] ) != 2 :
#                 # print('# ', ii, bonds[ii])
#                 pass
#             if len(bonds[ii] ) == 0 :
#                 no += 1
#             elif len(bonds[ii] ) == 1 :
#                 noh += 1
#             elif len(bonds[ii] ) == 2 :
#                 noh2 += 1
#             elif len(bonds[ii] ) == 3 :
#                 noh3 += 1
#             else :
#                 raise RuntimeError("unknow case: numb of H bonded to O > 3")
#     for ii in range(natoms) :
#         if atype[ii] == h_type :
#             if len(bonds[ii] ) != 1 :
#                 print('# ', ii, bonds[ii])
#             if len(bonds[ii] ) == 0 :
#                 nh += 1
#             elif len(bonds[ii] ) == 1 :
#                 pass
#             else :
#                 raise RuntimeError("unknow case: numb of O bonded to H > 1")
#     return no, noh, noh2, noh3, nh

def find_ions (atype, 
               bonds, 
               oh_sel = [0, 1], 
               ret_h2o = True) :
    no = []
    noh = []
    noh2 = []
    noh3 = []
    nh = []
    natoms = len(atype)
    o_type = oh_sel[0]
    h_type = oh_sel[1]
    for ii in range(natoms) :
        if atype[ii] == o_type :
            if len(bonds[ii] ) == 0 :
                no.append(ii)
            elif len(bonds[ii] ) == 1 :
                noh.append(ii)
            elif len(bonds[ii] ) == 2 :
                if ret_h2o :
                    noh2.append(ii)
            elif len(bonds[ii] ) == 3 :
                noh3.append(ii)
            else :
                raise RuntimeError("unknow case: numb of H bonded to O > 3")
    for ii in range(natoms) :
        if atype[ii] == h_type :
            if len(bonds[ii] ) == 0 :
                nh.append(ii)
            elif len(bonds[ii] ) == 1 :
                pass
            else :
                raise RuntimeError("unknow case: numb of O bonded to H > 1")
    return no, noh, noh2, noh3, nh
