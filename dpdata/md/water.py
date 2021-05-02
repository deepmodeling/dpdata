import numpy as np
from .pbc import posi_diff
from .pbc import posi_shift

def compute_bonds (box, 
                   posis,
                   atype, 
                   oh_sel = [0,1],
                   max_roh = 1.3, 
                   uniq_hbond = True):
    try :
        import ase
        import ase.neighborlist    
        # nlist implemented by ase
        return compute_bonds_ase(box, posis, atype, oh_sel, max_roh, uniq_hbond)        
    except ImportError:
        # nlist naivly implemented , scales as O(N^2)
        return compute_bonds_naive(box, posis, atype, oh_sel, max_roh, uniq_hbond)


def compute_bonds_ase(box, 
                      posis,
                      atype, 
                      oh_sel = [0,1],
                      max_roh = 1.3, 
                      uniq_hbond = True):
    natoms = len(posis)
    from ase import Atoms
    import ase.neighborlist    
    atoms = Atoms(positions=posis, cell=box, pbc=[1,1,1])
    nlist = ase.neighborlist.NeighborList(max_roh, self_interaction=False, bothways=True, primitive=ase.neighborlist.NewPrimitiveNeighborList)
    nlist.update(atoms)
    bonds = []
    o_type = oh_sel[0]
    h_type = oh_sel[1]
    for ii in range(natoms) :
        bonds.append([])
    for ii in range(natoms) :
        if atype[ii] == o_type :
            nn, ss = nlist.get_neighbors(ii)
            for jj in nn:
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
                      

def compute_bonds_naive(box,
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



def pbc_coords(cells, 
               coords, 
               atom_types, 
               oh_sel = [0, 1],
               max_roh = 1.3):
    bonds = compute_bonds(cells, coords, atom_types, oh_sel = oh_sel, max_roh = max_roh, uniq_hbond = True)

    new_coords = np.copy(coords)
    natoms = len(atom_types)
    o_type = oh_sel[0]
    h_type = oh_sel[1]
    for ii in range(natoms):
        if atom_types[ii] == o_type:
            assert(len(bonds[ii]) == 2), 'O has more than 2 bonded Hs, stop'
            for jj in bonds[ii]:
                assert(atom_types[jj] == h_type), 'The atom bonded to O is not H, stop'
                shift = posi_shift(cells, coords[jj], coords[ii])
                new_coords[jj] = coords[jj] + np.matmul(shift, cells)
    return new_coords

