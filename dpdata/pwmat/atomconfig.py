#!/usr/bin/python3 

import numpy as np

def _to_system_data_lower(lines) :
    system = {}
    natoms = int(lines[0].split()[0])
    cell = []
    for idx, ii in enumerate(lines):
        if 'lattice' in ii or 'Lattice' in ii or 'LATTICE' in ii:
            for kk in range(idx+1,idx+1+3):
                vector=[float(jj) for jj in lines[kk].split()[0:3]]
                cell.append(vector)
    system['cells'] = [np.array(cell)]
    coord = []
    atomic_number = []
    atom_numbs = []
    for idx, ii in enumerate(lines):
        if 'Position' in ii or 'POSITION' in ii or 'position' in ii:
            for kk in range(idx+1,idx+1+natoms):
                min = kk
                for jj in range(kk+1,idx+1+natoms):
                    if int(lines[jj].split()[0]) < int(lines[min].split()[0]):
                        min = jj
                        lines[min], lines[kk] = lines[kk],lines[min]
            for gg in range(idx+1,idx+1+natoms):
                tmpv = [float(jj) for jj in lines[gg].split()[1:4]]
                tmpv = np.matmul(np.array(tmpv), system['cells'][0])
                coord.append(tmpv)
                tmpn = int(lines[gg].split()[0])
                atomic_number.append(tmpn)
    for ii in np.unique(sorted(atomic_number)) :
        atom_numbs.append(atomic_number.count(ii))
    system['atom_numbs'] = [int(ii) for ii in atom_numbs]
    system['coords'] = [np.array(coord)]
    system['orig'] = np.zeros(3)
    atom_types = []
    for idx,ii in enumerate(system['atom_numbs']) :
        for jj in range(ii) :
            atom_types.append(idx)
    system['atom_types'] = np.array(atom_types, dtype = int)
    ELEMENTS=['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', \
            'Sc', 'Ti', 'V', 'Cr','Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', \
            'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',\
            'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', \
            'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', \
            'Md', 'No', 'Lr']

    system['atom_names'] = [ELEMENTS[ii-1] for ii in np.unique(sorted(atomic_number))]
    return system


def to_system_data(lines) :
    return _to_system_data_lower(lines)


def from_system_data(system, f_idx = 0, skip_zeros = True) :
    ret = ''
    natoms = sum(system['atom_numbs'])
    ret += '%d' % natoms
    ret += '\n'
    ret += 'LATTICE'
    ret += '\n'
    for ii in system['cells'][f_idx] :
        for jj in ii :
            ret += '%.16e ' % jj
        ret += '\n'
    ret += 'POSITION'
    ret += '\n'
    atom_numbs = system['atom_numbs']
    atom_names = system['atom_names']
    atype = system['atom_types']
    posis = system['coords'][f_idx]    
    # atype_idx = [[idx,tt] for idx,tt in enumerate(atype)]
    # sort_idx = np.argsort(atype, kind = 'mergesort')
    sort_idx = np.lexsort((np.arange(len(atype)), atype))
    atype = atype[sort_idx]
    posis = posis[sort_idx]
    symbal = []
    for ii, jj in zip(atom_numbs, atom_names):
        for kk in range(ii):
            symbal.append(jj) 
    ELEMENTS=['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', \
            'Sc', 'Ti', 'V', 'Cr','Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', \
            'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',\
            'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', \
            'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', \
            'Md', 'No', 'Lr']
    atomic_numbers = []
    for ii in symbal:
        atomic_numbers.append(ELEMENTS.index(ii)+1)
    posi_list = []
    for jj, ii in zip(atomic_numbers,posis) :
        ii = np.matmul(ii, np.linalg.inv(system['cells'][0]))
        posi_list.append('%d %15.10f %15.10f %15.10f 1 1 1' % \
                         (jj, ii[0], ii[1], ii[2])
        )
    for kk in range(len(posi_list)):
        min = kk
        for jj in range(kk,len(posi_list)):
            if int(posi_list[jj].split()[0]) < int(posi_list[min].split()[0]):
                min = jj
                posi_list[min], posi_list[kk] = posi_list[kk],posi_list[min]
    posi_list.append('')
    ret += '\n'.join(posi_list)
    return ret
