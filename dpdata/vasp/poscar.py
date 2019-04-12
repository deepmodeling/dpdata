#!/usr/bin/python3 

import numpy as np

def _to_system_data_lower(lines, cartesian = True) :
    '''
    Treat as cartesian poscar
    '''
    system = {}
    system['atom_names'] = [str(ii) for ii in lines[5].split()]
    system['atom_numbs'] = [int(ii) for ii in lines[6].split()]
    scale = float(lines[1])
    cell = []
    for ii in range(2,5) :
        boxv = [float(jj) for jj in lines[ii].split()]
        boxv = np.array(boxv) * scale
        cell.append(boxv)
    system['cells'] = [np.array(cell)]
    natoms = sum(system['atom_numbs'])
    coord = []
    for ii in range(8, 8+natoms) :
        tmpv = [float(jj) for jj in lines[ii].split()[:3]]
        if cartesian :
            tmpv = np.array(tmpv) * scale
        else :
            tmpv = np.matmul(np.array(tmpv), system['cells'][0])
        coord.append(tmpv)
    system['coords'] = [np.array(coord)]
    system['orig'] = np.zeros(3)
    atom_types = []
    for idx,ii in enumerate(system['atom_numbs']) :
        for jj in range(ii) :
            atom_types.append(idx)
    system['atom_types'] = np.array(atom_types, dtype = int)
    system['cells'] = np.array(system['cells'])
    system['coords'] = np.array(system['coords'])
    return system


def to_system_data(lines) :
    # remove the line that has 'selective dynamics'
    if lines[7][0] == 'S' or lines[7][0] == 's' :
        lines.pop(7)
    is_cartesian = (lines[7][0] in ['C', 'c', 'K', 'k'])
    if not is_cartesian :
        if not (lines[7][0] in ['d', 'D']) :
            raise RuntimeError('seem not to be a valid POSCAR of vasp 5.x, may be a POSCAR of vasp 4.x?')
    return _to_system_data_lower(lines, is_cartesian)


def from_system_data(system, f_idx = 0) :
    ret = ''
    for ii,name in zip(system['atom_numbs'], system['atom_names']) :
        ret += '%s%d ' % (name, ii)
    ret += '\n'
    ret += '1.0\n'
    for ii in system['cells'][f_idx] :
        for jj in ii :
            ret += '%.16e ' % jj
        ret += '\n'
    for ii in system['atom_names'] :
        ret += '%s ' % ii
    ret += '\n'
    for ii in system['atom_numbs'] :
        ret += '%d ' % ii
    ret += '\n'
    ret += 'cartesian\n'
    atype = system['atom_types']
    posis = system['coords'][f_idx]    
    # atype_idx = [[idx,tt] for idx,tt in enumerate(atype)]
    # sort_idx = np.argsort(atype, kind = 'mergesort')
    sort_idx = np.lexsort((np.arange(len(atype)), atype))
    atype = atype[sort_idx]
    posis = posis[sort_idx]
    posi_list = []
    for ii in posis :
        posi_list.append('%15.10f %15.10f %15.10f' % \
                         (ii[0], ii[1], ii[2])
        )
    posi_list.append('')
    ret += '\n'.join(posi_list)
    return ret
