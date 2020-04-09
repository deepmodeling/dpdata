#!/usr/bin/env python3

import numpy as np

nm2ang = 10.

def _get_line(line):
    atom_name = line[10:15].split()[0]
    atom_idx = int(line[15:20].split()[0])
    posis = [float(line[ii:ii+8]) for ii in range(20,44,8)]
    posis = np.array(posis) * nm2ang
    return atom_name, atom_idx, posis

def _get_cell(line):
    cell = np.zeros([3,3])
    lengths = [float(ii) for ii in line.split()]
    if len(lengths) >= 3:
        for dd in range(3):
            cell[dd][dd] = lengths[dd]
    else:
        raise RuntimeError('wrong box format: ', line)
    if len(lengths) == 9:
        cell[0][1] = lengths[3]
        cell[0][2] = lengths[4]
        cell[1][0] = lengths[5]
        cell[1][2] = lengths[6]
        cell[2][0] = lengths[7]
        cell[2][1] = lengths[8]
    cell = cell * nm2ang
    return cell

def file_to_system_data(fname):
    names = []
    idxs = []
    posis = []
    with open(fname) as fp:
        fp.readline()
        natoms = int(fp.readline())
        for ii in range(natoms):
            n, i, p = _get_line(fp.readline())
            names.append(n)
            idxs.append(i)
            posis.append(p)
        cell = _get_cell(fp.readline())
    posis = np.array(posis)    
    system = {}
    system['orig'] = np.array([0, 0, 0])    
    system['atom_names'] = list(set(names))
    system['atom_names'].sort()
    system['atom_numbs'] = [names.count(ii) for ii in system['atom_names']]
    system['atom_types'] = [system['atom_names'].index(ii) for ii in names]
    system['atom_types'] = np.array(system['atom_types'], dtype = int)
    system['coords'] = np.array([posis])
    system['cells'] = np.array([cell])
    return system
