#!/usr/bin/python3 

import numpy as np
import dpdata,warnings

ry2ev = 13.605693009
hartree2ev = 27.211386018
bohr2ang = 0.52917721067
kbar2evperang3 = 1e3 / 1.602176621e6

length_convert = bohr2ang
energy_convert = hartree2ev
force_convert = energy_convert / length_convert

def load_key (lines, key) :
    for ii in lines :
        if key in ii :
            words = ii.split(',')
            for jj in words :
                if key in jj :
                    return jj.split('=')[1]
    return None

def load_block(lines, key, nlines) :
    for idx,ii in enumerate(lines) :
        if key in ii :
            break
    return lines[idx+1:idx+1+nlines]

def convert_celldm(ibrav, celldm) :
    if ibrav == 1 :
        return celldm[0] * np.eye(3) 
    elif ibrav == 2 :
        return celldm[0] * 0.5 * np.array([[-1,0,1], [0,1,1], [-1,1,0]])
    elif ibrav == 3 :
        return celldm[0] * 0.5 * np.array([[1,1,1], [-1,1,1], [-1,-1,1]])
    elif ibrav == -3 :
        return celldm[0] * 0.5 * np.array([[-1,1,1], [1,-1,1], [1,1,-1]])
    else :        
        warnings.warn('unsupported ibrav ' + str(ibrav) + ' if no .cel file, the cell convertion may be wrong. ')
        return np.eye(3)
        #raise RuntimeError('unsupported ibrav ' + str(ibrav))

def load_cell_parameters(lines) :
    blk = load_block(lines, 'CELL_PARAMETERS', 3)
    ret = []
    for ii in blk :
        ret.append([float(jj) for jj in ii.split()[0:3]])
    return np.array(ret)


def load_atom_names(lines, ntypes) :
    blk = load_block(lines, 'ATOMIC_SPECIES', ntypes) 
    return [ii.split()[0] for ii in blk]


def load_celldm(lines) :
    celldm = np.zeros(6)
    for ii in range(6):
        key = 'celldm(%d)' % (ii+1)
        val = load_key(lines, key)
        if val is not None :
            celldm[ii] = float(val)    
    return celldm


def load_atom_types(lines, natoms, atom_names) :
    blk = load_block(lines, 'ATOMIC_POSITIONS', natoms)
    ret = []
    for ii in blk :
        ret.append(atom_names.index(ii.split()[0]))
    return np.array(ret, dtype = int)


def load_param_file(fname) :
    with open(fname) as fp:
        lines = fp.read().split('\n')
    natoms = int(load_key(lines, 'nat'))
    ntypes = int(load_key(lines, 'ntyp'))
    atom_names = load_atom_names(lines, ntypes)
    atom_types = load_atom_types(lines, natoms, atom_names)
    atom_numbs = []
    for ii in range(ntypes) :
        atom_numbs.append(np.sum(atom_types == ii))
    ibrav = int(load_key(lines, 'ibrav'))
    celldm = load_celldm(lines)
    if ibrav == 0 :
        cell = load_cell_parameters(lines) 
    else :
        cell = convert_celldm(ibrav, celldm)
    cell = cell * length_convert
    # print(atom_names)
    # print(atom_numbs)
    # print(atom_types)
    # print(cell)
    return atom_names, atom_numbs, atom_types, cell


def _load_pos_block(fp, natoms) :
    head = fp.readline()
    if not head:
        # print('get None')
        return None, None
    else :
        ss = head.split()[0]
        blk = []
        for ii in range(natoms) :
            newline = fp.readline()
            if not newline :
                return None, None
            blk.append([float(jj) for jj in newline.split()])
        return blk, ss


def load_data(fname, 
              natoms, 
              begin = 0, 
              step = 1,
              convert = 1.) :
    coords = []
    steps = []
    cc = 0
    with open(fname) as fp:
        while True:
            blk, ss = _load_pos_block(fp, natoms)
            if blk == None :
                break
            else :
                if cc >= begin and (cc - begin) % step == 0 :
                    coords.append(blk)
                    steps.append(ss)
            cc += 1
    coords= convert * np.array(coords)
    return coords, steps


# def load_pos(fname, natoms) :
#     coords = []
#     with open(fname) as fp:
#         while True:
#             blk = _load_pos_block(fp, natoms)
#             # print(blk)
#             if blk == None :
#                 break
#             else :
#                 coords.append(blk)
#     coords= length_convert * np.array(coords)
#     return coords


def load_energy(fname, begin = 0, step = 1) :
    data = np.loadtxt(fname)
    steps = []
    for ii in data[begin::step,0]:
        steps.append('%d'%ii)
    with open(fname) as fp:
        while True:
            line = fp.readline()
            if not line :
                return None
            if line.split()[0][0] != '#':
                nw = len(line.split())
                break
    data = np.reshape(data, [-1, nw])
    return energy_convert * data[begin::step,5], steps


# def load_force(fname, natoms) :
#     coords = []
#     with open(fname) as fp:
#         while True:
#             blk = _load_pos_block(fp, natoms)
#             # print(blk)
#             if blk == None :
#                 break
#             else :
#                 coords.append(blk)
#     coords= force_convert * np.array(coords)
#     return coords


def to_system_data(input_name, prefix, begin = 0, step = 1) :
    data = {}
    data['atom_names'], \
        data['atom_numbs'], \
        data['atom_types'], \
        cell \
        = load_param_file(input_name)
    data['coords'], csteps\
        = load_data(prefix + '.pos', 
                    np.sum(data['atom_numbs']), 
                    begin = begin, 
                    step = step, 
                    convert = length_convert)
    data['orig'] = np.zeros(3)
    try :
        data['cells'], tmp_steps \
            = load_data(prefix + '.cel', 
                        3, 
                        begin = begin, 
                        step = step, 
                        convert = length_convert)
        assert(csteps == tmp_steps), "the step key between files are not consistent"
    except FileNotFoundError :
        data['cells'] = np.tile(cell, (data['coords'].shape[0], 1, 1))
    return data, csteps


def to_system_label(input_name, prefix, begin = 0, step = 1) :
    atom_names, atom_numbs, atom_types, cell = load_param_file(input_name)
    energy, esteps = load_energy(prefix + '.evp', 
                                 begin = begin, 
                                 step = step)
    force, fsteps = load_data(prefix + '.for', 
                              np.sum(atom_numbs), 
                              begin = begin,
                              step = step, 
                              convert = force_convert)
    assert(esteps == fsteps), "the step key between files are not consistent "
    return energy, force, esteps


if __name__ == '__main__':
    prefix='nacl'
    atom_names, atom_numbs, atom_types, cell = load_param_file(prefix+'.in')
    coords = load_data(prefix+'.pos', np.sum(atom_numbs))
    cells = load_data(prefix+'.cel', 3)
    print(atom_names)
    print(atom_numbs)
    print(atom_types)
    print(cells)
    print(coords.shape)
    print(cells.shape)
