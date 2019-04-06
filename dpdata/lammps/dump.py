#!/usr/bin/env python3

import os, sys
import numpy as np
lib_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(lib_path)
import lmp

def _get_block (lines, key) :
    for idx in range(len(lines)) :
        if ('ITEM: ' + key) in lines[idx] :
            break
    idx_s = idx + 1
    for idx in range(idx_s, len(lines)) :
        if ('ITEM: ') in lines[idx] :
            break
    idx_e = idx
    if idx_e == len(lines)-1 :
        idx_e += 1
    return lines[idx_s:idx_e], lines[idx_s-1]

def get_atype(lines) :
    blk, head = _get_block(lines, 'ATOMS')
    keys = head.split()
    id_idx = keys.index('id') - 2
    tidx = keys.index('type') - 2
    atype = []
    for ii in blk :
        atype.append([int(ii.split()[id_idx]), int(ii.split()[tidx])])
    atype.sort()
    atype = np.array(atype, dtype = int)
    return atype[:,1]

def get_natoms(lines) :
    blk, head = _get_block(lines, 'NUMBER OF ATOMS')
    return int(blk[0])

def get_natomtypes(lines) :
    atype = get_atype(lines)
    return max(atype)

def get_natoms_vec(lines) :
    atype = get_atype(lines)
    natoms_vec = []
    natomtypes = get_natomtypes(lines)
    for ii in range(natomtypes) :
        natoms_vec.append(sum(atype == ii+1))
    assert (sum(natoms_vec) == get_natoms(lines))
    return natoms_vec

def get_posi(lines) :
    blk, head = _get_block(lines, 'ATOMS')
    keys = head.split()
    id_idx = keys.index('id') - 2
    xidx = keys.index('x') - 2
    yidx = keys.index('y') - 2
    zidx = keys.index('z') - 2
    sel = (xidx, yidx, zidx)
    posis = []
    for ii in blk :
        words = ii.split()
        posis.append([float(words[id_idx]), float(words[xidx]), float(words[yidx]), float(words[zidx])])
    posis.sort()
    posis = np.array(posis)     
    return posis[:,1:4]

def get_dumpbox(lines) :
    blk, h = _get_block(lines, 'BOX BOUNDS')
    bounds = np.zeros([3,2])
    tilt = np.zeros([3])
    for dd in range(3) :
        info = [float(jj) for jj in blk[dd].split()]
        bounds[dd][0] = info[0]
        bounds[dd][1] = info[1]
        tilt[dd] = info[2]
    return bounds, tilt

def dumpbox2box(bounds, tilt) :
    xy = tilt[0]
    xz = tilt[1]
    yz = tilt[2]
    xlo = bounds[0][0] - min(0.0,xy,xz,xy+xz)
    xhi = bounds[0][1] - max(0.0,xy,xz,xy+xz)
    ylo = bounds[1][0] - min(0.0,yz)
    yhi = bounds[1][1] - max(0.0,yz)
    zlo = bounds[2][0]
    zhi = bounds[2][1]
    info = [[xlo, xhi], [ylo, yhi], [zlo, zhi]]
    return lmp.lmpbox2box(info, tilt)

def box2dumpbox(orig, box) :
    lohi, tilt = lmp.box2lmpbox(orig, box)
    xy = tilt[0]
    xz = tilt[1]
    yz = tilt[2]
    bounds = np.zeros([3,2])
    bounds[0][0] = lohi[0][0] + min(0.0,xy,xz,xy+xz)
    bounds[0][1] = lohi[0][1] + max(0.0,xy,xz,xy+xz)
    bounds[1][0] = lohi[1][0] + min(0.0,yz)
    bounds[1][1] = lohi[1][1] + max(0.0,yz)
    bounds[2][0] = lohi[2][0]
    bounds[2][1] = lohi[2][1]
    return bounds, tilt

def system_data(lines) :
    system = {}
    system['atom_numbs'] = get_natoms_vec(lines)
    system['atom_names'] = []
    for ii in range(len(system['atom_numbs'])) :
        system['atom_names'].append('TYPE_%d' % ii)
    bounds, tilt = get_dumpbox(lines)
    orig, cell = dumpbox2box(bounds, tilt)
    system['orig'] = np.array(orig)
    system['cell'] = np.array(cell)
    natoms = sum(system['atom_numbs'])
    system['atom_types'] = get_atype(lines)
    system['coordinates'] = get_posi(lines)
    return system


def split_traj(dump_lines) :
    marks = []
    for idx,ii in enumerate(dump_lines) :
        if 'ITEM: TIMESTEP' in ii :
            marks.append(idx) 
    if len(marks) == 0 :
        return None
    elif len(marks) == 1 :
        return [dump_lines]
    else :
        block_size = marks[1] - marks[0]
        ret = []
        for ii in marks :
            ret.append(dump_lines[ii:ii+block_size])
        # for ii in range(len(marks)-1): 
        #     assert(marks[ii+1] - marks[ii] == block_size)
        return ret    
    return None


if __name__ == '__main__' :
    fname = 'dump.hti'
    lines = open(fname).read().split('\n')
    # print(get_natoms(lines))
    # print(get_natomtypes(lines))
    # print(get_natoms_vec(lines))
    posi = get_posi(lines)
    dbox, tilt = get_dumpbox(lines)
    orig, box = dumpbox2box(dbox, tilt)
    dbox1, tilt1 = box2dumpbox(orig, box)
    print(dbox - dbox1)
    print(tilt - tilt1)
    print(orig)
    print(box)
    np.savetxt('tmp.out', posi - orig, fmt='%.6f')
    print(system_data(lines))
