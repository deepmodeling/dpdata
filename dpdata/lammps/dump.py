#!/usr/bin/env python3

import os, sys
import numpy as np
lib_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(lib_path)
import lmp
import warnings
class UnwrapWarning(UserWarning):
    pass
warnings.simplefilter('once', UnwrapWarning)


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

def get_atype(lines, type_idx_zero = False) :
    blk, head = _get_block(lines, 'ATOMS')
    keys = head.split()
    id_idx = keys.index('id') - 2
    tidx = keys.index('type') - 2
    atype = []
    for ii in blk :
        atype.append([int(ii.split()[tidx]), int(ii.split()[id_idx])])
    # sort with type id
    atype.sort()
    atype = np.array(atype, dtype = int)    
    atype = atype[:, ::-1]
    if type_idx_zero :
        return atype[:,1] - 1
    else :
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

def get_coordtype_and_scalefactor(keys):
    # 4 types in total,with different scaling factor
    key_pc=['x','y','z'] # plain cartesian, sf = 1
    key_uc=['xu','yu','zu'] # unwraped cartesian, sf = 1
    key_s=['xs','ys','zs'] # scaled by lattice parameter, sf = lattice parameter
    key_su = ['xsu','ysu','zsu'] #scaled and unfolded,sf = lattice parameter
    lmp_coor_type = [key_pc,key_uc,key_s,key_su]
    sf = [0,0,1,1]
    uw = [0,1,0,1]  # unwraped or not
    for k in range(4):
        if all(i in keys for i in lmp_coor_type[k]):
            return lmp_coor_type[k], sf[k], uw[k]

def safe_get_posi(lines,cell,orig=np.zeros(3), unwrap=False) :
    blk, head = _get_block(lines, 'ATOMS')
    keys = head.split()
    coord_tp_and_sf = get_coordtype_and_scalefactor(keys)
    assert coord_tp_and_sf is not None, 'Dump file does not contain atomic coordinates!'
    coordtype, sf, uw = coord_tp_and_sf
    id_idx = keys.index('id') - 2
    tidx = keys.index('type') - 2
    xidx = keys.index(coordtype[0])-2
    yidx = keys.index(coordtype[1])-2
    zidx = keys.index(coordtype[2])-2
    sel = (xidx, yidx, zidx)
    posis = []
    for ii in blk :
        words = ii.split()
        posis.append([float(words[tidx]), float(words[id_idx]), float(words[xidx]), float(words[yidx]), float(words[zidx])])
    posis.sort()
    posis = np.array(posis)[:,2:5]
    if not sf:
        posis = (posis - orig) @ np.linalg.inv(cell)  # Convert to scaled coordinates for unscaled coordinates
    if uw and unwrap:
        return posis @ cell  # convert scaled coordinates back to Cartesien coordinates unwrap at the periodic boundaries
    else:
        if uw and not unwrap:
            warnings.warn(message='Your dump file contains unwrapped coordinates, but you did not specify unwrapping (unwrap = True). The default is wrapping at periodic boundaries (unwrap = False).\n',category=UnwrapWarning)
        return (posis % 1) @ cell  # Convert scaled coordinates back to Cartesien coordinates with wraping at periodic boundary conditions

def get_dumpbox(lines) :
    blk, h = _get_block(lines, 'BOX BOUNDS')
    bounds = np.zeros([3,2])
    tilt = np.zeros([3])
    load_tilt = 'xy xz yz' in h
    for dd in range(3) :
        info = [float(jj) for jj in blk[dd].split()]
        bounds[dd][0] = info[0]
        bounds[dd][1] = info[1]
        if load_tilt :
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


def load_file(fname, begin = 0, step = 1) :
    lines = []
    buff = []
    cc = -1
    with open(fname) as fp:
        while True:
            line = fp.readline().rstrip('\n')
            if not line :
                if cc >= begin and (cc - begin) % step == 0 :
                    lines += buff
                    buff = []
                cc += 1
                return lines
            if 'ITEM: TIMESTEP' in line :
                if cc >= begin and (cc - begin) % step == 0 :
                    lines += buff
                    buff = []
                cc += 1
            if cc >= begin and (cc - begin) % step == 0 :
                buff.append(line)


def system_data(lines, type_map = None, type_idx_zero = True, unwrap=False) :
    array_lines = split_traj(lines)
    lines = array_lines[0]
    system = {}
    system['atom_numbs'] = get_natoms_vec(lines)
    system['atom_names'] = []
    if type_map == None :
        for ii in range(len(system['atom_numbs'])) :
            system['atom_names'].append('TYPE_%d' % ii)
    else :
        assert(len(type_map) >= len(system['atom_numbs']))
        for ii in range(len(system['atom_numbs'])) :
            system['atom_names'].append(type_map[ii])
    bounds, tilt = get_dumpbox(lines)
    orig, cell = dumpbox2box(bounds, tilt)
    system['orig'] = np.array(orig) - np.array(orig)
    system['cells'] = [np.array(cell)]
    natoms = sum(system['atom_numbs'])
    system['atom_types'] = get_atype(lines, type_idx_zero = type_idx_zero)
    system['coords'] = [safe_get_posi(lines, cell, np.array(orig), unwrap)]
    for ii in range(1, len(array_lines)) :
        bounds, tilt = get_dumpbox(array_lines[ii])
        orig, cell = dumpbox2box(bounds, tilt)
        system['cells'].append(cell)
        system['coords'].append(safe_get_posi(array_lines[ii], cell, np.array(orig), unwrap))
    system['cells'] = np.array(system['cells'])
    system['coords'] = np.array(system['coords'])
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
    # fname = 'dump.hti'
    # lines = open(fname).read().split('\n')
    # # print(get_natoms(lines))
    # # print(get_natomtypes(lines))
    # # print(get_natoms_vec(lines))
    # posi = get_posi(lines)
    # dbox1, tilt1 = box2dumpbox(orig, box)
    # print(dbox - dbox1)
    # print(tilt - tilt1)
    # print(orig)
    # print(box)
    # np.savetxt('tmp.out', posi - orig, fmt='%.6f')
    # print(system_data(lines))
    lines = load_file('conf_unfold.dump', begin = 0, step = 1)
    al = split_traj(lines)
    s = system_data(lines,['O','H'])
    #l = np.linalg.norm(s['cells'][1],axis=1)
    #p = s['coords'][0] + l
    #np.savetxt('p',p,fmt='%1.10f')
