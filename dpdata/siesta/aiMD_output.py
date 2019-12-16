# !/usr/bin/python3

import numpy as np

ev2ev = 1
ang2ang = 1

#############################read output#####################################
def get_single_line_tail(fin, keyword, num=1):
    file = open(fin, 'r')
    part_res = []
    for value in file:
        if keyword in value:
            temp = len(value.split()) - num
            part_res.append(float(value.split()[temp]))

    file.close()
    return part_res

## atomnum: number of atoms,  row numbers
## begin_column: begin column num
## read_column_num: read column num
## column_num: the column number in nxet reading line
def extract_keyword(fout, keyword, down_line_num, begin_column, read_column_num, is_repeated_read, column_num):
    file = open(fout, 'r')
    ret = []
    part_ret = []
    flag = 0
    idx = 0
    extr_frame = 0
    length = obtain_nframe(fout)
    # for (num,value) in enumerate(file):
    for value in file:
        if keyword in value:
            flag = 1
            continue
        if flag == 1:
            if idx < down_line_num:
                idx += 1
            else:
                flag = 0
                part_ret.append(np.array(ret))
                ret = []
                extr_frame += 1
                if extr_frame == length:
                    file.close()
                    return part_ret
                ## is_repeated_read = 0: only read 1 time for SCF
                ## is_repeated_read = 1:  read all for aiMD --> get all frames
                if is_repeated_read:
                    idx = 0
                continue

            for i in range(begin_column, read_column_num):
                if len(value.split()) == column_num:
                    if not value.split()[i].isalpha():
                        ret.append(float(value.strip().split()[i]))
                    else:
                        ret.append(value.strip().split()[i])
            continue
    file.close()
    return part_ret

def obtain_nframe(fname):
    fp = open(fname, 'r')
    flag = False
    idx = 0
    temp = 0
    for ii in fp:
        if 'siesta: Stress tensor (static) (eV/Ang**3):' in ii:
            flag = True
            continue
        if flag:
            if not 'siesta: Pressure (static):' in ii:
                if len(ii.split()) == 3:
                    temp += 1
                    if temp == 3:
                        idx += 1
                        # print(idx)
                        flag = False
                        temp = 0
    fp.close()
    return idx

def get_atom_types(fout, atomnums):
    covert_type = extract_keyword(fout, 'outcoor: Atomic coordinates (Ang):', atomnums, 3, 4, 0, 6)[0]
    atomtype = []
    # print(covert_type)
    for i in range(0, len(covert_type)):
        atomtype.append(int(covert_type[i]) - 1)
    return atomtype

def get_atom_name(fout):
    file = open(fout, 'r')
    ret = []
    for value in file:
        if 'Species number:' in value:
            for j in range(len(value.split())):
                if value.split()[j] == 'Label:':
                    ret.append(value.split()[j+1])
                    break              
    file.close()
    return ret

def get_atom_numbs(atomtypes):
    atom_numbs = []
    for i in set(atomtypes):
        atom_numbs.append(atomtypes.count(i))
    return atom_numbs

def get_virial(fout, cell):
    viri = extract_keyword(fout, 'siesta: Stress tensor (static) (eV/Ang**3):', 3, 0, 3, 1, 3)
    vols = []
    length = obtain_nframe(fout)
    for ii in range(length):
        vols.append(np.linalg.det(cell[ii].reshape([3, 3])))
        for jj in range(len(viri[ii])):
            ## siesta: 1eV/A^3= 1.60217*10^11 Pa ,  ---> qe: kBar=10^8Pa
            # ii *= vols[idx] * 1e3 / 1.602176621e6 * (1.602176621e3)
            viri[ii][jj] *= vols[ii]
    return viri

def covert_dimension(arr, num):
    arr = np.array(arr)
    frames = len(arr)
    ret = np.zeros((frames, num, 3))
    for i in range(frames):
        ret[i] = arr[i].reshape(num, 3)
    return ret

def get_aiMD_frame(fname):
    NumberOfSpecies = int(get_single_line_tail(fname, 'redata: Number of Atomic Species')[0])
    atom_names = get_atom_name(fname)
    tot_natoms = int(get_single_line_tail(fname, 'Number of atoms', 3)[0])

    atom_types = get_atom_types(fname, tot_natoms)
    atom_numbs = get_atom_numbs(atom_types)
    assert (max(atom_types) + 1 == NumberOfSpecies)

    cell = extract_keyword(fname, 'outcell: Unit cell vectors (Ang):', 3, 0, 3, 1, 3)
    coord = extract_keyword(fname, 'outcoor: Atomic coordinates (Ang):', tot_natoms, 0, 3, 1, 6)
    energy = get_single_line_tail(fname, 'siesta: E_KS(eV) =')
    force = extract_keyword(fname, 'siesta: Atomic forces (eV/Ang):', tot_natoms, 1, 4, 1, 4)
    virial = get_virial(fname, cell)

    cells = covert_dimension(np.array(cell), 3)
    coords = covert_dimension(np.array(coord), tot_natoms)
    forces = covert_dimension(np.array(force), tot_natoms)
    virials = covert_dimension(np.array(virial), 3)
    return atom_names, atom_numbs, np.array(atom_types), cells, coords, np.array(energy), forces, virials
