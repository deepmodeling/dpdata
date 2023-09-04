#!/usr/bin/python3

import numpy as np

ev2ev = 1
ang2ang = 1


#############################read output#####################################
def get_single_line_tail(fin, keyword, num=1):
    file = open(fin)
    res = []
    for value in file:
        if keyword in value:
            temp = len(value.split()) - num
            res.append(float(value.split()[temp]))
            file.close()
            return res
    return res


## atomnum: number of atoms,  row numbers
## begin_column: begin column num
## column_num: read column num
def extract_keyword(fout, keyword, down_line_num, begin_column, column_num):
    file = open(fout)
    ret = []
    flag = 0
    idx = 0
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
                continue
            if len(value.split()) >= column_num:
                for i in range(begin_column, column_num):
                    if not value.split()[i].isalpha():
                        ret.append(float(value.strip().split()[i]))
                    else:
                        ret.append(value.strip().split()[i])
            ## compatible siesta-4.0.2 and siesta-4.1-b4
            else:
                flag = 0
                idx = 0
    file.close()
    return ret


def get_atom_types(fout, atomnums):
    covert_type = extract_keyword(
        fout, "outcoor: Atomic coordinates (Ang):", atomnums, 3, 4
    )
    atomtype = []
    for i in range(0, len(covert_type)):
        atomtype.append(int(covert_type[i]) - 1)
    return atomtype


def get_atom_name(fout):
    file = open(fout)
    ret = []
    for value in file:
        if "Species number:" in value:
            for j in range(len(value.split())):
                if value.split()[j] == "Label:":
                    ret.append(value.split()[j + 1])
                    break
    file.close()
    return ret


def get_atom_numbs(atomtypes):
    atom_numbs = []
    for i in set(atomtypes):
        atom_numbs.append(atomtypes.count(i))
    return atom_numbs


def get_virial(fout, cells):
    vols = []
    for ii in cells:
        ### calucate vol
        vols.append(np.linalg.det(ii.reshape([3, 3])))
    ret = extract_keyword(fout, "siesta: Stress tensor (static) (eV/Ang**3):", 3, 1, 4)
    ret = np.array([ret])
    for idx, ii in enumerate(ret):
        ## siesta: 1eV/A^3= 1.60217*10^11 Pa ,  ---> qe: kBar=10^8Pa
        # ii *= vols[idx] * 1e3 / 1.602176621e6 * (1.602176621e3)
        ii *= vols[idx]
    return ret


def obtain_frame(fname):
    NumberOfSpecies = int(
        get_single_line_tail(fname, "redata: Number of Atomic Species")[0]
    )
    atom_names = get_atom_name(fname)
    tot_natoms = int(get_single_line_tail(fname, "Number of atoms", 3)[0])
    atom_types = get_atom_types(fname, tot_natoms)
    atom_numbs = get_atom_numbs(atom_types)
    assert max(atom_types) + 1 == NumberOfSpecies
    cell = extract_keyword(fname, "outcell: Unit cell vectors (Ang):", 3, 0, 3)
    coord = extract_keyword(
        fname, "outcoor: Atomic coordinates (Ang):", tot_natoms, 0, 3
    )
    energy = get_single_line_tail(fname, "siesta: E_KS(eV) =")
    force = extract_keyword(fname, "siesta: Atomic forces (eV/Ang):", tot_natoms, 1, 4)
    virial = get_virial(fname, np.array([cell]))

    cell = np.array(cell).reshape(3, 3)
    coord = np.array(coord).reshape(tot_natoms, 3)
    force = np.array(force).reshape(tot_natoms, 3)
    virial = np.array(virial).reshape(3, 3)

    # data = {}
    # data['orig'] = np.array([0, 0, 0])
    # data['atom_names'] = atom_names
    # data['atom_numbs'] = atom_numbs
    # data['atom_types'] = np.array(atom_types)
    # data['cells'] = np.array([cell])
    # data['coords'] = np.array([coord])
    # data['energies'] = np.array([energy])
    # data['forces'] = np.array([force])
    # data['virials'] = virial
    # return data
    return (
        atom_names,
        atom_numbs,
        np.array(atom_types),
        np.array([cell]),
        np.array([coord]),
        np.array(energy),
        np.array([force]),
        np.array([virial]),
    )
