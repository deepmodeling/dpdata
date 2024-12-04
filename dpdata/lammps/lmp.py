#!/usr/bin/env python3
from __future__ import annotations

import numpy as np

ptr_float_fmt = "%15.10f"
ptr_int_fmt = "%6d"
ptr_key_fmt = "%15s"


def _get_block(lines, keys):
    for idx in range(len(lines)):
        if keys in lines[idx]:
            break
    if idx == len(lines) - 1:
        return None
    idx_s = idx + 2
    idx = idx_s
    ret = []
    while True:
        if idx == len(lines) or len(lines[idx].split()) == 0:
            break
        else:
            ret.append(lines[idx])
        idx += 1
    return ret


def lmpbox2box(lohi, tilt):
    xy = tilt[0]
    xz = tilt[1]
    yz = tilt[2]
    orig = np.array([lohi[0][0], lohi[1][0], lohi[2][0]])
    lens = []
    for dd in range(3):
        lens.append(lohi[dd][1] - lohi[dd][0])
    xx = [lens[0], 0, 0]
    yy = [xy, lens[1], 0]
    zz = [xz, yz, lens[2]]
    return orig, np.array([xx, yy, zz])


def box2lmpbox(orig, box):
    lohi = np.zeros([3, 2])
    for dd in range(3):
        lohi[dd][0] = orig[dd]
    tilt = np.zeros(3)
    tilt[0] = box[1][0]
    tilt[1] = box[2][0]
    tilt[2] = box[2][1]
    lens = np.zeros(3)
    lens[0] = box[0][0]
    lens[1] = box[1][1]
    lens[2] = box[2][2]
    for dd in range(3):
        lohi[dd][1] = lohi[dd][0] + lens[dd]
    return lohi, tilt


def get_atoms(lines):
    return _get_block(lines, "Atoms")


def get_natoms(lines):
    for ii in lines:
        if "atoms" in ii:
            return int(ii.split()[0])
    return None


def get_natomtypes(lines):
    for ii in lines:
        if "atom types" in ii:
            return int(ii.split()[0])
    return None


def _atom_info_mol(line):
    vec = line.split()
    # idx, mole_type, atom_type, charge, x, y, z
    return (
        int(vec[0]),
        int(vec[1]),
        int(vec[2]),
        float(vec[3]),
        float(vec[4]),
        float(vec[5]),
        float(vec[6]),
    )


def _atom_info_atom(line):
    vec = line.split()
    # idx, atom_type, x, y, z
    return int(vec[0]), int(vec[1]), float(vec[2]), float(vec[3]), float(vec[4])


def get_natoms_vec(lines):
    atype = get_atype(lines)
    natoms_vec = []
    natomtypes = get_natomtypes(lines)
    for ii in range(natomtypes):
        natoms_vec.append(sum(atype == ii + 1))
    assert sum(natoms_vec) == get_natoms(lines)
    return natoms_vec


def get_atype(lines, type_idx_zero=False):
    alines = get_atoms(lines)
    atype = []
    for ii in alines:
        # idx, mt, at, q, x, y, z = _atom_info_mol(ii)
        idx, at, x, y, z = _atom_info_atom(ii)
        if type_idx_zero:
            atype.append(at - 1)
        else:
            atype.append(at)
    return np.array(atype, dtype=int)


def get_posi(lines):
    atom_lines = get_atoms(lines)
    posis = []
    for ii in atom_lines:
        # posis.append([float(jj) for jj in ii.split()[4:7]])
        posis.append([float(jj) for jj in ii.split()[2:5]])
    return np.array(posis)


def get_spins(lines):
    atom_lines = get_atoms(lines)
    if len(atom_lines[0].split()) < 8:
        return None
    spins_ori = []
    spins_norm = []
    for ii in atom_lines:
        iis = ii.split()
        spins_ori.append([float(jj) for jj in iis[5:8]])
        spins_norm.append([float(iis[-1])])
    return np.array(spins_ori) * np.array(spins_norm)


def get_lmpbox(lines):
    box_info = []
    tilt = np.zeros(3)
    for ii in lines:
        if "xlo" in ii and "xhi" in ii:
            box_info.append([float(ii.split()[0]), float(ii.split()[1])])
            break
    for ii in lines:
        if "ylo" in ii and "yhi" in ii:
            box_info.append([float(ii.split()[0]), float(ii.split()[1])])
            break
    for ii in lines:
        if "zlo" in ii and "zhi" in ii:
            box_info.append([float(ii.split()[0]), float(ii.split()[1])])
            break
    for ii in lines:
        if "xy" in ii and "xz" in ii and "yz" in ii:
            tilt = np.array([float(jj) for jj in ii.split()[0:3]])
    return box_info, tilt


def system_data(lines, type_map=None, type_idx_zero=True):
    system = {}
    system["atom_numbs"] = get_natoms_vec(lines)
    system["atom_names"] = []
    if type_map is None:
        for ii in range(len(system["atom_numbs"])):
            system["atom_names"].append("Type_%d" % ii)  # noqa: UP031
    else:
        assert len(type_map) >= len(system["atom_numbs"])
        for ii in range(len(system["atom_numbs"])):
            system["atom_names"].append(type_map[ii])
    lohi, tilt = get_lmpbox(lines)
    orig, cell = lmpbox2box(lohi, tilt)
    system["orig"] = np.array(orig)
    system["cells"] = [np.array(cell)]
    natoms = sum(system["atom_numbs"])
    system["atom_types"] = get_atype(lines, type_idx_zero=type_idx_zero)
    system["coords"] = [get_posi(lines)]
    system["cells"] = np.array(system["cells"])
    system["coords"] = np.array(system["coords"])

    spins = get_spins(lines)
    if spins is not None:
        system["spins"] = np.array([spins])

    return system


def to_system_data(lines, type_map=None, type_idx_zero=True):
    return system_data(lines, type_map=type_map, type_idx_zero=type_idx_zero)


def from_system_data(system, f_idx=0):
    ret = ""
    ret += "\n"
    natoms = sum(system["atom_numbs"])
    ntypes = len(system["atom_numbs"])
    ret += "%d atoms\n" % natoms  # noqa: UP031
    ret += "%d atom types\n" % ntypes  # noqa: UP031
    ret += (ptr_float_fmt + " " + ptr_float_fmt + " xlo xhi\n") % (
        0,
        system["cells"][f_idx][0][0],
    )  # noqa: UP031
    ret += (ptr_float_fmt + " " + ptr_float_fmt + " ylo yhi\n") % (
        0,
        system["cells"][f_idx][1][1],
    )  # noqa: UP031
    ret += (ptr_float_fmt + " " + ptr_float_fmt + " zlo zhi\n") % (
        0,
        system["cells"][f_idx][2][2],
    )  # noqa: UP031
    ret += (
        ptr_float_fmt + " " + ptr_float_fmt + " " + ptr_float_fmt + " xy xz yz\n"
    ) % (
        system["cells"][f_idx][1][0],
        system["cells"][f_idx][2][0],
        system["cells"][f_idx][2][1],
    )  # noqa: UP031
    ret += "\n"
    ret += "Atoms # atomic\n"
    ret += "\n"
    coord_fmt = (
        ptr_int_fmt
        + " "
        + ptr_int_fmt
        + " "
        + ptr_float_fmt
        + " "
        + ptr_float_fmt
        + " "
        + ptr_float_fmt
        + "\n"
    )  # noqa: UP031

    if "spins" in system:
        coord_fmt = (
            coord_fmt.strip("\n")
            + " "
            + ptr_float_fmt
            + " "
            + ptr_float_fmt
            + " "
            + ptr_float_fmt
            + " "
            + ptr_float_fmt
            + "\n"
        )  # noqa: UP031
        spins_norm = np.linalg.norm(system["spins"][f_idx], axis=1)
    for ii in range(natoms):
        if "spins" in system:
            if spins_norm[ii] != 0:
                ret += coord_fmt % (
                    ii + 1,
                    system["atom_types"][ii] + 1,
                    system["coords"][f_idx][ii][0] - system["orig"][0],
                    system["coords"][f_idx][ii][1] - system["orig"][1],
                    system["coords"][f_idx][ii][2] - system["orig"][2],
                    system["spins"][f_idx][ii][0] / spins_norm[ii],
                    system["spins"][f_idx][ii][1] / spins_norm[ii],
                    system["spins"][f_idx][ii][2] / spins_norm[ii],
                    spins_norm[ii],
                )  # noqa: UP031
            else:
                ret += coord_fmt % (
                    ii + 1,
                    system["atom_types"][ii] + 1,
                    system["coords"][f_idx][ii][0] - system["orig"][0],
                    system["coords"][f_idx][ii][1] - system["orig"][1],
                    system["coords"][f_idx][ii][2] - system["orig"][2],
                    system["spins"][f_idx][ii][0],
                    system["spins"][f_idx][ii][1],
                    system["spins"][f_idx][ii][2] + 1,
                    spins_norm[ii],
                )  # noqa: UP031
        else:
            ret += coord_fmt % (
                ii + 1,
                system["atom_types"][ii] + 1,
                system["coords"][f_idx][ii][0] - system["orig"][0],
                system["coords"][f_idx][ii][1] - system["orig"][1],
                system["coords"][f_idx][ii][2] - system["orig"][2],
            )  # noqa: UP031
    return ret


if __name__ == "__main__":
    fname = "water-SPCE.data"
    lines = open(fname).read().split("\n")
    bonds, tilt = get_lmpbox(lines)
    # print(bonds, tilt)
    orig, box = lmpbox2box(bonds, tilt)
    # print(orig, box)
    bonds1, tilt1 = box2lmpbox(orig, box)
    # print(bonds1, tilt1)
    print(bonds1 - bonds)
    print(tilt1 - tilt)
    print(box)
    print(get_atype(lines))
    print(get_posi(lines))
