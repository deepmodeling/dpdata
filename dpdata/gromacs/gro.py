#!/usr/bin/env python3
import re

import numpy as np

from ..unit import LengthConversion

nm2ang = LengthConversion("nm", "angstrom").value()
ang2nm = LengthConversion("angstrom", "nm").value()
cell_idx_gmx2dp = [0, 4, 8, 1, 2, 3, 5, 6, 7]


def _format_atom_name(atom_name):
    patt = re.compile("[a-zA-Z]*")
    match = re.search(patt, atom_name)
    fmt_name = match.group().capitalize()
    return fmt_name


def _get_line(line, fmt_atom_name=True):
    atom_name = line[10:15].split()[0]
    if fmt_atom_name:
        atom_name = _format_atom_name(atom_name)
    atom_idx = int(line[15:20].split()[0])
    posis = [float(line[ii : ii + 8]) for ii in range(20, 44, 8)]
    posis = np.array(posis) * nm2ang
    return atom_name, atom_idx, posis


def _get_cell(line):
    cell = np.zeros([3, 3])
    lengths = [float(ii) for ii in line.split()]
    if len(lengths) >= 3:
        for dd in range(3):
            cell[dd][dd] = lengths[dd]
    else:
        raise RuntimeError("wrong box format: ", line)
    if len(lengths) == 9:
        cell[0][1] = lengths[3]
        cell[0][2] = lengths[4]
        cell[1][0] = lengths[5]
        cell[1][2] = lengths[6]
        cell[2][0] = lengths[7]
        cell[2][1] = lengths[8]
    cell = cell * nm2ang
    return cell


def file_to_system_data(fname, format_atom_name=True, **kwargs):
    system = {"coords": [], "cells": []}
    with open(fname) as fp:
        frame = 0
        while True:
            flag = fp.readline()
            if not flag:
                break
            else:
                frame += 1
                names = []
                idxs = []
                posis = []
                natoms = int(fp.readline())
                for ii in range(natoms):
                    n, i, p = _get_line(fp.readline(), fmt_atom_name=format_atom_name)
                    names.append(n)
                    idxs.append(i)
                    posis.append(p)
                cell = _get_cell(fp.readline())
                posis = np.array(posis)
                if frame == 1:
                    system["orig"] = np.zeros(3)
                    system["atom_names"] = list(set(names))
                    system["atom_numbs"] = [
                        names.count(ii) for ii in system["atom_names"]
                    ]
                    system["atom_types"] = [
                        system["atom_names"].index(ii) for ii in names
                    ]
                    system["atom_types"] = np.array(system["atom_types"], dtype=int)
                system["coords"].append(posis)
                system["cells"].append(cell)
    system["coords"] = np.array(system["coords"])
    system["cells"] = np.array(system["cells"])
    return system


def from_system_data(system, f_idx=0, **kwargs):
    resname = kwargs.get("resname", "MOL")
    shift = kwargs.get("shift", 0)
    ret = ""
    ret += " molecule" + "\n"
    n_atoms = sum(system["atom_numbs"])
    ret += " " + str(n_atoms) + "\n"
    for i in range(n_atoms):
        atom_type = system["atom_types"][i]
        atom_name = system["atom_names"][atom_type]
        coords = system["coords"][f_idx] * ang2nm
        ret += "{:>5d}{:<5s}{:>5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(
            1, resname, atom_name, i + shift + 1, *tuple(coords[i])
        )
    cell = (system["cells"][f_idx].flatten() * ang2nm)[cell_idx_gmx2dp]
    ret += " " + " ".join([f"{x:.3f}" for x in cell])

    return ret
