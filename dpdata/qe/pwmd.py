import os
import re
import sys

import numpy as np

from ..unit import (
    EnergyConversion,
    ForceConversion,
    LengthConversion,
    PressureConversion,
)

ry2ev = EnergyConversion("rydberg", "eV").value()
kbar2evperang3 = PressureConversion("kbar", "eV/angstrom^3").value()

length_convert = LengthConversion("bohr", "angstrom").value()
energy_convert = EnergyConversion("hartree", "eV").value()
force_convert = ForceConversion("hartree/bohr", "eV/angstrom").value()


def get_coords(lines, natoms):
    coord = []
    ret = []
    for i, string in enumerate(lines):
        if "ATOMIC_POSITIONS" in string:
            newlines = lines[i:]
            blk = get_block(newlines, "ATOMIC_POSITIONS")
            blk = blk[0 : sum(natoms)]
            for ii in blk:
                ret.append([float(jj) for jj in ii.split()[1:4]])
            coord.append(ret)
            ret = []
    coord = np.array(coord)
    return coord


def get_cell_vc(lines):
    cell = []
    ret = []
    for i, string in enumerate(lines):
        if "CELL_PARAMETERS" in string:
            newlines = lines[i:]
            blk = get_block(newlines, "CELL_PARAMETERS")
            for ii in blk:
                ret.append([float(jj) for jj in ii.split()[0:3]])
            cell.append(ret)
            ret = []
    cell = np.array(cell)
    return cell


def get_stress(lines):
    ret = []
    stress = []
    for i, string in enumerate(lines):
        if "total   stress" in string:
            newlines = lines[i:]
            blk = get_block(newlines, "total   stress")
            for ii in blk:
                ret.append([float(jj) for jj in ii.split()[3:6]])
            stress.append(ret)
            ret = []
    stress = np.array(stress)
    stress *= kbar2evperang3
    return stress


def get_energy(lines):
    energy = []
    for i, string in enumerate(lines):
        if "!    total energy" in string:
            energy.append(ry2ev * float(string.split("=")[1].split()[0]))
    energy = np.array(energy)

    return energy


def get_force(lines, natoms):
    ret = []
    force = []
    for i, string in enumerate(lines):
        if "Forces acting on atoms" in string:
            newlines = lines[i:]
            blk = get_block(newlines, "Forces acting on atoms", skip=1)
            blk = blk[0 : sum(natoms)]
            for ii in blk:
                ret.append([float(jj) for jj in ii.split("=")[1].split()])
            force.append(ret)
            ret = []
    force = np.array(force)
    force *= ry2ev / length_convert
    return force


def get_block(lines, keyword, skip=0):
    ret = []
    for idx, ii in enumerate(lines):
        if keyword in ii:
            blk_idx = idx + 1 + skip
            while len(lines[blk_idx]) == 0:
                blk_idx += 1
            while len(lines[blk_idx]) != 0 and blk_idx != len(lines):
                ret.append(lines[blk_idx])
                blk_idx += 1
            break
    return ret


def get_cell(lines):
    ret = []
    for idx, ii in enumerate(lines):
        if "ibrav" in ii:
            break
    blk = lines[idx : idx + 2]
    ibrav = int(blk[0].replace(",", "").split("=")[-1])
    if ibrav == 0:
        for iline in lines:
            if "CELL_PARAMETERS" in iline and "angstrom" not in iline.lower():
                raise RuntimeError(
                    "CELL_PARAMETERS must be written in Angstrom. Other units are not supported yet."
                )
        blk = get_block(lines, "CELL_PARAMETERS")
        for ii in blk:
            ret.append([float(jj) for jj in ii.split()[0:3]])
        ret = np.array(ret)
    elif ibrav == 1:
        a = None
        for iline in lines:
            line = iline.replace("=", " ").replace(",", "").split()
            if len(line) >= 2 and "a" == line[0]:
                # print("line = ", line)
                a = float(line[1])
            if len(line) >= 2 and "celldm(1)" == line[0]:
                a = float(line[1]) * length_convert
        # print("a = ", a)
        if not a:
            raise RuntimeError("parameter 'a' or 'celldm(1)' cannot be found.")
        ret = np.array([[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]])
    else:
        sys.exit("ibrav > 1 not supported yet.")
    return ret


def get_atoms(lines):
    atom_symbol_list = []
    for iline in lines:
        if "ATOMIC_POSITIONS" in iline:
            blk = get_block(lines, "ATOMIC_POSITIONS")
            for ii in blk:
                atom_symbol_list.append(ii.split()[0])

    atom_symbol_list = np.array(atom_symbol_list)
    tmp_names, symbol_idx = np.unique(atom_symbol_list, return_index=True)
    atom_types = []
    atom_numbs = []
    # preserve the atom_name order
    atom_names = atom_symbol_list[np.sort(symbol_idx)]
    for jj in atom_symbol_list:
        for idx, ii in enumerate(atom_names):
            if jj == ii:
                atom_types.append(idx)
    for idx in range(len(atom_names)):
        atom_numbs.append(atom_types.count(idx))
    atom_types = np.array(atom_types)
    return list(atom_names), atom_numbs, atom_types


def to_system_data(fname, begin=0, step=1):
    if type(fname) == str:
        path_out = fname
        outname = os.path.basename(path_out)
        # the name of the input file is assumed to be different from the output by 'in' and 'out'
        inname = outname.replace("out", "in")
        path_in = os.path.join(os.path.dirname(path_out), inname)
    elif type(fname) == list and len(fname) == 2:
        path_in = fname[0]
        path_out = fname[1]
    else:
        raise RuntimeError("invalid input")
    with open(path_out) as fp:
        outlines = fp.read().split("\n")
    with open(path_in) as fp:
        inlines = fp.read().split("\n")

    for i, string in enumerate(outlines):
        if "Program PWSCF" in string:
            checkpoint_index = i
    outlines = outlines[checkpoint_index:]
    # In case of output file from previous failed pw calculation was not deleted

    atom_names, atom_numbs, atom_types = get_atoms(inlines)

    coords = get_coords(outlines, atom_numbs)
    energies = get_energy(outlines)
    forces = get_force(outlines, atom_numbs)
    virials = get_stress(outlines)

    for line in inlines:
        # check calculation option in input file to find it is vc (variable cell) or not
        if "calculation" in line:
            calculation = line.strip()
    calculation = re.search(r"'([^']*)'", calculation)
    calculation = calculation.group(1)
    calculation = calculation.lower()
    if calculation == "md" or calculation == "relax":
        cells = get_cell(inlines)
        cells = np.tile(cells, (len(virials), 1, 1))
        for ii, temp in enumerate(virials):
            virials[ii] = virials[ii] * np.linalg.det(cells[ii])
    elif calculation == "vc-md" or calculation == "vc-relax":
        cells = get_cell_vc(outlines)
        for ii, temp in enumerate(cells):
            virials[ii] = virials[ii] * np.linalg.det(cells[ii])
    # because cell changes every step when variable cell calculation

    return (
        atom_names,
        atom_numbs,
        atom_types,
        cells,
        coords,
        energies,
        forces,
        virials,
    )
