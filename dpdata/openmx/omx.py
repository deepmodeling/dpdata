#!/usr/bin/python3
import warnings

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

from collections import OrderedDict

### iterout.c from OpenMX soure code: column numbers and physical quantities ###
# /* 1: */
# /* 2,3,4: */
# /* 5,6,7: force *   
# /* 8: x-component of velocity */
# /* 9: y-component of velocity */
# /* 10: z-component of velocity */
# /* 11: Net charge, electron charge is defined to be negative. */
# /* 12: magnetic moment (muB) */
# /* 13,14: angles of spin */  

# load atom_names, atom_numbs, atom_types, cells
def load_param_file(fname, mdname):
    with open(fname) as dat_file:
        lines = dat_file.readlines()
    atom_names = []
    atom_names_mode = False
    for line in lines:
        if "<Atoms.SpeciesAndCoordinates" in line:
            atom_names_mode = True
        elif "Atoms.SpeciesAndCoordinates>" in line:
            atom_names_mode = False
        elif atom_names_mode:
            parts = line.split()
            atom_names.append(parts[1])
    natoms = len(atom_names)
    atom_names_original=atom_names
    atom_names = list(OrderedDict.fromkeys(set(atom_names)))  # Python>=3.7
    atom_names = sorted(atom_names, key=atom_names_original.index) # Unique ordering of atomic species
    ntypes = len(atom_names)
    atom_numbs = [0] * ntypes
    atom_types = []
    atom_types_mode = False
    for line in lines:
        if "<Atoms.SpeciesAndCoordinates" in line:
            atom_types_mode = True
        elif "Atoms.SpeciesAndCoordinates>" in line:
            atom_types_mode = False
        elif atom_types_mode:
            parts = line.split()
            for i, atom_name in enumerate(atom_names):
                if parts[1] == atom_name:
                    atom_numbs[i] += 1
                    atom_types.append(i)
    atom_types=np.array(atom_types)

    with open(mdname) as md_file:
        lines = md_file.readlines()
    cell, cells = [], []
    for index, line in enumerate(lines):
        if "Cell_Vectors=" in line:
            parts = line.split()
            cell.append([float(parts[12]), float(parts[13]), float(parts[14])])
            cell.append([float(parts[15]), float(parts[16]), float(parts[17])])
            cell.append([float(parts[18]), float(parts[19]), float(parts[20])])
            cells.append(cell)
            cell = []
    cells = np.array(cells)
    return atom_names, atom_numbs, atom_types, cells

# load coords, energies, forces
def load_data(mdname, atom_names, natoms):
    with open(mdname) as md_file:
        lines = md_file.readlines()
    cnt = 0
    coord, coords = [], []
    for index, line in enumerate(lines):
        if "time=" in line:
            continue
        for atom_name in atom_names:
            atom_name += " "
            if atom_name in line:
                cnt += 1
                parts = line.split()
                for_line = [float(parts[1]), float(parts[2]), float(parts[3])]
                coord.append(for_line)
        if cnt == natoms:
            coords.append(coord)
            cnt = 0
            coord = []
    coords = np.array(coords)
    steps = [str(i) for i in range(1, coords.shape[0] + 1)]
    return coords, steps

def to_system_data(fname, mdname):
    data = {}
    (
        data["atom_names"],
        data["atom_numbs"],
        data["atom_types"],
        data["cells"],
    ) = load_param_file(fname, mdname)
    data["coords"], steps = load_data(
        mdname,
        data["atom_names"],
        np.sum(data["atom_numbs"]),
    )
    data["orig"] = np.zeros(3)
    return data, steps


def to_system_label(fname, mdname, data):
    atom_names, atom_numbs, atom_types, cells = load_param_file(fname, mdname)
    with open(mdname) as md_file:
        lines = md_file.readlines()
    cnt = 0
    energy = []
    field, fields = [], []
    for index, line in enumerate(lines):
        if "time=" in line:
            parts = line.split()
            ene_line = float(parts[4])  # Hartree
            energy.append(ene_line)
            continue
        for atom_name in atom_names:
            atom_name += " "
            if atom_name in line:
                cnt += 1
                parts = line.split()
                for_line = [float(parts[4]), float(parts[5]), float(parts[6])]
                field.append(for_line)
        if cnt == np.sum(data["atom_numbs"]):
            fields.append(field)
            cnt = 0
            field = []
    energy = energy_convert * np.array(energy)  # Hartree -> eV
    force = force_convert * np.array(fields)
    return energy, force


if __name__ == "__main__":
    label = 1
    file_name = "Methane"
    fname = f"../../{label}0.data/{file_name}.dat"
    mdname = f"../../{label}0.data/{file_name}.md"
    atom_names, atom_numbs, atom_types, cells = load_param_file(fname, mdname)
    coords, steps = load_data(mdname, atom_names, np.sum(atom_numbs))
    data, steps = to_system_data(fname, mdname)
    energy, force = to_system_label(fname, mdname, data)
    print(atom_names)
    print(atom_numbs)
    print(atom_types)
    # print(cells.shape)
    # print(coords.shape)
    # print(len(energy))
    # print(force.shape)