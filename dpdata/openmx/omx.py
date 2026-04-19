#!/usr/bin/python3
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from dpdata.utils import open_file

if TYPE_CHECKING:
    from dpdata.utils import FileType

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

import warnings
from collections import OrderedDict


def load_atom(lines):
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
    atom_names_original = atom_names
    atom_names = list(OrderedDict.fromkeys(set(atom_names)))
    atom_names = sorted(
        atom_names, key=atom_names_original.index
    )  # Unique ordering of atomic species
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
    atom_types = np.array(atom_types)
    return atom_names, atom_types, atom_numbs


def load_cells(lines):
    cells = []
    for line in lines:
        if "Cell_Vectors=" in line:
            part = line.split("Cell_Vectors=")[1]
            parts = part.split()
            values = list(map(float, parts[:9]))
            cell = [values[0:3], values[3:6], values[6:9]]
            cells.append(cell)
            # Checking SCF converged or not
            for token in line.split():
                if token.startswith("scf_conv="):
                    scf_conv = int(token.split("=")[1])
                    if scf_conv == 0:
                        warnings.warn("SCF not converged!", stacklevel=2)
    cells = np.array(cells)
    return cells


# load atom_names, atom_numbs, atom_types, cells
def load_param_file(fname: FileType, mdname: FileType):
    with open_file(fname) as dat_file:
        lines = dat_file.readlines()
    atom_names, atom_types, atom_numbs = load_atom(lines)

    with open_file(mdname) as md_file:
        lines = md_file.readlines()
    cells = load_cells(lines)
    return atom_names, atom_numbs, atom_types, cells


def load_coords(lines, atom_names, natoms):
    cnt = 0
    coord, coords = [], []
    for line in lines:
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
    return coords


def load_data(mdname: FileType, atom_names, natoms):
    with open_file(mdname) as md_file:
        lines = md_file.readlines()
    coords = load_coords(lines, atom_names, natoms)
    steps = [str(i) for i in range(1, coords.shape[0] + 1)]
    return coords, steps


def to_system_data(fname: FileType, mdname: FileType):
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


def load_energy(lines):
    energy = []
    for line in lines:
        if "time=" in line:
            parts = line.split()
            ene_line = float(parts[4])  # Hartree
            energy.append(ene_line)
            continue
    energy = energy_convert * np.array(energy)  # Hartree -> eV
    return energy


def load_force(lines, atom_names, atom_numbs):
    cnt = 0
    field, fields = [], []
    for line in lines:
        if "time=" in line:
            continue
        for atom_name in atom_names:
            atom_name += " "
            if atom_name in line:
                cnt += 1
                parts = line.split()
                for_line = [float(parts[4]), float(parts[5]), float(parts[6])]
                field.append(for_line)
        if cnt == np.sum(atom_numbs):
            fields.append(field)
            cnt = 0
            field = []
    force = force_convert * np.array(fields)
    return force


# load energy, force
def to_system_label(fname, mdname):
    atom_names, atom_numbs, atom_types, cells = load_param_file(fname, mdname)
    with open_file(mdname) as md_file:
        lines = md_file.readlines()
    energy = load_energy(lines)
    force = load_force(lines, atom_names, atom_numbs)
    return energy, force


if __name__ == "__main__":
    file_name = "Au111Surface"
    fname = f"{file_name}.dat"
    mdname = f"{file_name}.md"
    atom_names, atom_numbs, atom_types, cells = load_param_file(fname, mdname)
    coords, steps = load_data(mdname, atom_names, np.sum(atom_numbs))
    data, steps = to_system_data(fname, mdname)
    energy, force = to_system_label(fname, mdname)
    print(atom_names)
    print(atom_numbs)
    print(atom_types)
# print(cells.shape)
# print(coords.shape)
# print(len(energy))
# print(force.shape)
