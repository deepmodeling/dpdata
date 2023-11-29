#!/usr/bin/python3
import warnings

import numpy as np

from abc import ABC

from scipy import constants

AVOGADRO = constants.Avogadro  # Avagadro constant
ELE_CHG = constants.elementary_charge  # Elementary Charge, in C
BOHR = constants.value("atomic unit of length")  # Bohr, in m
HARTREE = constants.value("atomic unit of energy")  # Hartree, in Jole
RYDBERG = constants.Rydberg * constants.h * constants.c  # Rydberg, in Jole

# energy conversions
econvs = {
    "eV": 1.0,
    "hartree": HARTREE / ELE_CHG,
    "kJ_mol": 1 / (ELE_CHG * AVOGADRO / 1000),
    "kcal_mol": 1 / (ELE_CHG * AVOGADRO / 1000 / 4.184),
    "rydberg": RYDBERG / ELE_CHG,
    "J": 1 / ELE_CHG,
    "kJ": 1000 / ELE_CHG,
}

# length conversions
lconvs = {
    "angstrom": 1.0,
    "bohr": BOHR * 1e10,
    "nm": 10.0,
    "m": 1e10,
}


def check_unit(unit):
    if unit not in econvs.keys() and unit not in lconvs.keys():
        try:
            eunit = unit.split("/")[0]
            lunit = unit.split("/")[1]
            if eunit not in econvs.keys():
                raise RuntimeError(f"Invaild unit: {unit}")
            if lunit not in lconvs.keys():
                raise RuntimeError(f"Invalid unit: {unit}")
        except Exception:
            raise RuntimeError(f"Invalid unit: {unit}")


class Conversion(ABC):
    def __init__(self, unitA, unitB, check=True):
        """Parent class for unit conversion.

        Parameters
        ----------
        unitA : str
            unit to be converted
        unitB : str
            unit which unitA is converted to, i.e. `1 unitA = self._value unitB`
        check : bool
            whether to check unit validity

        Examples
        --------
        >>> conv = Conversion("foo", "bar", check=False)
        >>> conv.set_value("10.0")
        >>> print(conv)
        1 foo = 10.0 bar
        >>> conv.value()
        10.0
        """
        if check:
            check_unit(unitA)
            check_unit(unitB)
        self.unitA = unitA
        self.unitB = unitB
        self._value = 0.0

    def value(self):
        return self._value

    def set_value(self, value):
        self._value = value

    def __repr__(self):
        return f"1 {self.unitA} = {self._value} {self.unitB}"

    def __str__(self):
        return self.__repr__()


class EnergyConversion(Conversion):
    def __init__(self, unitA, unitB):
        """Class for energy conversion.

        Examples
        --------
        >>> conv = EnergyConversion("eV", "kcal_mol")
        >>> conv.value()
        23.06054783061903
        """
        super().__init__(unitA, unitB)
        self.set_value(econvs[unitA] / econvs[unitB])


class LengthConversion(Conversion):
    def __init__(self, unitA, unitB):
        """Class for length conversion.

        Examples
        --------
        >>> conv = LengthConversion("angstrom", "nm")
        >>> conv.value()
        0.1
        """
        super().__init__(unitA, unitB)
        self.set_value(lconvs[unitA] / lconvs[unitB])


class ForceConversion(Conversion):
    def __init__(self, unitA, unitB):
        """Class for force conversion.

        Parameters
        ----------
        unitA, unitB : str
            in format of "energy_unit/length_unit"

        Examples
        --------
        >>> conv = ForceConversion("kJ_mol/nm", "eV/angstrom")
        >>> conv.value()
        0.0010364269656262175
        """
        super().__init__(unitA, unitB)
        econv = EnergyConversion(unitA.split("/")[0], unitB.split("/")[0]).value()
        lconv = LengthConversion(unitA.split("/")[1], unitB.split("/")[1]).value()
        self.set_value(econv / lconv)


class PressureConversion(Conversion):
    def __init__(self, unitA, unitB):
        """Class for pressure conversion.

        Parameters
        ----------
        unitA, unitB : str
            in format of "energy_unit/length_unit^3", or in `["Pa", "pa", "kPa", "kpa", "bar", "kbar"]`

        Examples
        --------
        >>> conv = PressureConversion("kbar", "eV/angstrom^3")
        >>> conv.value()
        0.0006241509074460763
        """
        super().__init__(unitA, unitB, check=False)
        unitA, factorA = self._convert_unit(unitA)
        unitB, factorB = self._convert_unit(unitB)
        eunitA, lunitA = self._split_unit(unitA)
        eunitB, lunitB = self._split_unit(unitB)
        econv = EnergyConversion(eunitA, eunitB).value() * factorA / factorB
        lconv = LengthConversion(lunitA, lunitB).value()
        self.set_value(econv / lconv**3)

    def _convert_unit(self, unit):
        if unit == "Pa" or unit == "pa":
            return "J/m^3", 1
        elif unit == "kPa" or unit == "kpa":
            return "kJ/m^3", 1
        elif unit == "GPa" or unit == "Gpa":
            return "kJ/m^3", 1e6
        elif unit == "bar":
            return "J/m^3", 1e5
        elif unit == "kbar":
            return "kJ/m^3", 1e5
        else:
            return unit, 1

    def _split_unit(self, unit):
        eunit = unit.split("/")[0]
        lunit = unit.split("/")[1][:-2]
        return eunit, lunit


# from ..unit import (
#     EnergyConversion,
#     ForceConversion,
#     LengthConversion,
#     PressureConversion,
# )

ry2ev = EnergyConversion("rydberg", "eV").value()
kbar2evperang3 = PressureConversion("kbar", "eV/angstrom^3").value()

length_convert = LengthConversion("bohr", "angstrom").value()
energy_convert = EnergyConversion("hartree", "eV").value()
force_convert = ForceConversion("hartree/bohr", "eV/angstrom").value()

### operation start ###

import os, sys
from collections import OrderedDict

### commentted out on 2023/11/26 ###
# def load_param_file(fname):
#     # atom_names=load_atom(fname, "Atoms.SpeciesAndCoordinates")
#     # cell=load_cell(fname, "Atoms.UnitVectors")
#     ### future request: read from .md not .dat ###
#     with open(fname, "r") as dat_file:
#         lines = dat_file.readlines()
#         atom_names, cell = [], []
#         atom_names_mode, cell_mode = False, False
#         for line in lines:
#             if "<Atoms.SpeciesAndCoordinates" in line:
#                 atom_names_mode = True
#             elif "Atoms.SpeciesAndCoordinates>" in line:
#                 atom_names_mode = False
#             elif atom_names_mode:
#                 parts = line.split()
#                 atom_names.append(parts[1])
#             elif "<Atoms.UnitVectors" in line:
#                 cell_mode = True
#             elif "Atoms.UnitVectors>" in line:
#                 cell_mode = False
#             elif cell_mode:
#                 parts = line.split()
#                 cell.append(parts)
#         natoms=len(atom_names)
#         atom_names = list(OrderedDict.fromkeys(set(atom_names))) #注: Python3.7以降
#         # atom_names = list(set(atom_names)) #注: Python3.7以前
#         ntypes=len(atom_names)
#         cell=np.array(cell).astype(float)
#         atom_numbs = [0] * ntypes
#         atom_types = []
#         coords_mode = False
#         for line in lines:
#             if "<Atoms.SpeciesAndCoordinates" in line:
#                 coords_mode = True
#             elif "Atoms.SpeciesAndCoordinates>" in line:
#                 coords_mode = False
#             elif coords_mode:
#                 parts = line.split()
#                 for i, atom_name in enumerate(atom_names):
#                     if parts[1] == atom_name:
#                         atom_numbs[i]+=1
#                         atom_types.append(i)
#         if natoms != len(atom_types):
#             raise ValueError("Input file is incorrect.")
#         else:
#             atom_types = np.array(atom_types)
#         ## checking output ##
#         # atom_names=symbols
#         # atom_numbs=[4,1]
#         # atom_types=np.array([1,0,0,0,0])
#         # cell=10*np.eye(3)
#     return atom_names, atom_numbs, atom_types, cell


### modified on 2023/11/26 ###
def load_param_file(fname, mdname):
    ### future request: read from .md not .dat ###
    with open(fname, "r") as dat_file:
        lines = dat_file.readlines()
        atom_names = []
        atom_names_mode, cell_mode = False, False
        for line in lines:
            if "<Atoms.SpeciesAndCoordinates" in line:
                atom_names_mode = True
            elif "Atoms.SpeciesAndCoordinates>" in line:
                atom_names_mode = False
            elif atom_names_mode:
                parts = line.split()
                atom_names.append(parts[1])
        natoms = len(atom_names)
        atom_names = list(OrderedDict.fromkeys(set(atom_names)))  # 注: Python3.7以降
        # atom_names = list(set(atom_names)) #注: Python3.7以前
        ntypes = len(atom_names)
        atom_numbs = [0] * ntypes
        atom_types = []
        coords_mode = False
        for line in lines:
            if "<Atoms.SpeciesAndCoordinates" in line:
                coords_mode = True
            elif "Atoms.SpeciesAndCoordinates>" in line:
                coords_mode = False
            elif coords_mode:
                parts = line.split()
                for i, atom_name in enumerate(atom_names):
                    if parts[1] == atom_name:
                        atom_numbs[i] += 1
                        atom_types.append(i)

    if natoms != len(atom_types):
        raise ValueError("Input file is incorrect.")
    else:
        atom_types = np.array(atom_types)

    with open(mdname, "r") as md_file:
        lines = md_file.readlines()
        cnt = 0
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

    ## Raise error if files are inconsistent ##
    # if (cells[0] != cell_init).all():
    # raise ValueError("There is a consistency between input files.")

    ## checking output ##
    # atom_names=['H','C']
    # atom_numbs=[4,1]
    # atom_types=np.array([1,0,0,0,0])
    # cell=10*np.eye(3)
    return atom_names, atom_numbs, atom_types, cells


def load_data(mdname, atom_names, natoms, begin=0, step=1, convert=1.0):
    with open(mdname, "r") as md_file:
        lines = md_file.readlines()
        cnt = 0
        coord, coords = [], []
        for index, line in enumerate(lines):
            # future request
            if "time" in line:
                continue
            ## depend on atom_numbs ##
            for atom_name in atom_names:
                atom_name += " "
                if atom_name in line:
                    # elif f"{atom_names[0]} " in line or f"{atom_names[1]} " in line:
                    # elif f"{atom_names[0]} " in line:
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
        ## checking output ##
        # coords = np.random.rand(200,5,3)
        # steps = [str(i) for i in range(1, coords.shape[0]+1)]
    return coords, steps


def to_system_data(fname, mdname, begin=0, step=1):
    data = {}
    (
        data["atom_names"],
        data["atom_numbs"],
        data["atom_types"],
        data["cells"],
    ) = load_param_file(fname, mdname)
    data["coords"], csteps = load_data(
        mdname,
        data["atom_names"],
        np.sum(data["atom_numbs"]),
        begin=begin,
        step=step,
        convert=length_convert,
    )
    data["orig"] = np.zeros(3)
    # data["cells"]= np.array([cell for _ in range(len(csteps))])
    return data, csteps


def to_system_label(fname, mdname, data, begin=0, step=1):
    atom_names, atom_numbs, atom_types, cells = load_param_file(fname, mdname)
    with open(mdname, "r") as md_file:
        lines = md_file.readlines()
        cnt = 0
        energy = []
        field, fields = [], []
        for index, line in enumerate(lines):
            if "time" in line:
                parts = line.split()
                evp_line = float(parts[4])  # Hartree
                energy.append(evp_line)
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
        energy = energy_convert * np.array(energy)  # Hartree->eV
        force = force_convert * np.array(fields)
        ## checking output ##
        # energy=np.array([-29933]*200)
        # force=np.random.rand(200,5,3)
        # esteps=[str(i) for i in range(1, 201)]
    return energy, force


if __name__ == "__main__":
    label = 1
    file_name = "Methane"
    fname = f"../../{label}0.data/{file_name}.dat"
    mdname = f"../../{label}0.data/{file_name}.md"
    atom_names, atom_numbs, atom_types, cells = load_param_file(fname, mdname)
    coords, steps = load_data(mdname, atom_names, np.sum(atom_numbs))
    data, csteps = to_system_data(fname, mdname, begin=0, step=1)
    energy, force = to_system_label(fname, mdname, data, begin=0, step=1)
    print(atom_names)
    print(atom_numbs)
    print(atom_types)
    # print(cells)
    print(cells.shape)
    print(coords.shape)
    # print(data["coords"].shape)
    # print(data["cells"])
    # print(data["cells"].shape)
    # print(csteps)
    # print(energy)
    # print(len(energy))
    # print(force)
    # print(force.shape)
