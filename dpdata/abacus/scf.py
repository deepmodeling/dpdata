from __future__ import annotations

import os
import re
import warnings

import numpy as np

from dpdata.utils import open_file

from ..unit import LengthConversion, PressureConversion
from .stru import get_frame_from_stru

bohr2ang = LengthConversion("bohr", "angstrom").value()
kbar2evperang3 = PressureConversion("kbar", "eV/angstrom^3").value()


def CheckFile(ifile):
    if not os.path.isfile(ifile):
        print(f"Can not find file {ifile}")
        return False
    return True


def get_geometry_in(fname, inlines):
    geometry_path_in = os.path.join(fname, "STRU")
    for line in inlines:
        if "stru_file" in line and "stru_file" == line.split()[0]:
            atom_file = line.split()[1]
            geometry_path_in = os.path.join(fname, atom_file)
            break
    return geometry_path_in


def get_path_out(fname, inlines):
    path_out = os.path.join(fname, "OUT.ABACUS/running_scf.log")
    for line in inlines:
        if "suffix" in line and "suffix" == line.split()[0]:
            suffix = line.split()[1]
            path_out = os.path.join(fname, f"OUT.{suffix}/running_scf.log")
            break
    return path_out


def get_energy(outlines):
    Etot = None
    for line in reversed(outlines):
        if "final etot is" in line:
            Etot = float(line.split()[-2])  # in eV
            return Etot, True
        elif "convergence has NOT been achieved!" in line:
            return Etot, False
        elif "convergence has not been achieved" in line:
            return Etot, False

    return Etot, False


def collect_force(outlines):
    force = []
    for i, line in enumerate(outlines):
        if "TOTAL-FORCE (eV/Angstrom)" in line:
            value_pattern = re.compile(
                r"^\s*[A-Z][a-z]?[1-9][0-9]*\s+[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s+[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s+[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s*$"
            )
            j = i
            # find the first line of force
            noforce = False
            while not value_pattern.match(outlines[j]):
                j += 1
                if (
                    j >= i + 10
                ):  # if can not find the first line of force in 10 lines, then stop
                    warnings.warn("Warning: can not find the first line of force")
                    noforce = True
                    break
            if noforce:
                break

            force.append([])
            while value_pattern.match(outlines[j]):
                force[-1].append([float(ii) for ii in outlines[j].split()[1:4]])
                j += 1
    return force  # only return the last force


def get_force(outlines, natoms):
    force = collect_force(outlines)
    if len(force) == 0:
        return None
    else:
        return np.array(force[-1])  # only return the last force


def collect_stress(outlines):
    stress = []
    for i, line in enumerate(outlines):
        if "TOTAL-STRESS (KBAR)" in line:
            value_pattern = re.compile(
                r"^\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s+[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s+[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s*$"
            )
            j = i
            nostress = False
            while not value_pattern.match(outlines[j]):
                j += 1
                if (
                    j >= i + 10
                ):  # if can not find the first line of stress in 10 lines, then stop
                    warnings.warn("Warning: can not find the first line of stress")
                    nostress = True
                    break
            if nostress:
                break

            stress.append([])
            while value_pattern.match(outlines[j]):
                stress[-1].append(
                    list(map(lambda x: float(x), outlines[j].split()[0:3]))
                )
                j += 1
    return stress


def get_stress(outlines):
    stress = collect_stress(outlines)
    if len(stress) == 0:
        return None
    else:
        return np.array(stress[-1]) * kbar2evperang3  # only return the last stress


def get_mag_force(outlines):
    """Read atomic magmom and magnetic force from OUT.ABACUS/running_scf.log.

    Returns
    -------
    magmom: list of list of atomic magnetic moments (three dimensions: ION_STEP * NATOMS * 1/3)
    magforce: list of list of atomic magnetic forces (three dimensions: ION_STEP * NATOMS * 1/3)
    e.g.:
    -------------------------------------------------------------------------------------------
    Total Magnetism (uB)
    -------------------------------------------------------------------------------------------
        Fe         0.0000000001         0.0000000000         3.0000000307
        Fe        -0.0000000000        -0.0000000000         3.0000001151
    -------------------------------------------------------------------------------------------
    -------------------------------------------------------------------------------------------
    Magnetic force (eV/uB)
    -------------------------------------------------------------------------------------------
        Fe         0.0000000000         0.0000000000        -1.2117698671
        Fe         0.0000000000         0.0000000000        -1.2117928796
    -------------------------------------------------------------------------------------------

    """
    mags = []
    magforces = []
    for i, line in enumerate(outlines):
        if "Total Magnetism (uB)" in line:
            j = i + 2
            mag = []
            while "-------------------------" not in outlines[j]:
                imag = [float(ii) for ii in outlines[j].split()[1:]]
                if len(imag) == 1:
                    imag = [0, 0, imag[0]]
                mag.append(imag)
                j += 1
            mags.append(mag)
        if "Magnetic force (eV/uB)" in line:
            j = i + 2
            magforce = []
            while "-------------------------" not in outlines[j]:
                imagforce = [float(ii) for ii in outlines[j].split()[1:]]
                if len(imagforce) == 1:
                    imagforce = [0, 0, imagforce[0]]
                magforce.append(imagforce)
                j += 1
            magforces.append(magforce)
    return np.array(mags), np.array(magforces)


def get_frame(fname):
    data = {
        "atom_names": [],
        "atom_numbs": [],
        "atom_types": [],
        "cells": np.array([]),
        "coords": np.array([]),
        "energies": np.array([]),
        "forces": np.array([]),
    }

    if isinstance(fname, str):
        # if the input parameter is only one string, it is assumed that it is the
        # base directory containing INPUT file;
        path_in = os.path.join(fname, "INPUT")
    else:
        raise RuntimeError("invalid input")

    if not CheckFile(path_in):
        return data

    with open_file(path_in) as fp:
        inlines = fp.read().split("\n")

    geometry_path_in = get_geometry_in(fname, inlines)

    # get OUT.ABACUS/running_scf.log
    path_out = get_path_out(fname, inlines)
    if not (CheckFile(geometry_path_in) and CheckFile(path_out)):
        return data
    with open_file(path_out) as fp:
        outlines = fp.read().split("\n")

    # get energy
    energy, converge = get_energy(outlines)
    if not converge:
        return data

    # read STRU file
    data = get_frame_from_stru(geometry_path_in)
    natoms = sum(data["atom_numbs"])
    # should remove spins from STRU file
    if "spins" in data:
        data.pop("spins")
    move = data.pop("move", None)

    # get magmom and magforce, force and stress
    magmom, magforce = get_mag_force(outlines)
    if len(magmom) > 0:
        magmom = magmom[-1:]
    if len(magforce) > 0:
        magforce = magforce[-1:]

    force = get_force(outlines, natoms)
    stress = get_stress(outlines)

    data["energies"] = np.array(energy)[np.newaxis]
    data["forces"] = np.empty((0,)) if force is None else force[np.newaxis, :, :]
    data["orig"] = np.zeros(3)
    if stress is not None:
        cell = data["cells"][0]
        stress *= np.abs(np.linalg.det(cell))
        data["virials"] = stress[np.newaxis, :, :]

    if len(magmom) > 0:
        data["spins"] = magmom
    if len(magforce) > 0:
        data["force_mags"] = magforce
    if move is not None:
        data["move"] = move
    return data
