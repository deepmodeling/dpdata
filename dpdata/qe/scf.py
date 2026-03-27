#!/usr/bin/env python3
from __future__ import annotations

import os

import numpy as np

from dpdata.utils import open_file

from .traj import (
    kbar2evperang3,
    ry2ev,
)
from .traj import (
    length_convert as bohr2ang,
)

_QE_BLOCK_KEYWORDS = [
    "ATOMIC_SPECIES",
    "ATOMIC_POSITIONS",
    "K_POINTS",
    "ADDITIONAL_K_POINTS",
    "CELL_PARAMETERS",
    "CONSTRAINTS",
    "OCCUPATIONS",
    "ATOMIC_VELOCITIES",
    "ATOMIC_FORCES",
    "SOLVENTS",
    "HUBBARD",
]


def get_block(lines, keyword, skip=0):
    ret = []
    for idx, ii in enumerate(lines):
        if keyword in ii:
            blk_idx = idx + 1 + skip
            while len(lines[blk_idx].split()) == 0:
                blk_idx += 1
            while (
                len(lines[blk_idx].split()) != 0
                and (lines[blk_idx].split()[0] not in _QE_BLOCK_KEYWORDS)
            ) and blk_idx != len(lines):
                ret.append(lines[blk_idx])
                blk_idx += 1
            break
    return ret


def _parse_lattice_parameters(lines):
    """Parse lattice parameters from QE input lines."""
    params = {}
    for iline in lines:
        line = iline.replace("=", " ").replace(",", "").split()
        if len(line) >= 2:
            if line[0] == "a":
                params["a"] = float(line[1])
            elif line[0] == "b":
                params["b"] = float(line[1])
            elif line[0] == "c":
                params["c"] = float(line[1])
            elif line[0] == "cosab":
                params["cosab"] = float(line[1])
            elif line[0] == "cosac":
                params["cosac"] = float(line[1])
            elif line[0] == "cosbc":
                params["cosbc"] = float(line[1])
            elif line[0].startswith("celldm("):
                # Extract celldm index from celldm(1), celldm(2), etc.
                idx = int(line[0][7:-1])  # Extract number from celldm(n)
                if "celldm" not in params:
                    params["celldm"] = {}
                params["celldm"][idx] = float(line[1])
    return params


def _convert_ibrav_to_cell(ibrav, params):
    """Convert ibrav and lattice parameters to cell matrix."""
    # Extract parameters
    a = params.get("a")
    b = params.get("b") 
    c = params.get("c")
    cosab = params.get("cosab", 0.0)
    cosac = params.get("cosac", 0.0)
    cosbc = params.get("cosbc", 0.0)
    celldm = params.get("celldm", {})
    
    # Convert celldm parameters if present (celldm is in atomic units)
    if 1 in celldm:
        a = celldm[1] * bohr2ang if a is None else a
    if 2 in celldm and a is not None:
        b = celldm[2] * a if b is None else b
    if 3 in celldm and a is not None:
        c = celldm[3] * a if c is None else c
    if 4 in celldm:
        cosab = celldm[4] if cosab == 0.0 else cosab
    if 5 in celldm:
        cosac = celldm[5] if cosac == 0.0 else cosac
    if 6 in celldm:
        cosbc = celldm[6] if cosbc == 0.0 else cosbc
    
    # Set defaults only for ibrav types that don't require explicit b,c
    if ibrav in [1, 2, 3, -3]:  # Cubic lattices
        if b is None:
            b = a
        if c is None:
            c = a
    elif ibrav in [6, 7]:  # Tetragonal lattices
        if b is None:
            b = a
        # c must be specified explicitly
    
    # Validate required parameters
    if a is None:
        raise RuntimeError("parameter 'a' or 'celldm(1)' cannot be found.")
    
    # Generate cell matrix based on ibrav
    if ibrav == 1:  # simple cubic
        return np.array([[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, a]])
    
    elif ibrav == 2:  # face-centered cubic
        return a * 0.5 * np.array([[-1, 0, 1], [0, 1, 1], [-1, 1, 0]])
    
    elif ibrav == 3:  # body-centered cubic
        return a * 0.5 * np.array([[1, 1, 1], [-1, 1, 1], [-1, -1, 1]])
        
    elif ibrav == -3:  # reverse body-centered cubic
        return a * 0.5 * np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]])
    
    elif ibrav == 4:  # hexagonal
        if c is None:
            raise RuntimeError("parameter 'c' or 'celldm(3)' required for ibrav=4")
        return np.array([[a, 0, 0], 
                        [-a/2, a*np.sqrt(3)/2, 0], 
                        [0, 0, c]])
    
    elif ibrav == 6:  # simple tetragonal  
        if c is None:
            raise RuntimeError("parameter 'c' or 'celldm(3)' required for ibrav=6")
        return np.array([[a, 0, 0], [0, a, 0], [0, 0, c]])
    
    elif ibrav == 7:  # body-centered tetragonal
        if c is None:
            raise RuntimeError("parameter 'c' or 'celldm(3)' required for ibrav=7")
        return np.array([[a/2, -a/2, c/2], 
                        [a/2, a/2, c/2], 
                        [-a/2, -a/2, c/2]])
    
    elif ibrav == 8:  # simple orthorhombic
        if b is None:
            raise RuntimeError("parameter 'b' or 'celldm(2)' required for ibrav=8")
        if c is None:
            raise RuntimeError("parameter 'c' or 'celldm(3)' required for ibrav=8")
        return np.array([[a, 0, 0], [0, b, 0], [0, 0, c]])
    
    elif ibrav == 9:  # base-centered orthorhombic
        if b is None:
            raise RuntimeError("parameter 'b' or 'celldm(2)' required for ibrav=9")
        if c is None:
            raise RuntimeError("parameter 'c' or 'celldm(3)' required for ibrav=9")
        return np.array([[a/2, b/2, 0], 
                        [-a/2, b/2, 0], 
                        [0, 0, c]])
    
    elif ibrav == -9:  # reverse base-centered orthorhombic
        if b is None:
            raise RuntimeError("parameter 'b' or 'celldm(2)' required for ibrav=-9")
        if c is None:
            raise RuntimeError("parameter 'c' or 'celldm(3)' required for ibrav=-9")
        return np.array([[a/2, -b/2, 0], 
                        [a/2, b/2, 0], 
                        [0, 0, c]])
    
    elif ibrav == 10:  # face-centered orthorhombic
        if b is None:
            raise RuntimeError("parameter 'b' or 'celldm(2)' required for ibrav=10")
        if c is None:
            raise RuntimeError("parameter 'c' or 'celldm(3)' required for ibrav=10")
        return np.array([[a/2, 0, c/2], 
                        [a/2, b/2, 0], 
                        [0, b/2, c/2]])
    
    elif ibrav == 11:  # body-centered orthorhombic
        if b is None:
            raise RuntimeError("parameter 'b' or 'celldm(2)' required for ibrav=11")
        if c is None:
            raise RuntimeError("parameter 'c' or 'celldm(3)' required for ibrav=11")
        return np.array([[a/2, b/2, c/2], 
                        [-a/2, b/2, c/2], 
                        [-a/2, -b/2, c/2]])
    
    elif ibrav == 12:  # simple monoclinic
        if b is None:
            raise RuntimeError("parameter 'b' or 'celldm(2)' required for ibrav=12")
        if c is None:
            raise RuntimeError("parameter 'c' or 'celldm(3)' required for ibrav=12")
        sinab = np.sqrt(1 - cosab**2)
        return np.array([[a, 0, 0], 
                        [b*cosab, b*sinab, 0], 
                        [0, 0, c]])
    
    elif ibrav == -12:  # reverse monoclinic 
        if b is None:
            raise RuntimeError("parameter 'b' or 'celldm(2)' required for ibrav=-12")
        if c is None:
            raise RuntimeError("parameter 'c' or 'celldm(3)' required for ibrav=-12")
        sinac = np.sqrt(1 - cosac**2)
        return np.array([[a, 0, 0], 
                        [0, b, 0], 
                        [c*cosac, 0, c*sinac]])
    
    elif ibrav == 13:  # base-centered monoclinic
        if b is None:
            raise RuntimeError("parameter 'b' or 'celldm(2)' required for ibrav=13")
        if c is None:
            raise RuntimeError("parameter 'c' or 'celldm(3)' required for ibrav=13")
        sinab = np.sqrt(1 - cosab**2)
        return np.array([[a/2, 0, -c/2], 
                        [b*cosab, b*sinab, 0], 
                        [a/2, 0, c/2]])
    
    elif ibrav == -13:  # reverse base-centered monoclinic
        if b is None:
            raise RuntimeError("parameter 'b' or 'celldm(2)' required for ibrav=-13")
        if c is None:
            raise RuntimeError("parameter 'c' or 'celldm(3)' required for ibrav=-13")
        sinac = np.sqrt(1 - cosac**2)
        return np.array([[a/2, -b/2, 0], 
                        [a/2, b/2, 0], 
                        [c*cosac, 0, c*sinac]])
    
    elif ibrav == 14:  # triclinic
        if b is None:
            raise RuntimeError("parameter 'b' or 'celldm(2)' required for ibrav=14")
        if c is None:
            raise RuntimeError("parameter 'c' or 'celldm(3)' required for ibrav=14")
        sinab = np.sqrt(1 - cosab**2)
        sinac = np.sqrt(1 - cosac**2)
        cosbc_prime = (cosbc - cosab*cosac) / (sinab*sinac)
        sinbc_prime = np.sqrt(1 - cosbc_prime**2)
        return np.array([[a, 0, 0], 
                        [b*cosab, b*sinab, 0], 
                        [c*cosac, c*sinac*cosbc_prime, c*sinac*sinbc_prime]])
    
    else:
        raise RuntimeError(f"ibrav = {ibrav} is not supported yet.")


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
    else:
        # Parse lattice parameters and convert based on ibrav
        params = _parse_lattice_parameters(lines)
        ret = _convert_ibrav_to_cell(ibrav, params)
    
    return ret


def get_coords(lines, cell):
    coord = []
    atom_symbol_list = []
    for iline in lines:
        if "ATOMIC_POSITIONS" in iline and (
            "angstrom" not in iline.lower() and "crystal" not in iline.lower()
        ):
            raise RuntimeError(
                "ATOMIC_POSITIONS must be written in Angstrom or crystal. Other units are not supported yet."
            )
        if "ATOMIC_POSITIONS" in iline and "angstrom" in iline.lower():
            blk = get_block(lines, "ATOMIC_POSITIONS")
            for ii in blk:
                coord.append([float(jj) for jj in ii.split()[1:4]])
                atom_symbol_list.append(ii.split()[0])
            coord = np.array(coord)
        elif "ATOMIC_POSITIONS" in iline and "crystal" in iline.lower():
            blk = get_block(lines, "ATOMIC_POSITIONS")
            for ii in blk:
                coord.append([float(jj) for jj in ii.split()[1:4]])
                atom_symbol_list.append(ii.split()[0])
            coord = np.array(coord)
            coord = np.matmul(coord, cell)
    atom_symbol_list = np.array(atom_symbol_list)
    tmp_names, symbol_idx = np.unique(atom_symbol_list, return_index=True)
    atom_types = []
    atom_numbs = []
    # preserve the atom_name order
    atom_names = atom_symbol_list[np.sort(symbol_idx, kind="stable")]
    for jj in atom_symbol_list:
        for idx, ii in enumerate(atom_names):
            if jj == ii:
                atom_types.append(idx)
    for idx in range(len(atom_names)):
        atom_numbs.append(atom_types.count(idx))
    atom_types = np.array(atom_types)

    return list(atom_names), atom_numbs, atom_types, coord


def get_energy(lines):
    energy = None
    for ii in lines:
        if "!    total energy" in ii:
            energy = ry2ev * float(ii.split("=")[1].split()[0])
    return energy


def get_force(lines, natoms):
    blk = get_block(lines, "Forces acting on atoms", skip=1)
    ret = []
    blk = blk[0 : sum(natoms)]
    for ii in blk:
        ret.append([float(jj) for jj in ii.split("=")[1].split()])
    ret = np.array(ret)
    ret *= ry2ev / bohr2ang
    return ret


def get_stress(lines):
    blk = get_block(lines, "total   stress")
    if len(blk) == 0:
        return None
    ret = []
    for ii in blk:
        ret.append([float(jj) for jj in ii.split()[3:6]])
    ret = np.array(ret)
    ret *= kbar2evperang3
    return ret


def get_frame(fname):
    if isinstance(fname, str):
        path_out = fname
        outname = os.path.basename(path_out)
        # the name of the input file is assumed to be different from the output by 'in' and 'out'
        inname = outname.replace("out", "in")
        path_in = os.path.join(os.path.dirname(path_out), inname)
    elif isinstance(fname, list) and len(fname) == 2:
        path_in = fname[0]
        path_out = fname[1]
    else:
        raise RuntimeError("invalid input")
    with open_file(path_out) as fp:
        outlines = fp.read().split("\n")
    with open_file(path_in) as fp:
        inlines = fp.read().split("\n")
    cell = get_cell(inlines)
    atom_names, natoms, types, coords = get_coords(inlines, cell)
    energy = get_energy(outlines)
    force = get_force(outlines, natoms)
    stress = get_stress(outlines)
    if stress is not None:
        stress = (stress * np.linalg.det(cell))[np.newaxis, :, :]
    return (
        atom_names,
        natoms,
        types,
        cell[np.newaxis, :, :],
        coords[np.newaxis, :, :],
        np.array(energy)[np.newaxis],
        force[np.newaxis, :, :],
        stress,
    )
