from __future__ import annotations

import os
import re

import numpy as np

from dpdata.amber.mask import pick_by_amber_mask
from dpdata.unit import EnergyConversion
from dpdata.utils import open_file

from ..periodic_table import ELEMENTS

kcalmol2eV = EnergyConversion("kcal_mol", "eV").value()
symbols = ["X"] + ELEMENTS

energy_convert = kcalmol2eV
force_convert = energy_convert


def cell_lengths_angles_to_cell(
    cell_lengths: np.ndarray, cell_angles: np.ndarray
) -> np.ndarray:
    """Convert cell lengths and angles to cell vectors.

    Parameters
    ----------
    cell_lengths
        Cell lengths with shape ``(..., 3)`` where the last dimension is
        ``a, b, c``.
    cell_angles
        Cell angles in degrees with shape ``(..., 3)`` where the last dimension
        is ``alpha, beta, gamma``.

    Returns
    -------
    np.ndarray
        Cell vectors with shape ``(..., 3, 3)``.
    """
    alpha = np.deg2rad(cell_angles[..., 0])
    beta = np.deg2rad(cell_angles[..., 1])
    gamma = np.deg2rad(cell_angles[..., 2])

    a = cell_lengths[..., 0]
    b = cell_lengths[..., 1]
    c = cell_lengths[..., 2]

    cos_alpha = np.cos(alpha)
    cos_beta = np.cos(beta)
    cos_gamma = np.cos(gamma)
    sin_gamma = np.sin(gamma)
    if np.any(np.isclose(sin_gamma, 0.0)):
        raise RuntimeError("Invalid AMBER cell angles")

    z_factor = (
        1
        - cos_alpha**2
        - cos_beta**2
        - cos_gamma**2
        + 2 * cos_alpha * cos_beta * cos_gamma
    )
    if np.any(z_factor <= 0.0):
        raise RuntimeError("Invalid AMBER cell angles")

    z = np.sqrt(z_factor) / sin_gamma

    shape = cell_lengths.shape[:-1] + (3, 3)
    cells = np.zeros(shape)
    cells[..., 0, 0] = a
    cells[..., 1, 0] = b * cos_gamma
    cells[..., 1, 1] = b * sin_gamma
    cells[..., 2, 0] = c * cos_beta
    cells[..., 2, 1] = c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma
    cells[..., 2, 2] = c * z
    return cells


def read_amber_traj(
    parm7_file,
    nc_file,
    mdfrc_file=None,
    mden_file=None,
    mdout_file=None,
    use_element_symbols=None,
    labeled=True,
):
    """The amber trajectory includes:
    * nc, NetCDF format, stores coordinates
    * mdfrc, NetCDF format, stores forces
    * mden (optional), text format, stores energies
    * mdout (optional), text format, may store energies if there is no mden_file
    * parm7, text format, stores types.

    Parameters
    ----------
    parm7_file, nc_file, mdfrc_file, mden_file, mdout_file:
        filenames
    use_element_symbols : None or list or str
        If use_element_symbols is a list of atom indexes, these atoms will use element symbols
        instead of amber types. For example, a ligand will use C, H, O, N, and so on
        instead of h1, hc, o, os, and so on.
        IF use_element_symbols is str, it will be considered as Amber mask.
    labeled : bool
        Whether to return labeled data
    """
    from scipy.io import netcdf_file

    flag_atom_type = False
    flag_atom_numb = False
    amber_types = []
    atomic_number = []
    with open_file(parm7_file) as f:
        for line in f:
            if line.startswith("%FLAG"):
                flag_atom_type = line.startswith("%FLAG AMBER_ATOM_TYPE")
                flag_atom_numb = (use_element_symbols is not None) and line.startswith(
                    "%FLAG ATOMIC_NUMBER"
                )
            elif flag_atom_type or flag_atom_numb:
                if line.startswith("%FORMAT"):
                    fmt = re.findall(r"\d+", line)
                    fmt0 = int(fmt[0])
                    fmt1 = int(fmt[1])
                else:
                    for ii in range(fmt0):
                        start_index = ii * fmt1
                        end_index = (ii + 1) * fmt1
                        if end_index >= len(line):
                            continue
                        content = line[start_index:end_index].strip()
                        if flag_atom_type:
                            amber_types.append(content)
                        elif flag_atom_numb:
                            atomic_number.append(int(content))
    if use_element_symbols is not None:
        if isinstance(use_element_symbols, str):
            use_element_symbols = pick_by_amber_mask(parm7_file, use_element_symbols)
        for ii in use_element_symbols:
            amber_types[ii] = symbols[atomic_number[ii]]

    with netcdf_file(nc_file, "r") as f:
        coords = np.array(f.variables["coordinates"][:])
        cell_lengths = np.array(f.variables["cell_lengths"][:])
        cell_angles = np.array(f.variables["cell_angles"][:])
        cells = cell_lengths_angles_to_cell(cell_lengths, cell_angles)

    if labeled:
        with netcdf_file(mdfrc_file, "r") as f:
            forces = np.array(f.variables["forces"][:])

        # load energy from mden_file or mdout_file
        energies = []
        if mden_file is not None and os.path.isfile(mden_file):
            with open_file(mden_file) as f:
                for line in f:
                    if line.startswith("L6"):
                        s = line.split()
                        if s[2] != "E_pot":
                            energies.append(float(s[2]))
        elif mdout_file is not None and os.path.isfile(mdout_file):
            with open_file(mdout_file) as f:
                for line in f:
                    if "EPtot" in line:
                        s = line.split()
                        energies.append(float(s[-1]))
        else:
            raise RuntimeError("Please provide one of mden_file and mdout_file")

    atom_names, atom_types, atom_numbs = np.unique(
        amber_types, return_inverse=True, return_counts=True
    )

    data = {}
    data["atom_names"] = list(atom_names)
    data["atom_numbs"] = list(atom_numbs)
    data["atom_types"] = atom_types
    if labeled:
        data["forces"] = forces * force_convert
        data["energies"] = np.array(energies) * energy_convert
    data["coords"] = coords
    data["cells"] = cells
    data["orig"] = np.array([0, 0, 0])
    return data
