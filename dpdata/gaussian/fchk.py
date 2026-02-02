from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from dpdata.utils import open_file

if TYPE_CHECKING:
    from dpdata.utils import FileType

from ..periodic_table import ELEMENTS
from ..unit import (
    EnergyConversion,
    ForceConversion,
    HessianConversion,
    LengthConversion,
)

length_convert = LengthConversion("bohr", "angstrom").value()
energy_convert = EnergyConversion("hartree", "eV").value()
force_convert = ForceConversion("hartree/bohr", "eV/angstrom").value()
hessian_convert = HessianConversion("hartree/bohr^2", "eV/angstrom^2").value()


def create_full_hessian(hessian_raw: list | np.ndarray, natoms: int) -> np.ndarray:
    """
    Reconstructs the full, symmetric Hessian matrix from a 1D array
    containing its lower triangular elements.

    Args:
        hessian_raw (list | np.ndarray): A 1D list or NumPy array containing the
                                         lower triangular elements (including the
                                         diagonal) of the Hessian matrix.
        natoms (int): The number of atoms in the system.

    Returns
    -------
    np.ndarray: A full, symmetric (3*natoms, 3*natoms) Hessian matrix.

    Raises
    ------
    ValueError: If the number of elements in `hessian_raw` does not match
        the expected number for the lower triangle of a
        (3*natoms, 3*natoms) matrix.
    """
    # Convert input to a NumPy array in case it's a list
    hessian_block = np.array(hessian_raw)

    # Calculate the dimension of the final matrix
    dim = 3 * natoms

    # Validate that the input data has the correct length
    # A lower triangle of an n x n matrix has n*(n+1)/2 elements
    expected_length = dim * (dim + 1) // 2
    if hessian_block.size != expected_length:
        raise ValueError(
            f"Input length {hessian_block.size} != expected {expected_length}"
        )

    # Create a zero matrix, then fill the lower triangle
    hessian_full = np.zeros((dim, dim), dtype=hessian_block.dtype)
    lower_triangle_indices = np.tril_indices(dim)
    hessian_full[lower_triangle_indices] = hessian_block

    # This is done by copying the lower triangle to the upper triangle
    # M_full = M_lower + M_lower.T - diag(M_lower)
    hessian_full = hessian_full + hessian_full.T - np.diag(np.diag(hessian_full))

    return hessian_full


def to_system_data(file_name: FileType, has_forces=True, has_hessian=True):
    """Read Gaussian fchk file.

    Parameters
    ----------
    file_name : str
        file name
    has_forces : bool, default True
        whether to read force
        Note: Cartesian Gradient in fchk file is converted to forces by taking negative sign
    has_hessian : bool, default True
        whether to read hessian

    Returns
    -------
    data : dict
        system data, including hessian if has_hessian is True
    """
    data = {}
    natoms = 0
    atom_numbers = []
    coords_t = []
    energy_t = []
    forces_t = []
    hessian_t = []
    # Read fchk file
    with open_file(file_name) as fp:
        for line in fp:
            if isinstance(line, bytes):
                line = line.decode(errors="ignore")
            if "Number of atoms" in line:
                natoms = int(line.split()[-1])
            elif "Atomic numbers" in line and "I" in line:
                n = int(line.split()[-1])
                atom_numbers = []
                while len(atom_numbers) < n:
                    next_line = next(fp)
                    if isinstance(next_line, bytes):
                        next_line = next_line.decode(errors="ignore")
                    atom_numbers += [int(x) for x in next_line.split()]
            elif "Current cartesian coordinates" in line and "R" in line:
                n = int(line.split()[-1])
                coords_raw = []
                while len(coords_raw) < n:
                    next_line = next(fp)
                    if isinstance(next_line, bytes):
                        next_line = next_line.decode(errors="ignore")
                    coords_raw += [float(x) for x in next_line.split()]
                coords = np.array(coords_raw).reshape(-1, 3) * length_convert
                coords_t.append(coords)
            elif "Total Energy" in line:
                energy = float(line.split()[-1]) * energy_convert
                energy_t.append(energy)
            elif "Cartesian Gradient" in line:
                n = int(line.split()[-1])
                forces_raw = []
                while len(forces_raw) < n:
                    next_line = next(fp)
                    if isinstance(next_line, bytes):
                        next_line = next_line.decode(errors="ignore")
                    forces_raw += [float(x) for x in next_line.split()]
                # Cartesian Gradient is the negative of forces: F = -∇E
                forces = -np.array(forces_raw).reshape(-1, 3) * force_convert
                forces_t.append(forces)
            elif "Cartesian Force Constants" in line and "R" in line:
                n = int(line.split()[-1])
                hessian_raw = []
                while len(hessian_raw) < n:
                    next_line = next(fp)
                    if isinstance(next_line, bytes):
                        next_line = next_line.decode(errors="ignore")
                    hessian_raw += [float(x) for x in next_line.split()]
                hessian_full = (
                    create_full_hessian(hessian_raw, natoms) * hessian_convert
                )
                # store as (natoms, 3, natoms, 3) to align with registered shape
                hessian_t.append(hessian_full.reshape(natoms, 3, natoms, 3))
    # Assert key data
    assert coords_t, "cannot find coords"
    assert energy_t, "cannot find energy"
    if has_forces:
        assert forces_t, "cannot find forces"
    if has_hessian:
        assert hessian_t, "cannot find hessian"
    # Assemble data
    atom_symbols = [ELEMENTS[z - 1] for z in atom_numbers]
    atom_names, atom_types, atom_numbs = np.unique(
        atom_symbols, return_inverse=True, return_counts=True
    )
    data["atom_names"] = list(atom_names)
    data["atom_numbs"] = list(atom_numbs)
    data["atom_types"] = atom_types
    data["coords"] = np.array(coords_t).reshape(-1, natoms, 3)
    data["orig"] = np.zeros(3)
    data["cells"] = np.array([np.eye(3) * 100])
    data["nopbc"] = True
    if energy_t:
        data["energies"] = np.array(energy_t)
    if has_forces and forces_t:
        data["forces"] = np.array(forces_t)
    if has_hessian and hessian_t:
        data["hessian"] = np.array(hessian_t)
    return data
