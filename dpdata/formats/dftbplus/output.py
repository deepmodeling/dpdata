from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from dpdata.utils import open_file

if TYPE_CHECKING:
    from dpdata.utils import FileType


def read_dftb_plus(
    fn_1: FileType, fn_2: FileType
) -> tuple[str, np.ndarray, float, np.ndarray]:
    """Read from DFTB+ input and output.

    Parameters
    ----------
    fn_1 : str
        DFTB+ input file name
    fn_2 : str
        DFTB+ output file name

    Returns
    -------
    str
        atomic symbols
    np.ndarray
        atomic coordinates
    float
        total potential energy
    np.ndarray
        atomic forces

    """
    coord = None
    symbols = None
    forces = None
    energy = None
    natoms = None
    with open_file(fn_1) as f:
        lines = iter(f)
        for line in lines:
            if not line.startswith("Geometry"):
                continue
            # GenFormat declares the atom count on the first line after the
            # opening brace and the element table on the next line.  The old
            # state machine hard-coded four geometry rows, so ammonia-sized
            # fixtures passed while larger molecules were truncated.
            count_line = next(lines).split()
            natoms = int(count_line[0])
            coordinate_mode = count_line[1].upper()
            components = next(lines).split()
            coord = []
            symbols = []
            for _ in range(natoms):
                s = next(lines).split()
                symbols.append(components[int(s[1]) - 1])
                coord.append([float(s[2]), float(s[3]), float(s[4])])
            if coordinate_mode == "F":
                # Fractional GenFormat coordinates are expressed in the
                # lattice-vector basis that follows the atom records.
                origin = np.array([float(value) for value in next(lines).split()])
                lattice = np.array(
                    [[float(value) for value in next(lines).split()] for _ in range(3)]
                )
                coord = (np.asarray(coord) @ lattice + origin).tolist()
            elif coordinate_mode not in {"C", "S"}:
                raise ValueError(
                    f"unsupported GenFormat coordinate mode: {coordinate_mode}"
                )
            break
    if natoms is None:
        raise ValueError("GenFormat Geometry block not found in DFTB+ input")
    with open_file(fn_2) as f:
        lines = iter(f)
        for line in lines:
            if line.startswith("Total Forces"):
                forces = []
                for _ in range(natoms):
                    s = next(lines).split()
                    forces.append([float(s[1]), float(s[2]), float(s[3])])
            elif line.startswith("Total energy:"):
                s = line.split()
                energy = float(s[2])

    symbols = np.array(symbols)
    forces = np.array(forces)
    coord = np.array(coord)
    assert coord.shape == forces.shape

    return symbols, coord, energy, forces
