from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from dpdata.utils import open_file

if TYPE_CHECKING:
    from dpdata.utils import FileType


def read_orca_sp_output(
    fn: FileType,
) -> tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    """Read from ORCA output.

    Note that both the energy and the gradient should be printed.

    Parameters
    ----------
    fn : str
        file name

    Returns
    -------
    np.ndarray
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
    with open_file(fn) as f:
        flag = 0
        for line in f:
            if flag in (1, 3, 4):
                flag += 1
            elif flag == 2:
                s = line.split()
                if not len(s):
                    flag = 0
                else:
                    symbols.append(s[0].capitalize())
                    coord.append([float(s[1]), float(s[2]), float(s[3])])
            elif flag == 5:
                s = line.split()
                if not len(s):
                    flag = 0
                else:
                    forces.append([float(s[3]), float(s[4]), float(s[5])])
            elif line.startswith("CARTESIAN COORDINATES (ANGSTROEM)"):
                # coord
                flag = 1
                coord = []
                symbols = []
            elif line.startswith("CARTESIAN GRADIENT"):
                flag = 3
                forces = []
            elif line.startswith("FINAL SINGLE POINT ENERGY"):
                energy = float(line.split()[-1])
    symbols = np.array(symbols)
    forces = -np.array(forces)
    coord = np.array(coord)
    assert coord.shape == forces.shape

    return symbols, coord, energy, forces
