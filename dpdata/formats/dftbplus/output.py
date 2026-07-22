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
    with open_file(fn_1) as f:
        flag = 0
        n_atoms = 0
        for line in f:
            if flag == 1:
                # Header line: "N C" where N is atom count, C is coordinate type
                header = line.split()
                n_atoms = int(header[0])
                flag += 1
            elif flag == 2:
                components = line.split()
                flag += 1
            elif line.startswith("Geometry"):
                flag = 1
                coord = []
                symbols = []
                n_atoms = 0
            elif flag == 3 and n_atoms > 0:
                # Parse exactly n_atoms coordinate rows
                s = line.split()
                components_num = int(s[1])
                symbols.append(components[components_num - 1])
                coord.append([float(s[2]), float(s[3]), float(s[4])])
                n_atoms -= 1
                if n_atoms == 0:
                    flag = 0
    with open_file(fn_2) as f:
        flag = 0
        n_forces = 0
        for line in f:
            if line.startswith("Total Forces"):
                flag = 8
                forces = []
                n_forces = 0
            elif flag == 8:
                # First line after "Total Forces" header contains force data
                # We don't know the count yet, so parse until we hit a non-force line
                s = line.split()
                if len(s) >= 4:
                    try:
                        # Try to parse as force row: "index fx fy fz"
                        int(s[0])
                        float(s[1])
                        forces.append([float(s[1]), float(s[2]), float(s[3])])
                        n_forces += 1
                    except (ValueError, IndexError):
                        # Not a force row, we're done
                        flag = 0
                else:
                    flag = 0
            elif line.startswith("Total energy:"):
                s = line.split()
                energy = float(s[2])
                flag = 0

    symbols = np.array(symbols)
    forces = np.array(forces)
    coord = np.array(coord)
    assert coord.shape == forces.shape

    return symbols, coord, energy, forces
