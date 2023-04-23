from typing import Tuple

import numpy as np


def read_psi4_output(fn: str) -> Tuple[str, np.ndarray, float, np.ndarray]:
    """Read from Psi4 output.

    Note that both the energy and the gradient should be printed.

    Parameters
    ----------
    fn : str
        file name

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
    with open(fn) as f:
        flag = 0
        for line in f:
            if flag in (1, 3, 4, 5, 6):
                flag += 1
            elif flag == 2:
                s = line.split()
                if not len(s):
                    flag = 0
                else:
                    symbols.append(s[0].capitalize())
                    coord.append([float(s[1]), float(s[2]), float(s[3])])
            elif flag == 7:
                s = line.split()
                if not len(s):
                    flag = 0
                else:
                    forces.append([float(s[1]), float(s[2]), float(s[3])])
            elif line.startswith(
                "       Center              X                  Y                   Z               Mass"
            ):
                # coord
                flag = 1
                coord = []
                symbols = []
            elif line.startswith("  ## Total Gradient"):
                flag = 3
                forces = []
            elif line.startswith("    Total Energy ="):
                energy = float(line.split()[-1])
    symbols = np.array(symbols)
    forces = -np.array(forces)
    coord = np.array(coord)
    assert coord.shape == forces.shape

    return symbols, coord, energy, forces
