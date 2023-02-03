from typing import Tuple

import numpy as np


def coord_to_xyz(coord: np.ndarray, types: list) -> str:
    """Convert coordinates and types to xyz format.

    Parameters
    ----------
    coord : np.ndarray
        coordinates, Nx3 array
    types : list
        list of types

    Returns
    -------
    str
        xyz format string

    Examples
    --------
    >>> coord_to_xyz(np.ones((1,3)), ["C"])
    1

    C 1.000000 1.000000 1.000000
    """
    buff = [str(len(types)), ""]
    for at, cc in zip(types, coord):
        buff.append("{} {:.6f} {:.6f} {:.6f}".format(at, *cc))
    return "\n".join(buff)


def xyz_to_coord(xyz: str) -> Tuple[np.ndarray, list]:
    """Convert xyz format to coordinates and types.

    Parameters
    ----------
    xyz : str
        xyz format string

    Returns
    -------
    coords : np.ndarray
        coordinates, Nx3 array
    types : list
        list of types
    """
    symbols = []
    coords = []
    for ii, line in enumerate(xyz.split("\n")):
        if ii == 0:
            natoms = int(line.strip())
        elif 2 <= ii <= 1 + natoms:
            # symbol x y z
            symbol, x, y, z = line.split()
            coords.append((float(x), float(y), float(z)))
            symbols.append(symbol)
    return np.array(coords), symbols
