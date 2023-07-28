from typing import Tuple

import numpy as np


def read_dftb_plus(fn: str) -> Tuple[str, np.ndarray, float, np.ndarray]:
    with open(fn) as f:
        flag = 0
        for line in f:
            if flag == 1:
                flag += 1
            elif flag == 2:
                components = line.split()
                flag += 1
            elif line.startswith("Geometry"):
                flag = 1
                coord = []
                symbols = []
            elif flag in (3, 4, 5, 6):
                s = line.split()
                components_num = int(s[1])
                symbols.append(components[components_num - 1])
                coord.append([float(s[2]), float(s[3]), float(s[4])])
                flag += 1
                if flag == 7:
                    flag = 0
            elif line.startswith("Total Forces"):
                flag = 8
                forces = []
            elif flag in (8, 9, 10, 11):
                s = line.split()
                forces.append([float(s[1]), float(s[2]), float(s[3])])
                flag += 1
                if flag == 12:
                    flag = 0
            elif line.startswith("Total energy:"):
                s = line.split()
                energy = float(s[2])

    symbols = np.array(symbols)
    forces = -np.array(forces)
    coord = np.array(coord)
    assert coord.shape == forces.shape

    return symbols, coord, energy, forces
