import numpy as np
from dpdata.constant import hartree2ev, bohr2ang, symbols

length_convert = bohr2ang
energy_convert = hartree2ev
force_convert = energy_convert / length_convert

def to_system_data(file_name, md=False):
    data = {}
    # read from log lines
    flag = 0
    energy_t = []
    coords_t = []
    atom_symbols = []
    forces_t = []

    with open(file_name) as fp:
        for line in fp:
            if line.startswith(" SCF Done"):
                # energies
                energy = float(line.split()[4])
            elif line.startswith(" Center     Atomic                   Forces (Hartrees/Bohr)"):
                flag = 1
                forces = []
            elif line.startswith("                          Input orientation:") or line.startswith("                         Z-Matrix orientation:"):
                flag = 5
                coords = []
                atom_symbols = []

            if 1 <= flag <= 3 or 5 <= flag <= 9:
                flag += 1
            elif flag == 4:
                # forces
                if line.startswith(" -------"):
                    forces_t.append(forces)
                    energy_t.append(energy)
                    coords_t.append(coords)
                    flag = 0
                else:
                    s = line.split()
                    forces.append([float(x) for x in s[2:5]])
            elif flag == 10:
                # atom_symbols and coords
                if line.startswith(" -------"):
                    flag = 0
                else:
                    s = line.split()
                    coords.append([float(x) for x in s[3:6]])
                    atom_symbols.append(symbols[int(s[1])])

    assert(coords_t), "cannot find coords"
    assert(energy_t), "cannot find energies"
    assert(forces_t), "cannot find forces"

    atom_names, data['atom_types'], atom_numbs = np.unique(atom_symbols, return_inverse=True, return_counts=True)
    data['atom_names'] = list(atom_names)
    data['atom_numbs'] = list(atom_numbs)
    if not md:
        forces_t = forces_t[-1:]
        energy_t = energy_t[-1:]
        coords_t = coords_t[-1:]
    data['forces'] = np.array(forces_t) * force_convert
    data['energies'] = np.array(energy_t) * energy_convert
    data['coords'] = np.array(coords_t)
    data['orig'] = np.array([0, 0, 0])
    data['cells'] = np.array([[[100., 0., 0.], [0., 100., 0.], [0., 0., 100.]] for _ in energy_t])
    data['nopbc'] = True
    return data
