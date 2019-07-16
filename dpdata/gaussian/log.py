import numpy as np


hartree2ev = 27.211386018
bohr2ang = 0.52917721067

length_convert = bohr2ang
energy_convert = hartree2ev
force_convert = energy_convert / length_convert

symbols = ['X', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

def to_system_data(file_name):
    data = {}
    # read from log lines
    flag = 0
    energy = 0
    coords = []
    atom_symbols = []
    forces = []

    with open(file_name) as fp:
        for line in fp:
            if line.startswith(" SCF Done"):
                # energies
                energy = float(line.split()[4])
            elif line.startswith(" Center     Atomic                   Forces (Hartrees/Bohr)"):
                flag = 1
                forces = []
            elif line.startswith("                          Input orientation:"):
                flag = 5
                coords = []
                atom_symbols = []

            if 1 <= flag <= 3 or 5 <= flag <= 9:
                flag += 1
            elif flag == 4:
                # forces
                if line.startswith(" -------"):
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

    assert(coords), "cannot find coords"
    assert(energy), "cannot find energies"
    assert(forces), "cannot find forces"

    atom_names, data['atom_types'], atom_numbs = np.unique(atom_symbols, return_inverse=True, return_counts=True)
    data['atom_names'] = list(atom_names)
    data['atom_numbs'] = list(atom_numbs)
    data['forces'] = np.array([forces]) * force_convert
    data['energies'] = np.array([energy]) * energy_convert
    data['coords'] = np.array([coords])
    data['orig'] = np.array([0, 0, 0])
    return data
