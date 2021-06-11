import numpy as np


START = 0
READ_CHARGE = 1
READ_CHARGE_SUCCESS = 2
READ_COORDS_START = 3
READ_COORDS = 6

def parse_sqm_out(fname):
    '''
        Read atom symbols, charges and coordinates from ambertools sqm.out file
    '''
    atom_symbols = []
    coords = []
    charges = []
    with open(fname) as f:
        flag = 0
        for line in f:
            if line.startswith("  Atom    Element       Mulliken Charge"):
                flag = READ_CHARGE
            elif line.startswith(" Total Mulliken Charge"):
                flag = READ_CHARGE_SUCCESS
            elif line.startswith(" Final Structure"):
                flag = READ_COORDS_START
            elif flag == READ_CHARGE:
                ls = line.strip().split()
                atom_symbols.append(ls[-2])
                charges.append(float(ls[-1]))
            elif READ_COORDS_START <= flag < READ_COORDS:
                flag += 1
            elif flag == READ_COORDS:
                ls = line.strip()
                if not ls:
                    break
                else:
                    symbol = line.strip().split()[-4]
                    coord = list(map(float, line.strip().split()[-3:]))
                    coords.append(coord)
    return atom_symbols, charges, np.array(coords)


def to_system_data(fname):
    data = {}
    atom_symbols, charges, coords = parse_sqm_out(fname)
    atom_names, data['atom_types'], atom_numbs = np.unique(atom_symbols, return_inverse=True, return_counts=True)
    data['atom_names'] = list(atom_names)
    data['atom_numbs'] = list(atom_numbs)
    data['charges'] = np.array([charges])
    data['coords'] = np.array([coords])
    data['orig'] = np.array([0, 0, 0])
    data['cells'] = np.array([[[100., 0., 0.], [0., 100., 0.], [0., 0., 100.]]])
    data['nopbc'] = True
    return data
