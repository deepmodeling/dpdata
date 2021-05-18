import os
import numpy as np

def load_type(folder, type_map = None) :
    data = {}
    data['atom_types'] \
        = np.loadtxt(os.path.join(folder, 'type.raw'), ndmin=1).astype(int)
    ntypes = np.max(data['atom_types']) + 1
    data['atom_numbs'] = []
    for ii in range (ntypes) :
        data['atom_numbs'].append(np.count_nonzero(data['atom_types'] == ii))
    data['atom_names'] = []
    # if find type_map.raw, use it
    if os.path.isfile(os.path.join(folder, 'type_map.raw')) :
        with open(os.path.join(folder, 'type_map.raw')) as fp:
            my_type_map = fp.read().split()
    # else try to use arg type_map 
    elif type_map is not None:
        my_type_map = type_map
    # in the last case, make artificial atom names
    else:
        my_type_map = []
        for ii in range(ntypes) :
            my_type_map.append('Type_%d' % ii)
    assert(len(my_type_map) >= len(data['atom_numbs']))
    for ii in range(len(data['atom_numbs'])) :
        data['atom_names'].append(my_type_map[ii])

    return data


def to_system_data(folder, type_map = None, labels = True) :
    if os.path.isdir(folder) :
        data = load_type(folder, type_map = type_map)
        data['orig'] = np.zeros([3])        
        data['cells'] = np.loadtxt(os.path.join(folder, 'box.raw'))
        data['coords'] = np.loadtxt(os.path.join(folder, 'coord.raw'))
        data['cells'] = np.reshape(data['cells'], [-1, 3, 3])
        nframes = data['cells'].shape[0]
        data['cells'] = np.reshape(data['cells'], [nframes, 3, 3])
        data['coords'] = np.reshape(data['coords'], [nframes, -1, 3])
        if labels :
            if os.path.exists(os.path.join(folder, 'energy.raw')) :        
                data['energies'] = np.loadtxt(os.path.join(folder, 'energy.raw'))
                data['energies'] = np.reshape(data['energies'], [nframes])
            if os.path.exists(os.path.join(folder, 'force.raw')) :
                data['forces'] = np.loadtxt(os.path.join(folder, 'force.raw'))
                data['forces'] = np.reshape(data['forces'], [nframes, -1, 3])
            if os.path.exists(os.path.join(folder, 'virial.raw')) :
                data['virials'] = np.loadtxt(os.path.join(folder, 'virial.raw'))
                data['virials'] = np.reshape(data['virials'], [nframes, 3, 3])
        if os.path.isfile(os.path.join(folder, "nopbc")):
            data['nopbc'] = True
        return data
    else :        
        raise RuntimeError('not dir ' + folder)


def dump (folder, data) :
    os.makedirs(folder, exist_ok = True)
    nframes = data['cells'].shape[0]
    np.savetxt(os.path.join(folder, 'type.raw'),    data['atom_types'], fmt = '%d')
    np.savetxt(os.path.join(folder, 'type_map.raw'),    data['atom_names'], fmt = '%s')
    np.savetxt(os.path.join(folder, 'box.raw'),     np.reshape(data['cells'],    [nframes,  9]))
    np.savetxt(os.path.join(folder, 'coord.raw'),   np.reshape(data['coords'],   [nframes, -1]))
    # BondOrder System
    if "bonds" in data:
        np.savetxt(os.path.join(folder, "bonds.raw"), data['bonds'], header="begin_atom, end_atom, bond_order")
    if "formal_charges" in data:
        np.savetxt(os.path.join(folder, "formal_charges.raw"), data['formal_charges'])
    # Labeled System
    if 'energies' in data :
        np.savetxt(os.path.join(folder, 'energy.raw'),  np.reshape(data['energies'], [nframes,  1]))
    if 'forces' in data :
        np.savetxt(os.path.join(folder, 'force.raw'),   np.reshape(data['forces'],   [nframes, -1]))
    if 'virials' in data :
        np.savetxt(os.path.join(folder, 'virial.raw'), np.reshape(data['virials'], [nframes, 9]))
    try:
        os.remove(os.path.join(folder, "nopbc"))
    except OSError:
        pass
    if data.get("nopbc", False):
        with open(os.path.join(folder, "nopbc"), "w") as fw_nopbc:
            pass

