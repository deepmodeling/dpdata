import os
import numpy as np

def load_type(folder, type_map = None) :
    data = {}
    data['atom_types'] \
        = np.loadtxt(os.path.join(folder, 'type.raw')).astype(int)
    ntypes = np.max(data['atom_types']) + 1
    data['atom_numbs'] = []
    for ii in range (ntypes) :
        data['atom_numbs'].append(np.count_nonzero(data['atom_types'] == ii))
    data['atom_names'] = []
    if type_map == None :
        for ii in range(ntypes) :
            data['atom_names'].append('Type_%d' % ii)
    else :
        assert(len(type_map) >= len(data['atom_numbs']))
        for ii in range(len(data['atom_numbs'])) :
            data['atom_names'].append(type_map[ii])
    return data


def to_system_data(folder, type_map = None) :
    if os.path.isdir(folder) :
        data = load_type(folder, type_map = type_map)
        data['orig'] = np.zeros([3])        
        data['virials'] = []
        data['cells'] = np.loadtxt(os.path.join(folder, 'box.raw'))
        data['coords'] = np.loadtxt(os.path.join(folder, 'coord.raw'))
        data['energies'] = np.loadtxt(os.path.join(folder, 'energy.raw'))
        data['forces'] = np.loadtxt(os.path.join(folder, 'force.raw'))
        data['cells'] = np.reshape(data['cells'], [-1, 3, 3])
        nframes = data['cells'].shape[0]
        data['cells'] = np.reshape(data['cells'], [nframes, 3, 3])
        data['coords'] = np.reshape(data['coords'], [nframes, -1, 3])
        data['energies'] = np.reshape(data['energies'], [nframes])
        data['forces'] = np.reshape(data['forces'], [nframes, -1, 3])
        if os.path.exists(os.path.join(folder, 'virial.raw')) :
            data['virials'] = np.loadtxt(os.path.join(folder, 'virial.raw'))
            data['virials'] = np.reshape(data['virials'], [nframes, 3, 3])
        return data
    else :        
        raise RuntimeError('not dir ' + folder)


def dump (folder, data) :
    os.makedirs(folder, exist_ok = True)
    nframes = data['cells'].shape[0]
    np.savetxt(os.path.join(folder, 'type.raw'),    data['atom_types'], fmt = '%d')
    np.savetxt(os.path.join(folder, 'box.raw'),     np.reshape(data['cells'],    [nframes,  9]))
    np.savetxt(os.path.join(folder, 'coord.raw'),   np.reshape(data['coords'],   [nframes, -1]))
    np.savetxt(os.path.join(folder, 'energy.raw'),  np.reshape(data['energies'], [nframes,  1]))
    np.savetxt(os.path.join(folder, 'force.raw'),   np.reshape(data['forces'],   [nframes, -1]))
    if len(data['virials']) != 0 :            
        np.savetxt(os.path.join(folder, 'virial.raw'), np.reshape(data['virials'], [nframes, 9]))
    

