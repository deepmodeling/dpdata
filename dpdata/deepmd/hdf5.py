"""Utils for deepmd/hdf5 format."""
from typing import Union

import h5py
import numpy as np

from wcmatch.glob import globfilter


__all__ = ['to_system_data', 'dump']

def to_system_data(f: Union[h5py.File, h5py.Group],
                   folder: str, 
                   type_map: list = None,
                   labels: bool = True) :
    """Load a HDF5 file.

    Parameters
    ----------
    f : h5py.File or h5py.Group
        HDF5 file or group object
    folder : str
        path in the HDF5 file
    type_map : list
        type map
    labels : bool
        labels
    """
    g = f[folder] if folder else f

    data = {}
    data['atom_types'] = g['type.raw'][:]
    ntypes = np.max(data['atom_types']) + 1
    natoms = data['atom_types'].size
    data['atom_numbs'] = []
    for ii in range (ntypes) :
        data['atom_numbs'].append(np.count_nonzero(data['atom_types'] == ii))
    data['atom_names'] = []
    # if find type_map.raw, use it
    if 'type_map.raw' in g.keys():
        my_type_map = list(np.char.decode(g['type_map.raw'][:]))
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

    data['orig'] = np.zeros([3])
    if 'nopbc' in g.keys():
        data['nopbc'] = True
    sets = globfilter(g.keys(), 'set.*')

    data_types = {
        'cells': {'fn': 'box', 'labeled': False, 'shape': (3,3), 'required': 'nopbc' not in data},
        'coords': {'fn': 'coord', 'labeled': False, 'shape': (natoms,3), 'required': True},
        'energies': {'fn': 'energy', 'labeled': True, 'shape': tuple(), 'required': False},
        'forces': {'fn': 'force', 'labeled': True, 'shape': (natoms,3), 'required': False},
        'virials': {'fn': 'virial', 'labeled': True, 'shape': (3,3), 'required': False},
    }
    
    for dt, prop in data_types.items():
        all_data = []

        for ii in sets:
            set = g[ii]
            fn = '%s.npy' % prop['fn']
            if fn in set.keys():
                dd = set[fn][:]
                nframes = dd.shape[0]
                all_data.append(np.reshape(dd, (nframes, *prop['shape'])))
            elif prop['required']:
                raise RuntimeError("%s/%s/%s not found" % (folder, ii, fn))
        
        if len(all_data) > 0 :
            data[dt] = np.concatenate(all_data, axis = 0)
    if 'cells' not in data:
        nframes = data['coords'].shape[0]
        data['cells'] = np.zeros((nframes, 3, 3))
    return data

def dump(f: Union[h5py.File, h5py.Group],
         folder: str,
         data: dict,
         set_size = 5000,
         comp_prec = np.float32,
         ) -> None:
    """Dump data to a HDF5 file.

    Parameters
    ----------
    f : h5py.File or h5py.Group
        HDF5 file or group object
    folder : str
        path in the HDF5 file
    data : dict
        System or LabeledSystem data
    set_size : int, default: 5000
        size of a set
    comp_prec : np.dtype, default: np.float32
        precision of data
    """
    # if folder is None, use the root of the file
    if folder:
        if folder in f:
            del f[folder]
        g = f.create_group(folder)
    else:
        g = f
    # dump raw (array in fact)
    g.create_dataset('type.raw', data=data['atom_types'])    
    g.create_dataset('type_map.raw', data=np.array(data['atom_names'], dtype='S'))
    # BondOrder System
    if "bonds" in data:
        g.create_dataset("bonds.raw", data=data['bonds'])
    if "formal_charges" in data:
        g.create_dataset("formal_charges.raw", data=data['formal_charges'])
    # reshape frame properties and convert prec
    nframes = data['cells'].shape[0]

    nopbc = data.get("nopbc", False)
    reshaped_data = {}

    data_types = {
        'cells': {'fn': 'box', 'shape': (nframes, 9), 'dump': not nopbc},
        'coords': {'fn': 'coord',  'shape': (nframes, -1), 'dump': True},
        'energies': {'fn': 'energy', 'shape': (nframes,), 'dump': True},
        'forces': {'fn': 'force', 'shape': (nframes, -1), 'dump': True},
        'virials': {'fn': 'virial', 'shape': (nframes, 9), 'dump': True},
    }
    for dt, prop in data_types.items():
        if dt in data:
            if prop['dump']:
                reshaped_data[dt] = np.reshape(data[dt], prop['shape']).astype(comp_prec)

    # dump frame properties: cell, coord, energy, force and virial
    nsets = nframes // set_size
    if set_size * nsets < nframes :
        nsets += 1
    for ii in range(nsets) :
        set_stt = ii * set_size
        set_end = (ii+1) * set_size
        set_folder = g.create_group('set.%03d' % ii)
        for dt, prop in data_types.items():
            if dt in reshaped_data:
                set_folder.create_dataset('%s.npy' % prop['fn'], data=reshaped_data[dt][set_stt:set_end])

    if nopbc:
       g.create_dataset("nopbc", data=True)
