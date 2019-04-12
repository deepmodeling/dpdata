import os,glob,shutil
import numpy as np
from .raw import load_type

def _load_set(folder) :
    cells = np.load(os.path.join(folder, 'box.npy'))
    coords = np.load(os.path.join(folder, 'coord.npy'))
    eners = np.load(os.path.join(folder, 'energy.npy'))
    forces = np.load(os.path.join(folder, 'force.npy'))
    virs = None
    if os.path.isfile(os.path.join(folder, 'virial.npy')) :
        virs = np.load(os.path.join(folder, 'virial.npy'))
    return cells, coords, eners, forces, virs

def to_system_data(folder, type_map = None) :
    data = load_type(folder, type_map = type_map)
    data['orig'] = np.zeros([3])
    data['virials'] = []
    sets = glob.glob(os.path.join(folder, 'set.*'))
    all_cells = []
    all_coords = []
    all_eners = []
    all_forces = []
    all_virs = []
    for ii in sets :
        cells, coords, eners, forces, virs = _load_set(ii)
        nframes = np.reshape(cells, [-1,3,3]).shape[0]
        all_cells.append(np.reshape(cells, [nframes,3,3]))
        all_coords.append(np.reshape(coords, [nframes,-1,3]))
        all_eners.append(np.reshape(eners, [nframes]))
        all_forces.append(np.reshape(forces, [nframes,-1,3]))
        if len(virs) > 0:
            virs = all_virs.append(np.reshape(virs, [nframes,3,3]))
    data['cells'] = np.concatenate(all_cells, axis = 0)
    data['coords'] = np.concatenate(all_coords, axis = 0)
    data['energies'] = np.concatenate(all_eners, axis = 0)
    data['forces'] = np.concatenate(all_forces, axis = 0)
    if len(all_virs) > 0:
        data['virials'] = np.concatenate(all_virs, axis = 0)
    return data


def dump(folder, 
         data, 
         set_size = 5000, 
         comp_prec = np.float32,
         remove_sets = True) :
    os.makedirs(folder, exist_ok = True)
    sets = glob.glob(os.path.join(folder, 'set.*'))
    if len(sets) > 0:
        if remove_sets :
            for ii in sets :
                shutil.rmtree(ii)
        else :
            raise RuntimeError('found ' + str(sets) + ' in ' + folder + 'not a clean deepmd raw dir. please firstly clean set.* then try compress')
    # dump raw 
    np.savetxt(os.path.join(folder, 'type.raw'), data['atom_types'], fmt = '%d')    
    # reshape frame properties and convert prec
    nframes = data['cells'].shape[0]
    cells  = np.reshape(data['cells'],    [nframes,  9]).astype(comp_prec)
    coords = np.reshape(data['coords'],   [nframes, -1]).astype(comp_prec)
    eners  = np.reshape(data['energies'], [nframes    ]).astype(comp_prec)
    forces = np.reshape(data['forces'],   [nframes, -1]).astype(comp_prec)
    if len(data['virials']) > 0 :
        virials = np.reshape(data['virials'],   [nframes, 9]).astype(comp_prec)
    else :
        virials = []
    # dump frame properties: cell, coord, energy, force and virial
    nsets = nframes // set_size
    if set_size * nsets < nframes :
        nsets += 1
    for ii in range(nsets) :
        set_stt = ii * set_size
        set_end = (ii+1) * set_size
        set_folder = os.path.join(folder, 'set.%03d' % ii)
        os.makedirs(set_folder)
        np.save(os.path.join(set_folder, 'box'),        cells  [set_stt:set_end])
        np.save(os.path.join(set_folder, 'coord'),      coords [set_stt:set_end])
        np.save(os.path.join(set_folder, 'energy'),     eners  [set_stt:set_end])
        np.save(os.path.join(set_folder, 'force'),      forces [set_stt:set_end])
        if len(virials) > 0:
            np.save(os.path.join(set_folder, 'virial'), virials[set_stt:set_end])
        
        
