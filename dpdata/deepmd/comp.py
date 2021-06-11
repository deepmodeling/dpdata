import os,glob,shutil
import numpy as np
from .raw import load_type

def _cond_load_data(fname) :
    tmp = None
    if os.path.isfile(fname) :
        tmp = np.load(fname)
    return tmp

def _load_set(folder) :
    cells = np.load(os.path.join(folder, 'box.npy'))
    coords = np.load(os.path.join(folder, 'coord.npy'))
    eners  = _cond_load_data(os.path.join(folder, 'energy.npy'))
    forces = _cond_load_data(os.path.join(folder, 'force.npy'))
    virs   = _cond_load_data(os.path.join(folder, 'virial.npy'))
    return cells, coords, eners, forces, virs

def to_system_data(folder, 
                   type_map = None, 
                   labels = True) :
    # data is empty
    data = load_type(folder, type_map = type_map)
    data['orig'] = np.zeros([3])
    sets = sorted(glob.glob(os.path.join(folder, 'set.*')))
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
        if eners is not None:
            eners = np.reshape(eners, [nframes])
        if labels:
            if eners is not None and eners.size > 0:
                all_eners.append(np.reshape(eners, [nframes]))
            if forces is not None and forces.size > 0:
                all_forces.append(np.reshape(forces, [nframes,-1,3]))
            if virs is not None and virs.size > 0:
                all_virs.append(np.reshape(virs, [nframes,3,3]))
    data['cells'] = np.concatenate(all_cells, axis = 0)
    data['coords'] = np.concatenate(all_coords, axis = 0)    
    if len(all_eners) > 0 :
        data['energies'] = np.concatenate(all_eners, axis = 0)
    if len(all_forces) > 0 :
        data['forces'] = np.concatenate(all_forces, axis = 0)
    if len(all_virs) > 0:
        data['virials'] = np.concatenate(all_virs, axis = 0)
    if os.path.isfile(os.path.join(folder, "nopbc")):
        data['nopbc'] = True
    return data


def dump(folder, 
         data, 
         set_size = 5000, 
         comp_prec = np.float32,
         remove_sets = True) :
    os.makedirs(folder, exist_ok = True)
    sets = sorted(glob.glob(os.path.join(folder, 'set.*')))
    if len(sets) > 0:
        if remove_sets :
            for ii in sets :
                shutil.rmtree(ii)
        else :
            raise RuntimeError('found ' + str(sets) + ' in ' + folder + 'not a clean deepmd raw dir. please firstly clean set.* then try compress')
    # dump raw 
    np.savetxt(os.path.join(folder, 'type.raw'), data['atom_types'], fmt = '%d')    
    np.savetxt(os.path.join(folder, 'type_map.raw'),    data['atom_names'], fmt = '%s')
    # BondOrder System
    if "bonds" in data:
        np.savetxt(os.path.join(folder, "bonds.raw"), data['bonds'], header="begin_atom, end_atom, bond_order")
    if "formal_charges" in data:
        np.savetxt(os.path.join(folder, "formal_charges.raw"), data['formal_charges'])
    # reshape frame properties and convert prec
    nframes = data['cells'].shape[0]
    cells  = np.reshape(data['cells'],    [nframes,  9]).astype(comp_prec)
    coords = np.reshape(data['coords'],   [nframes, -1]).astype(comp_prec)
    eners = None
    forces = None
    virials = None
    if 'energies' in data:
        eners   = np.reshape(data['energies'], [nframes    ]).astype(comp_prec)    
    if 'forces' in data:
        forces  = np.reshape(data['forces'],   [nframes, -1]).astype(comp_prec)
    if 'virials' in data :
        virials = np.reshape(data['virials'],  [nframes,  9]).astype(comp_prec)
    if 'atom_pref' in data:
        atom_pref = np.reshape(data['atom_pref'], [nframes, -1]).astype(comp_prec)
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
        if eners is not None:
            np.save(os.path.join(set_folder, 'energy'), eners  [set_stt:set_end])
        if forces is not None:
            np.save(os.path.join(set_folder, 'force'),  forces [set_stt:set_end])
        if virials is not None:
            np.save(os.path.join(set_folder, 'virial'), virials[set_stt:set_end])
        if 'atom_pref' in data:
            np.save(os.path.join(set_folder, "atom_pref"), atom_pref[set_stt:set_end])
    try:
        os.remove(os.path.join(folder, "nopbc"))
    except OSError:
        pass
    if data.get("nopbc", False):
        with open(os.path.join(folder, "nopbc"), "w") as fw_nopbc:
            pass
        
