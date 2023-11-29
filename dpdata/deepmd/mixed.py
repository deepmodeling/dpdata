import glob
import os
import shutil

import numpy as np


def load_type(folder):
    data = {}
    data["atom_names"] = []
    # if find type_map.raw, use it
    assert os.path.isfile(
        os.path.join(folder, "type_map.raw")
    ), "Mixed type system must have type_map.raw!"
    with open(os.path.join(folder, "type_map.raw")) as fp:
        data["atom_names"] = fp.read().split()

    return data


def formula(atom_names, atom_numbs):
    """Return the formula of this system, like C3H5O2."""
    return "".join([f"{symbol}{numb}" for symbol, numb in zip(atom_names, atom_numbs)])


def _cond_load_data(fname):
    tmp = None
    if os.path.isfile(fname):
        tmp = np.load(fname)
    return tmp


def _load_set(folder, nopbc: bool):
    coords = np.load(os.path.join(folder, "coord.npy"))
    if nopbc:
        cells = np.zeros((coords.shape[0], 3, 3))
    else:
        cells = np.load(os.path.join(folder, "box.npy"))
    eners = _cond_load_data(os.path.join(folder, "energy.npy"))
    forces = _cond_load_data(os.path.join(folder, "force.npy"))
    virs = _cond_load_data(os.path.join(folder, "virial.npy"))
    real_atom_types = np.load(os.path.join(folder, "real_atom_types.npy"))
    return cells, coords, eners, forces, virs, real_atom_types


def to_system_data(folder, type_map=None, labels=True):
    # data is empty
    data = load_type(folder)
    old_type_map = data["atom_names"].copy()
    if type_map is not None:
        assert isinstance(type_map, list)
        missing_type = [i for i in old_type_map if i not in type_map]
        assert (
            not missing_type
        ), f"These types are missing in selected type_map: {missing_type} !"
        index_map = np.array([type_map.index(i) for i in old_type_map])
        data["atom_names"] = type_map.copy()
    else:
        index_map = None
    data["orig"] = np.zeros([3])
    if os.path.isfile(os.path.join(folder, "nopbc")):
        data["nopbc"] = True
    sets = sorted(glob.glob(os.path.join(folder, "set.*")))
    all_cells = []
    all_coords = []
    all_eners = []
    all_forces = []
    all_virs = []
    all_real_atom_types = []
    for ii in sets:
        cells, coords, eners, forces, virs, real_atom_types = _load_set(
            ii, data.get("nopbc", False)
        )
        nframes = np.reshape(cells, [-1, 3, 3]).shape[0]
        all_cells.append(np.reshape(cells, [nframes, 3, 3]))
        all_coords.append(np.reshape(coords, [nframes, -1, 3]))
        if index_map is None:
            all_real_atom_types.append(np.reshape(real_atom_types, [nframes, -1]))
        else:
            all_real_atom_types.append(
                np.reshape(index_map[real_atom_types], [nframes, -1])
            )
        if eners is not None:
            eners = np.reshape(eners, [nframes])
        if labels:
            if eners is not None and eners.size > 0:
                all_eners.append(np.reshape(eners, [nframes]))
            if forces is not None and forces.size > 0:
                all_forces.append(np.reshape(forces, [nframes, -1, 3]))
            if virs is not None and virs.size > 0:
                all_virs.append(np.reshape(virs, [nframes, 3, 3]))
    all_cells_concat = np.concatenate(all_cells, axis=0)
    all_coords_concat = np.concatenate(all_coords, axis=0)
    all_real_atom_types_concat = np.concatenate(all_real_atom_types, axis=0)
    all_eners_concat = None
    all_forces_concat = None
    all_virs_concat = None
    if len(all_eners) > 0:
        all_eners_concat = np.concatenate(all_eners, axis=0)
    if len(all_forces) > 0:
        all_forces_concat = np.concatenate(all_forces, axis=0)
    if len(all_virs) > 0:
        all_virs_concat = np.concatenate(all_virs, axis=0)
    data_list = []
    while True:
        if all_real_atom_types_concat.size == 0:
            break
        temp_atom_numbs = [
            np.count_nonzero(all_real_atom_types_concat[0] == i)
            for i in range(len(data["atom_names"]))
        ]
        # temp_formula = formula(data['atom_names'], temp_atom_numbs)
        temp_idx = np.arange(all_real_atom_types_concat.shape[0])[
            (all_real_atom_types_concat == all_real_atom_types_concat[0]).all(-1)
        ]
        rest_idx = np.arange(all_real_atom_types_concat.shape[0])[
            (all_real_atom_types_concat != all_real_atom_types_concat[0]).any(-1)
        ]
        temp_data = data.copy()
        temp_data["atom_names"] = data["atom_names"].copy()
        temp_data["atom_numbs"] = temp_atom_numbs
        temp_data["atom_types"] = all_real_atom_types_concat[0]
        all_real_atom_types_concat = all_real_atom_types_concat[rest_idx]
        temp_data["cells"] = all_cells_concat[temp_idx]
        all_cells_concat = all_cells_concat[rest_idx]
        temp_data["coords"] = all_coords_concat[temp_idx]
        all_coords_concat = all_coords_concat[rest_idx]
        if labels:
            if all_eners_concat is not None and all_eners_concat.size > 0:
                temp_data["energies"] = all_eners_concat[temp_idx]
                all_eners_concat = all_eners_concat[rest_idx]
            if all_forces_concat is not None and all_forces_concat.size > 0:
                temp_data["forces"] = all_forces_concat[temp_idx]
                all_forces_concat = all_forces_concat[rest_idx]
            if all_virs_concat is not None and all_virs_concat.size > 0:
                temp_data["virials"] = all_virs_concat[temp_idx]
                all_virs_concat = all_virs_concat[rest_idx]
        data_list.append(temp_data)
    return data_list


def dump(folder, data, set_size=2000, comp_prec=np.float32, remove_sets=True):
    os.makedirs(folder, exist_ok=True)
    sets = sorted(glob.glob(os.path.join(folder, "set.*")))
    if len(sets) > 0:
        if remove_sets:
            for ii in sets:
                shutil.rmtree(ii)
        else:
            raise RuntimeError(
                "found "
                + str(sets)
                + " in "
                + folder
                + "not a clean deepmd raw dir. please firstly clean set.* then try compress"
            )
    # if not converted to mixed
    if "real_atom_types" not in data:
        from dpdata import LabeledSystem, System

        if "energies" in data:
            temp_sys = LabeledSystem(data=data)
        else:
            temp_sys = System(data=data)
        temp_sys.convert_to_mixed_type()
    # dump raw
    np.savetxt(os.path.join(folder, "type.raw"), data["atom_types"], fmt="%d")
    np.savetxt(os.path.join(folder, "type_map.raw"), data["real_atom_names"], fmt="%s")
    # BondOrder System
    if "bonds" in data:
        np.savetxt(
            os.path.join(folder, "bonds.raw"),
            data["bonds"],
            header="begin_atom, end_atom, bond_order",
        )
    if "formal_charges" in data:
        np.savetxt(os.path.join(folder, "formal_charges.raw"), data["formal_charges"])
    # reshape frame properties and convert prec
    nframes = data["cells"].shape[0]
    cells = np.reshape(data["cells"], [nframes, 9]).astype(comp_prec)
    coords = np.reshape(data["coords"], [nframes, -1]).astype(comp_prec)
    eners = None
    forces = None
    virials = None
    real_atom_types = None
    if "energies" in data:
        eners = np.reshape(data["energies"], [nframes]).astype(comp_prec)
    if "forces" in data:
        forces = np.reshape(data["forces"], [nframes, -1]).astype(comp_prec)
    if "virials" in data:
        virials = np.reshape(data["virials"], [nframes, 9]).astype(comp_prec)
    if "atom_pref" in data:
        atom_pref = np.reshape(data["atom_pref"], [nframes, -1]).astype(comp_prec)
    if "real_atom_types" in data:
        real_atom_types = np.reshape(data["real_atom_types"], [nframes, -1]).astype(
            np.int64
        )
    # dump frame properties: cell, coord, energy, force and virial
    nsets = nframes // set_size
    if set_size * nsets < nframes:
        nsets += 1
    for ii in range(nsets):
        set_stt = ii * set_size
        set_end = (ii + 1) * set_size
        set_folder = os.path.join(folder, "set.%06d" % ii)
        os.makedirs(set_folder)
        np.save(os.path.join(set_folder, "box"), cells[set_stt:set_end])
        np.save(os.path.join(set_folder, "coord"), coords[set_stt:set_end])
        if eners is not None:
            np.save(os.path.join(set_folder, "energy"), eners[set_stt:set_end])
        if forces is not None:
            np.save(os.path.join(set_folder, "force"), forces[set_stt:set_end])
        if virials is not None:
            np.save(os.path.join(set_folder, "virial"), virials[set_stt:set_end])
        if real_atom_types is not None:
            np.save(
                os.path.join(set_folder, "real_atom_types"),
                real_atom_types[set_stt:set_end],
            )
        if "atom_pref" in data:
            np.save(os.path.join(set_folder, "atom_pref"), atom_pref[set_stt:set_end])
    try:
        os.remove(os.path.join(folder, "nopbc"))
    except OSError:
        pass
    if data.get("nopbc", False):
        with open(os.path.join(folder, "nopbc"), "w") as fw_nopbc:
            pass


def mix_system(*system, type_map, **kwargs):
    """Mix the systems into mixed_type ones according to the unified given type_map.

    Parameters
    ----------
    *system : System
        The systems to mix
    type_map : list of str
        Maps atom type to name
    **kwargs : dict
        Other parameters

    Returns
    -------
    mixed_systems: dict
        dict of mixed system with key 'atom_numbs'
    """
    mixed_systems = {}
    temp_systems = {}
    atom_numbs_frame_index = {}  # index of frames in cur sys
    for sys in system:
        tmp_sys = sys.copy()
        natom = tmp_sys.get_natoms()
        tmp_sys.convert_to_mixed_type(type_map=type_map)
        if str(natom) not in atom_numbs_frame_index:
            atom_numbs_frame_index[str(natom)] = 0
        atom_numbs_frame_index[str(natom)] += tmp_sys.get_nframes()
        if str(natom) not in temp_systems or not temp_systems[str(natom)]:
            temp_systems[str(natom)] = tmp_sys
        else:
            temp_systems[str(natom)].append(tmp_sys)
    for natom in temp_systems:
        if atom_numbs_frame_index[natom] > 0:
            sys_name = f"{natom}"
            mixed_systems[sys_name] = temp_systems[natom]
    return mixed_systems


def split_system(sys, split_num=10000):
    rest = sys.get_nframes() - split_num
    if rest <= 0:
        return sys, None, 0
    else:
        split_sys = sys.sub_system(range(split_num))
        rest_sys = sys.sub_system(range(split_num, sys.get_nframes()))
        return split_sys, rest_sys, rest
