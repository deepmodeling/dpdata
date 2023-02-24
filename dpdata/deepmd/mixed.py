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
    """
    Return the formula of this system, like C3H5O2
    """
    return "".join(
        ["{}{}".format(symbol, numb) for symbol, numb in zip(atom_names, atom_numbs)]
    )


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
    data["orig"] = np.zeros([3])
    if os.path.isfile(os.path.join(folder, "nopbc")):
        data["nopbc"] = True
    sets = sorted(glob.glob(os.path.join(folder, "set.*")))
    assert len(sets) == 1, "Mixed type must have only one set!"
    cells, coords, eners, forces, virs, real_atom_types = _load_set(
        sets[0], data.get("nopbc", False)
    )
    nframes = np.reshape(cells, [-1, 3, 3]).shape[0]
    cells = np.reshape(cells, [nframes, 3, 3])
    coords = np.reshape(coords, [nframes, -1, 3])
    real_atom_types = np.reshape(real_atom_types, [nframes, -1])
    natom = real_atom_types.shape[1]
    if labels:
        if eners is not None and eners.size > 0:
            eners = np.reshape(eners, [nframes])
        if forces is not None and forces.size > 0:
            forces = np.reshape(forces, [nframes, -1, 3])
        if virs is not None and virs.size > 0:
            virs = np.reshape(virs, [nframes, 3, 3])
    data_list = []
    while True:
        if real_atom_types.size == 0:
            break
        temp_atom_numbs = [
            np.count_nonzero(real_atom_types[0] == i)
            for i in range(len(data["atom_names"]))
        ]
        # temp_formula = formula(data['atom_names'], temp_atom_numbs)
        temp_idx = np.arange(real_atom_types.shape[0])[
            (real_atom_types == real_atom_types[0]).all(-1)
        ]
        rest_idx = np.arange(real_atom_types.shape[0])[
            (real_atom_types != real_atom_types[0]).any(-1)
        ]
        temp_data = data.copy()
        temp_data["atom_names"] = data["atom_names"].copy()
        temp_data["atom_numbs"] = temp_atom_numbs
        temp_data["atom_types"] = real_atom_types[0]
        real_atom_types = real_atom_types[rest_idx]
        temp_data["cells"] = cells[temp_idx]
        cells = cells[rest_idx]
        temp_data["coords"] = coords[temp_idx]
        coords = coords[rest_idx]
        if labels:
            if eners is not None and eners.size > 0:
                temp_data["energies"] = eners[temp_idx]
                eners = eners[rest_idx]
            if forces is not None and forces.size > 0:
                temp_data["forces"] = forces[temp_idx]
                forces = forces[rest_idx]
            if virs is not None and virs.size > 0:
                temp_data["virials"] = virs[temp_idx]
                virs = virs[rest_idx]
        data_list.append(temp_data)
    return data_list


def dump(folder, data, comp_prec=np.float32, remove_sets=True):
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
    set_folder = os.path.join(folder, "set.%03d" % 0)
    os.makedirs(set_folder)
    np.save(os.path.join(set_folder, "box"), cells)
    np.save(os.path.join(set_folder, "coord"), coords)
    if eners is not None:
        np.save(os.path.join(set_folder, "energy"), eners)
    if forces is not None:
        np.save(os.path.join(set_folder, "force"), forces)
    if virials is not None:
        np.save(os.path.join(set_folder, "virial"), virials)
    if real_atom_types is not None:
        np.save(os.path.join(set_folder, "real_atom_types"), real_atom_types)
    if "atom_pref" in data:
        np.save(os.path.join(set_folder, "atom_pref"), atom_pref)
    try:
        os.remove(os.path.join(folder, "nopbc"))
    except OSError:
        pass
    if data.get("nopbc", False):
        with open(os.path.join(folder, "nopbc"), "w") as fw_nopbc:
            pass


def mix_system(*system, type_map, split_num=200, **kwargs):
    """Mix the systems into mixed_type ones

    Parameters
    ----------
    *system : System
        The systems to mix
    type_map : list of str
        Maps atom type to name
    split_num : int
        Number of frames in each system

    Returns
    -------
    mixed_systems: dict
        dict of mixed system with key '{atom_numbs}/sys.xxx'
    """
    mixed_systems = {}
    temp_systems = {}
    atom_numbs_sys_index = {}  # index of sys
    atom_numbs_frame_index = {}  # index of frames in cur sys
    for sys in system:
        tmp_sys = sys.copy()
        natom = tmp_sys.get_natoms()
        tmp_sys.convert_to_mixed_type(type_map=type_map)
        if str(natom) not in atom_numbs_sys_index:
            atom_numbs_sys_index[str(natom)] = 0
        if str(natom) not in atom_numbs_frame_index:
            atom_numbs_frame_index[str(natom)] = 0
        atom_numbs_frame_index[str(natom)] += tmp_sys.get_nframes()
        if str(natom) not in temp_systems or not temp_systems[str(natom)]:
            temp_systems[str(natom)] = tmp_sys
        else:
            temp_systems[str(natom)].append(tmp_sys)
        if atom_numbs_frame_index[str(natom)] >= split_num:
            while True:
                sys_split, temp_systems[str(natom)], rest_num = split_system(
                    temp_systems[str(natom)], split_num=split_num
                )
                sys_name = (
                    f"{str(natom)}/sys." + "%.6d" % atom_numbs_sys_index[str(natom)]
                )
                mixed_systems[sys_name] = sys_split
                atom_numbs_sys_index[str(natom)] += 1
                if rest_num < split_num:
                    atom_numbs_frame_index[str(natom)] = rest_num
                    break
    for natom in temp_systems:
        if atom_numbs_frame_index[natom] > 0:
            sys_name = f"{natom}/sys." + "%.6d" % atom_numbs_sys_index[natom]
            mixed_systems[sys_name] = temp_systems[natom]
    return mixed_systems


def split_system(sys, split_num=100):
    rest = sys.get_nframes() - split_num
    if rest <= 0:
        return sys, None, 0
    else:
        split_sys = sys.sub_system(range(split_num))
        rest_sys = sys.sub_system(range(split_num, sys.get_nframes()))
        return split_sys, rest_sys, rest
