from __future__ import annotations

import copy

import numpy as np

from .comp import dump as comp_dump
from .comp import to_system_data as comp_to_system_data


def to_system_data(folder, type_map=None, labels=True):
    data = comp_to_system_data(folder, type_map, labels)
    # data is empty
    old_type_map = data["atom_names"].copy()
    if type_map is not None:
        assert isinstance(type_map, list)
        missing_type = [i for i in old_type_map if i not in type_map]
        assert not missing_type, (
            f"These types are missing in selected type_map: {missing_type} !"
        )
        index_map = np.array([type_map.index(i) for i in old_type_map])
        data["atom_names"] = type_map.copy()
    else:
        index_map = None
    all_real_atom_types_concat = data.pop("real_atom_types").astype(int)
    if index_map is not None:
        all_real_atom_types_concat = index_map[all_real_atom_types_concat]
    all_cells_concat = data["cells"]
    all_coords_concat = data["coords"]
    if labels:
        all_eners_concat = data.get("energies")
        all_forces_concat = data.get("forces")
        all_virs_concat = data.get("virials")

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
    # if not converted to mixed
    if "real_atom_types" not in data:
        from dpdata import LabeledSystem, System

        # not change the original content
        data = copy.deepcopy(data)

        if "energies" in data:
            temp_sys = LabeledSystem(data=data)
        else:
            temp_sys = System(data=data)
        temp_sys.convert_to_mixed_type()

    data = data.copy()
    data["atom_names"] = data.pop("real_atom_names")
    comp_dump(folder, data, set_size, comp_prec, remove_sets)


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
