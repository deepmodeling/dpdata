from __future__ import annotations

import copy
import math

import numpy as np

import dpdata
from dpdata.data_type import Axis

from .comp import dump as comp_dump
from .comp import to_system_data as comp_to_system_data


def _pad_to(sys_data, target_natoms, dtypes):
    """Pad system data dict so that NATOMS dimension becomes target_natoms.

    Virtual atoms get real_atom_types = -1, and all other per-atom data is
    padded with zeros.

    Parameters
    ----------
    sys_data : dict
        System data dict, already in mixed-type format.
    target_natoms : int
        Target number of atoms after padding.
    dtypes : tuple[DataType, ...]
        Registered data types to iterate for generic per-atom padding.
    """
    natoms = sys_data["atom_types"].shape[0]
    npad = target_natoms - natoms
    if npad <= 0:
        return
    nframes = sys_data["coords"].shape[0]

    # Pad atom_types (all MIXED_TOKEN = 0)
    sys_data["atom_types"] = np.concatenate(
        [sys_data["atom_types"], np.zeros(npad, dtype=int)]
    )
    sys_data["atom_numbs"] = [target_natoms]

    # Pad real_atom_types with -1 (virtual atom sentinel)
    sys_data["real_atom_types"] = np.concatenate(
        [
            sys_data["real_atom_types"],
            -np.ones((nframes, npad), dtype=sys_data["real_atom_types"].dtype),
        ],
        axis=1,
    )

    # Pad coords and all other per-atom data generically
    reserved = {
        "atom_numbs",
        "atom_names",
        "atom_types",
        "orig",
        "cells",
        "real_atom_names",
        "real_atom_types",
        "nopbc",
    }
    for dtype in dtypes:
        if dtype.name in reserved:
            continue
        if dtype.name not in sys_data:
            continue
        if not (
            len(dtype.shape) >= 2
            and dtype.shape[0] == Axis.NFRAMES
            and Axis.NATOMS in dtype.shape
        ):
            continue
        axis_natoms = list(dtype.shape).index(Axis.NATOMS)
        arr = sys_data[dtype.name]
        pad_width = [(0, 0)] * len(arr.shape)
        pad_width[axis_natoms] = (0, npad)
        sys_data[dtype.name] = np.pad(
            arr, pad_width, mode="constant", constant_values=0
        )


def _strip_virtual_atoms(atom_types_row, coords, extra_data, dtypes):
    """Strip virtual atoms (type -1) from a group of frames.

    Parameters
    ----------
    atom_types_row : np.ndarray
        1-D array of atom type indices for the group (same for all frames).
    coords : np.ndarray
        Coordinates array, shape (nframes, natoms_padded, 3).
    extra_data : dict
        Dict of {name: array} for this group, arrays already frame-sliced.
    dtypes : tuple[DataType, ...]
        Registered data types.

    Returns
    -------
    atom_types : np.ndarray
        Atom types with virtual atoms removed.
    coords : np.ndarray
        Coords with virtual atoms removed.
    extra_data : dict
        Extra data with virtual atoms removed.
    """
    real_mask = atom_types_row >= 0
    if real_mask.all():
        return atom_types_row, coords, extra_data

    atom_types = atom_types_row[real_mask]
    coords = coords[:, real_mask, :]

    stripped = {}
    for name, arr in extra_data.items():
        for dtype in dtypes:
            if dtype.name == name and Axis.NATOMS in dtype.shape:
                axis_natoms = list(dtype.shape).index(Axis.NATOMS)
                idx = [slice(None)] * len(arr.shape)
                idx[axis_natoms] = real_mask
                arr = arr[tuple(idx)]
                break
        stripped[name] = arr

    return atom_types, coords, stripped


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
        # Preserve -1 (virtual atom sentinel) during remapping
        valid = all_real_atom_types_concat >= 0
        remapped = np.full_like(all_real_atom_types_concat, -1)
        remapped[valid] = index_map[all_real_atom_types_concat[valid]]
        all_real_atom_types_concat = remapped
    all_cells_concat = data["cells"]
    all_coords_concat = data["coords"]

    # handle custom registered data types
    if labels:
        dtypes = dpdata.system.LabeledSystem.DTYPES
    else:
        dtypes = dpdata.system.System.DTYPES
    reserved = {
        "atom_numbs",
        "atom_names",
        "atom_types",
        "real_atom_names",
        "real_atom_types",
        "cells",
        "coords",
        "orig",
        "nopbc",
    }
    extra_data = {}
    for dtype in dtypes:
        name = dtype.name
        if name in reserved:
            continue
        if not (len(dtype.shape) and dtype.shape[0] == dpdata.system.Axis.NFRAMES):
            continue
        if name in data:
            extra_data[name] = data.pop(name)

    data_list = []
    while True:
        if all_real_atom_types_concat.size == 0:
            break
        # temp_formula = formula(data['atom_names'], temp_atom_numbs)
        temp_idx = np.arange(all_real_atom_types_concat.shape[0])[
            (all_real_atom_types_concat == all_real_atom_types_concat[0]).all(-1)
        ]
        rest_idx = np.arange(all_real_atom_types_concat.shape[0])[
            (all_real_atom_types_concat != all_real_atom_types_concat[0]).any(-1)
        ]

        # Extract data for this group
        group_atom_types = all_real_atom_types_concat[0]
        group_coords = all_coords_concat[temp_idx]
        group_extra = {}
        for name in extra_data:
            group_extra[name] = extra_data[name][temp_idx]
            extra_data[name] = extra_data[name][rest_idx]

        # Strip virtual atoms (type -1) introduced by padding
        group_atom_types, group_coords, group_extra = _strip_virtual_atoms(
            group_atom_types, group_coords, group_extra, dtypes
        )

        temp_atom_numbs = [
            np.count_nonzero(group_atom_types == i)
            for i in range(len(data["atom_names"]))
        ]

        temp_data = data.copy()
        temp_data["atom_names"] = data["atom_names"].copy()
        temp_data["atom_numbs"] = temp_atom_numbs
        temp_data["atom_types"] = group_atom_types
        all_real_atom_types_concat = all_real_atom_types_concat[rest_idx]
        temp_data["cells"] = all_cells_concat[temp_idx]
        all_cells_concat = all_cells_concat[rest_idx]
        temp_data["coords"] = group_coords
        all_coords_concat = all_coords_concat[rest_idx]

        for name in group_extra:
            temp_data[name] = group_extra[name]

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


def mix_system(*system, type_map, atom_numb_pad=None, **kwargs):
    """Mix the systems into mixed_type ones according to the unified given type_map.

    Parameters
    ----------
    *system : System
        The systems to mix
    type_map : list of str
        Maps atom type to name
    atom_numb_pad : int, optional
        If provided, pad atom counts to the next multiple of this number
        using virtual atoms (type -1 in real_atom_types). This reduces the
        number of subdirectories when systems have many different atom counts.
        For example, atom_numb_pad=8 groups systems into multiples of 8.
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
    # Use LabeledSystem DTYPES as superset for generic per-atom padding
    dtypes = dpdata.system.LabeledSystem.DTYPES
    for sys in system:
        tmp_sys = sys.copy()
        natom = tmp_sys.get_natoms()
        tmp_sys.convert_to_mixed_type(type_map=type_map)
        if atom_numb_pad is not None and atom_numb_pad > 1:
            padded_natom = math.ceil(natom / atom_numb_pad) * atom_numb_pad
            _pad_to(tmp_sys.data, padded_natom, dtypes)
            group_key = str(padded_natom)
        else:
            group_key = str(natom)
        if group_key not in atom_numbs_frame_index:
            atom_numbs_frame_index[group_key] = 0
        atom_numbs_frame_index[group_key] += tmp_sys.get_nframes()
        if group_key not in temp_systems or not temp_systems[group_key]:
            temp_systems[group_key] = tmp_sys
        else:
            temp_systems[group_key].append(tmp_sys)
    for natom_key in temp_systems:
        if atom_numbs_frame_index[natom_key] > 0:
            mixed_systems[natom_key] = temp_systems[natom_key]
    return mixed_systems


def split_system(sys, split_num=10000):
    rest = sys.get_nframes() - split_num
    if rest <= 0:
        return sys, None, 0
    else:
        split_sys = sys.sub_system(range(split_num))
        rest_sys = sys.sub_system(range(split_num, sys.get_nframes()))
        return split_sys, rest_sys, rest
