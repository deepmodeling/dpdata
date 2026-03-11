from __future__ import annotations

import glob
import os
import shutil
import warnings

import numpy as np

import dpdata
from dpdata.utils import open_file

from .raw import load_type


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
    return cells, coords


def to_system_data(folder, type_map=None, labels=True):
    # data is empty
    data = load_type(folder, type_map=type_map)
    data["orig"] = np.zeros([3])
    if os.path.isfile(os.path.join(folder, "nopbc")):
        data["nopbc"] = True
    sets = sorted(glob.glob(os.path.join(folder, "set.*")))
    all_cells = []
    all_coords = []
    for ii in sets:
        cells, coords = _load_set(ii, data.get("nopbc", False))
        nframes = np.reshape(cells, [-1, 3, 3]).shape[0]
        all_cells.append(np.reshape(cells, [nframes, 3, 3]))
        all_coords.append(np.reshape(coords, [nframes, -1, 3]))
    data["cells"] = np.concatenate(all_cells, axis=0)
    data["coords"] = np.concatenate(all_coords, axis=0)
    # allow custom dtypes
    if labels:
        dtypes = dpdata.system.LabeledSystem.DTYPES
    else:
        dtypes = dpdata.system.System.DTYPES

    for dtype in dtypes:
        if dtype.name in (
            "atom_numbs",
            "atom_names",
            "atom_types",
            "orig",
            "cells",
            "coords",
            "real_atom_names",
            "nopbc",
        ):
            # skip as these data contains specific rules
            continue
        if not (len(dtype.shape) and dtype.shape[0] == dpdata.system.Axis.NFRAMES):
            warnings.warn(
                f"Shape of {dtype.name} is not (nframes, ...), but {dtype.shape}. This type of data will not converted from deepmd/npy format."
            )
            continue
        natoms = data["atom_types"].shape[0]
        shape = [
            natoms if xx == dpdata.system.Axis.NATOMS else xx for xx in dtype.shape[1:]
        ]
        all_data = []
        for ii in sets:
            tmp = _cond_load_data(os.path.join(ii, dtype.deepmd_name + ".npy"))
            if tmp is not None:
                all_data.append(np.reshape(tmp, [tmp.shape[0], *shape]))
        if len(all_data) > 0:
            data[dtype.name] = np.concatenate(all_data, axis=0)
    return data


def dump(folder, data, set_size=5000, comp_prec=np.float32, remove_sets=True):
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
    # dump raw
    np.savetxt(os.path.join(folder, "type.raw"), data["atom_types"], fmt="%d")
    np.savetxt(os.path.join(folder, "type_map.raw"), data["atom_names"], fmt="%s")
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
    # dump frame properties: cell, coord, energy, force and virial
    nsets = nframes // set_size
    if set_size * nsets < nframes:
        nsets += 1
    for ii in range(nsets):
        set_stt = ii * set_size
        set_end = (ii + 1) * set_size
        set_folder = os.path.join(folder, "set.%03d" % ii)  # noqa: UP031
        os.makedirs(set_folder)
    try:
        os.remove(os.path.join(folder, "nopbc"))
    except OSError:
        pass
    if data.get("nopbc", False):
        with open_file(os.path.join(folder, "nopbc"), "w") as fw_nopbc:
            pass
    # allow custom dtypes
    labels = "energies" in data
    if labels:
        dtypes = dpdata.system.LabeledSystem.DTYPES
    else:
        dtypes = dpdata.system.System.DTYPES
    for dtype in dtypes:
        if dtype.name in (
            "atom_numbs",
            "atom_names",
            "atom_types",
            "orig",
            "real_atom_names",
            "nopbc",
        ):
            # skip as these data contains specific rules
            continue
        if dtype.name not in data:
            continue
        if not (len(dtype.shape) and dtype.shape[0] == dpdata.system.Axis.NFRAMES):
            warnings.warn(
                f"Shape of {dtype.name} is not (nframes, ...), but {dtype.shape}. This type of data will not converted to deepmd/npy format."
            )
            continue
        ddata = np.reshape(data[dtype.name], [nframes, -1])
        if np.issubdtype(ddata.dtype, np.floating):
            ddata = ddata.astype(comp_prec)
        for ii in range(nsets):
            set_stt = ii * set_size
            set_end = (ii + 1) * set_size
            set_folder = os.path.join(folder, "set.%03d" % ii)  # noqa: UP031
            np.save(os.path.join(set_folder, dtype.deepmd_name), ddata[set_stt:set_end])
