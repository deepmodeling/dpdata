"""Utils for deepmd/hdf5 format."""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import numpy as np

import dpdata

if TYPE_CHECKING:
    import h5py

__all__ = ["to_system_data", "dump"]


def to_system_data(
    f: h5py.File | h5py.Group,
    folder: str,
    type_map: list | None = None,
    labels: bool = True,
):
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
    from wcmatch.glob import globfilter

    g = f[folder] if folder else f

    data = {}
    # ignore empty files or groups
    if "type.raw" not in g.keys():
        return data
    data["atom_types"] = g["type.raw"][:]
    ntypes = np.max(data["atom_types"]) + 1
    natoms = data["atom_types"].size
    data["atom_numbs"] = []
    for ii in range(ntypes):
        data["atom_numbs"].append(np.count_nonzero(data["atom_types"] == ii))
    data["atom_names"] = []
    # if find type_map.raw, use it
    if "type_map.raw" in g.keys():
        my_type_map = list(np.char.decode(g["type_map.raw"][:]))
    # else try to use arg type_map
    elif type_map is not None:
        my_type_map = type_map
    # in the last case, make artificial atom names
    else:
        my_type_map = []
        for ii in range(ntypes):
            my_type_map.append("Type_%d" % ii)  # noqa: UP031
    assert len(my_type_map) >= len(data["atom_numbs"])
    for ii in range(len(data["atom_numbs"])):
        data["atom_names"].append(my_type_map[ii])

    data["orig"] = np.zeros([3])
    if "nopbc" in g.keys():
        data["nopbc"] = True
    sets = globfilter(g.keys(), "set.*")

    data_types = {}
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
            "real_atom_types",
            "real_atom_names",
            "nopbc",
        ):
            # skip as these data contains specific rules
            continue
        if not (len(dtype.shape) and dtype.shape[0] == dpdata.system.Axis.NFRAMES):
            warnings.warn(
                f"Shape of {dtype.name} is not (nframes, ...), but {dtype.shape}. This type of data will not converted from deepmd/hdf5 format."
            )
            continue
        shape = [
            natoms if xx == dpdata.system.Axis.NATOMS else xx for xx in dtype.shape[1:]
        ]

        data_types[dtype.name] = {
            "fn": dtype.deepmd_name,
            "shape": shape,
            "required": dtype.required
            and not (dtype.name == "cells" and data.get("nopbc", False)),
        }

    for dt, prop in data_types.items():
        all_data = []

        for ii in sets:
            set = g[ii]
            fn = "{}.npy".format(prop["fn"])
            if fn in set.keys():
                dd = set[fn][:]
                nframes = dd.shape[0]
                all_data.append(np.reshape(dd, (nframes, *prop["shape"])))
            elif prop["required"]:
                raise RuntimeError(f"{folder}/{ii}/{fn} not found")

        if len(all_data) > 0:
            data[dt] = np.concatenate(all_data, axis=0)
    if "cells" not in data:
        nframes = data["coords"].shape[0]
        data["cells"] = np.zeros((nframes, 3, 3))
    return data


def dump(
    f: h5py.File | h5py.Group,
    folder: str,
    data: dict,
    set_size=5000,
    comp_prec=np.float32,
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
    # ignore empty systems
    if not len(data["coords"]):
        return
    # dump raw (array in fact)
    g.create_dataset("type.raw", data=data["atom_types"])
    g.create_dataset("type_map.raw", data=np.array(data["atom_names"], dtype="S"))
    # BondOrder System
    if "bonds" in data:
        g.create_dataset("bonds.raw", data=data["bonds"])
    if "formal_charges" in data:
        g.create_dataset("formal_charges.raw", data=data["formal_charges"])
    # reshape frame properties and convert prec
    nframes = data["cells"].shape[0]

    nopbc = data.get("nopbc", False)
    reshaped_data = {}

    data_types = {}

    labels = "energies" in data
    if labels:
        dtypes = dpdata.system.LabeledSystem.DTYPES
    else:
        dtypes = dpdata.system.System.DTYPES
    # allow custom dtypes
    for dtype in dtypes:
        if dtype.name in (
            "atom_numbs",
            "atom_names",
            "atom_types",
            "orig",
            "real_atom_types",
            "real_atom_names",
            "nopbc",
        ):
            # skip as these data contains specific rules
            continue
        if not (len(dtype.shape) and dtype.shape[0] == dpdata.system.Axis.NFRAMES):
            warnings.warn(
                f"Shape of {dtype.name} is not (nframes, ...), but {dtype.shape}. This type of data will not converted to deepmd/hdf5 format."
            )
            continue

        data_types[dtype.name] = {
            "fn": dtype.deepmd_name,
            "shape": (nframes, -1),
            "dump": not (dtype.name == "cells" and nopbc),
        }

    for dt, prop in data_types.items():
        if dt in data:
            if prop["dump"]:
                ddata = np.reshape(data[dt], prop["shape"])
                if np.issubdtype(ddata.dtype, np.floating):
                    ddata = ddata.astype(comp_prec)
                reshaped_data[dt] = ddata

    # dump frame properties: cell, coord, energy, force and virial
    nsets = nframes // set_size
    if set_size * nsets < nframes:
        nsets += 1
    for ii in range(nsets):
        set_stt = ii * set_size
        set_end = (ii + 1) * set_size
        set_folder = g.create_group("set.%03d" % ii)  # noqa: UP031
        for dt, prop in data_types.items():
            if dt in reshaped_data:
                set_folder.create_dataset(
                    "{}.npy".format(prop["fn"]), data=reshaped_data[dt][set_stt:set_end]
                )

    if nopbc:
        g.create_dataset("nopbc", data=True)
