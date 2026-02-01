from __future__ import annotations

import os

import lmdb
import msgpack
import msgpack_numpy as m
import numpy as np

import dpdata
from dpdata.format import Format

m.patch()


@Format.register("lmdb")
class LMDBFormat(Format):
    """
    Class for handling the LMDB format, which stores atomic configurations in a
    Lightning Memory-Mapped Database (LMDB).

    This format is optimized for machine learning workflows where fast, random
    access to a large number of frames is required. All frames from multiple
    systems (with potentially different numbers of atoms) are stored in a
    single LMDB database file.

    For single systems, the standard ``dpdata.System.to('lmdb', ...)`` and
    ``dpdata.System('...', fmt='lmdb')`` APIs can be used.

    .. note::

        The standard ``dpdata.MultiSystems.to()`` and ``dpdata.MultiSystems()``
        constructor are not supported for this format. This is due to an
        architectural limitation, as those APIs are designed for a "directory
        of files" paradigm, whereas this format creates a single, unified
        database file. To save and load multiple systems, you must call the
        format's methods directly, as shown in the examples below.

    Examples
    --------
    **Saving a single LabeledSystem**

    >>> import dpdata
    >>> system = dpdata.LabeledSystem("path/to/input.vasp", fmt="vasp/outcar")
    >>> system.to("lmdb", "my_single_system.lmdb")

    **Loading a single LabeledSystem**

    >>> loaded_system = dpdata.LabeledSystem("my_single_system.lmdb", fmt="lmdb")

    **Saving multiple systems to a single LMDB database**

    >>> from dpdata.plugins.lmdb import LMDBFormat
    >>> system_1 = dpdata.LabeledSystem("path/to/system1/OUTCAR", fmt="vasp/outcar")
    >>> system_2 = dpdata.LabeledSystem("path/to/system2/OUTCAR", fmt="vasp/outcar")
    >>> multi_systems_obj = dpdata.MultiSystems(system_1, system_2)
    >>> lmdb_formatter = LMDBFormat()
    >>> lmdb_formatter.to_multi_systems(
    ...     list(multi_systems_obj.systems.values()), "my_multi_system_db.lmdb"
    ... )

    **Loading multiple systems from a single LMDB database**

    >>> from dpdata.plugins.lmdb import LMDBFormat
    >>> lmdb_formatter = LMDBFormat()
    >>> loaded_multi_systems = lmdb_formatter.from_multi_systems("my_multi_system_db.lmdb")
    """

    def to_multi_systems(self, systems, file_name, **kwargs):
        """Save multiple systems to a single LMDB database.

        Parameters
        ----------
        systems : list of dpdata.System
            A list of System objects to be saved.
        file_name : str
            The path to the LMDB database directory. It will be created if it
            doesn't exist.
        """
        os.makedirs(file_name, exist_ok=True)
        env = lmdb.open(file_name, map_size=1000000000)
        global_frame_idx = 0
        system_info = []

        with env.begin(write=True) as txn:
            for system_obj in systems:
                data = system_obj.data
                nframes = data["coords"].shape[0]
                natoms = data["atom_numbs"]
                formula = system_obj.formula

                system_info.append(
                    {
                        "formula": formula,
                        "natoms": natoms,
                        "nframes": nframes,
                        "start_idx": global_frame_idx,
                    }
                )

                for i in range(nframes):
                    frame_data = {
                        "atom_names": data["atom_names"],
                        "atom_numbs": data["atom_numbs"],
                        "atom_types": data["atom_types"],
                        "orig": data["orig"],
                        "coords": data["coords"][i],
                        "cells": data["cells"][i],
                    }
                    if "energies" in data:
                        frame_data["energies"] = data["energies"][i]
                    if "forces" in data:
                        frame_data["forces"] = data["forces"][i]
                    if "virials" in data:
                        frame_data["virials"] = data["virials"][i]
                    if "nopbc" in data:
                        frame_data["nopbc"] = data["nopbc"]

                    key = f"{global_frame_idx:012d}".encode("ascii")
                    value = msgpack.packb(frame_data, use_bin_type=True)
                    txn.put(key, value)
                    global_frame_idx += 1

            metadata = {"nframes": global_frame_idx, "system_info": system_info}
            txn.put(b"__metadata__", msgpack.packb(metadata, use_bin_type=True))

    def to_labeled_system(self, data, file_name, **kwargs):
        """Save a single LabeledSystem to an LMDB database.

        Parameters
        ----------
        data : dict
            The data dictionary of a LabeledSystem.
        file_name : str
            The path to the LMDB database directory.
        """
        from dpdata.system import LabeledSystem

        self.to_multi_systems([LabeledSystem(data=data)], file_name, **kwargs)

    def to_system(self, data, file_name, **kwargs):
        """Save a single System to an LMDB database.

        Parameters
        ----------
        data : dict
            The data dictionary of a System.
        file_name : str
            The path to the LMDB database directory.
        """
        from dpdata.system import System

        self.to_multi_systems([System(data=data)], file_name, **kwargs)

    def from_multi_systems(self, file_name, **kwargs):
        """Load multiple systems from a single LMDB database.

        Parameters
        ----------
        file_name : str
            The path to the LMDB database directory.

        Returns
        -------
        dpdata.MultiSystems
            A MultiSystems object containing all systems stored in the LMDB.
        """
        from dpdata.system import LabeledSystem, System

        systems = []
        env = lmdb.open(file_name, readonly=True)

        with env.begin() as txn:
            metadata_packed = txn.get(b"__metadata__")
            metadata = msgpack.unpackb(metadata_packed, raw=False)

            for sys_info in metadata["system_info"]:
                system_frames = []
                start_idx = sys_info["start_idx"]
                nframes = sys_info["nframes"]

                for i in range(start_idx, start_idx + nframes):
                    key = f"{i:012d}".encode("ascii")
                    value = txn.get(key)
                    frame_data = msgpack.unpackb(value, raw=False)
                    system_frames.append(frame_data)

                # Aggregate data for one system
                agg_data = {
                    "atom_names": system_frames[0]["atom_names"],
                    "atom_numbs": system_frames[0]["atom_numbs"],
                    "atom_types": system_frames[0]["atom_types"],
                    "orig": system_frames[0]["orig"],
                    "coords": np.array([d["coords"] for d in system_frames]),
                    "cells": np.array([d["cells"] for d in system_frames]),
                }
                is_labeled = "energies" in system_frames[0]
                if is_labeled:
                    agg_data["energies"] = np.array(
                        [d["energies"] for d in system_frames]
                    )
                    agg_data["forces"] = np.array([d["forces"] for d in system_frames])
                    if "virials" in system_frames[0]:
                        agg_data["virials"] = np.array(
                            [d.get("virials") for d in system_frames]
                        )
                    if "nopbc" in system_frames[0]:
                        agg_data["nopbc"] = system_frames[0]["nopbc"]
                    systems.append(LabeledSystem(data=agg_data))
                else:
                    if "nopbc" in system_frames[0]:
                        agg_data["nopbc"] = system_frames[0]["nopbc"]
                    systems.append(System(data=agg_data))

        return dpdata.MultiSystems(*systems)

    def from_labeled_system(self, file_name, **kwargs):
        """Load data for a single LabeledSystem from an LMDB database.

        Parameters
        ----------
        file_name : str
            The path to the LMDB database directory.

        Returns
        -------
        dict
            The data dictionary for the loaded LabeledSystem.
        """
        # from_multi_systems returns a MultiSystems object
        multisystems_obj = self.from_multi_systems(file_name, **kwargs)
        # We need the data dictionary of the first (and only) system for from_labeled_system
        return multisystems_obj[0].data

    def from_system(self, file_name, **kwargs):
        """Load data for a single System from an LMDB database.

        Parameters
        ----------
        file_name : str
            The path to the LMDB database directory.

        Returns
        -------
        dict
            The data dictionary for the loaded System.
        """
        # from_multi_systems returns a MultiSystems object
        multisystems_obj = self.from_multi_systems(file_name, **kwargs)
        # We need the data dictionary of the first (and only) system for from_system
        return multisystems_obj[0].data
