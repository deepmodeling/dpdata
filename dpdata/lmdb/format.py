from __future__ import annotations

import os

import lmdb
import msgpack
import msgpack_numpy as m
import numpy as np

import dpdata
from dpdata.format import Format

m.patch()


class LMDBError(Exception):
    """Base class for LMDB errors."""


class LMDBMetadataError(LMDBError):
    """Metadata not found in LMDB."""


class LMDBFrameError(LMDBError):
    """Frame data not found in LMDB."""


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

    def to_multi_systems(
        self, systems, file_name, map_size=1000000000, frame_idx_fmt="012d", **kwargs
    ):
        """Save multiple systems to a single LMDB database.

        Parameters
        ----------
        systems : list of dpdata.System
            A list of System objects to be saved.
        file_name : str
            The path to the LMDB database directory. It will be created if it
            doesn't exist.
        map_size : int, optional
            Maximum size of the LMDB database in bytes. Default is 1GB.
        frame_idx_fmt : str, optional
            The format string used to encode the frame index as a key. Default is "012d".
        """
        from dpdata.data_type import Axis

        os.makedirs(file_name, exist_ok=True)
        with lmdb.open(file_name, map_size=map_size) as env:
            global_frame_idx = 0
            system_info = []

            with env.begin(write=True) as txn:
                for system_obj in systems:
                    data = system_obj.data
                    nframes = system_obj.get_nframes()
                    formula = system_obj.formula

                    # Identify symbolic shapes and frame-dependent keys
                    data_shapes = {}
                    frame_dependent_keys = []
                    for dt in type(system_obj).DTYPES:
                        if dt.name in data:
                            if dt.shape is not None:
                                data_shapes[dt.name] = [
                                    s.value if isinstance(s, Axis) else s
                                    for s in dt.shape
                                ]
                                if Axis.NFRAMES in dt.shape:
                                    frame_dependent_keys.append(dt.name)
                            else:
                                data_shapes[dt.name] = None

                    system_info.append(
                        {
                            "formula": formula,
                            "natoms": system_obj.get_atom_numbs(),
                            "nframes": nframes,
                            "start_idx": global_frame_idx,
                            "data_shapes": data_shapes,
                            "frame_dependent_keys": frame_dependent_keys,
                        }
                    )

                    for i in range(nframes):
                        frame_data = {}
                        for key, val in data.items():
                            if key in frame_dependent_keys:
                                frame_data[key] = val[i]
                            else:
                                frame_data[key] = val

                        key = f"{global_frame_idx:{frame_idx_fmt}}".encode("ascii")
                        value = msgpack.packb(frame_data, use_bin_type=True)
                        txn.put(key, value)
                        global_frame_idx += 1

                metadata = {
                    "nframes": global_frame_idx,
                    "system_info": system_info,
                    "frame_idx_fmt": frame_idx_fmt,
                }
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

    def from_multi_systems(self, file_name, map_size=1000000000, **kwargs):
        """Load multiple systems from a single LMDB database.

        Parameters
        ----------
        file_name : str
            The path to the LMDB database directory.
        map_size : int, optional
            Maximum size of the LMDB database in bytes. This parameter is included
            for consistency with `to_multi_systems` but is generally ignored
            when opening in `readonly=True` mode.

        Returns
        -------
        dpdata.MultiSystems
            A MultiSystems object containing all systems stored in the LMDB.
        """
        from dpdata.data_type import Axis, DataType
        from dpdata.system import LabeledSystem, System

        systems = []
        with lmdb.open(file_name, readonly=True) as env:
            with env.begin() as txn:
                metadata_packed = txn.get(b"__metadata__")
                if metadata_packed is None:
                    raise LMDBMetadataError("LMDB database does not contain metadata.")
                metadata = msgpack.unpackb(metadata_packed, raw=False)
                frame_idx_fmt = metadata.get("frame_idx_fmt", "012d")

                for sys_info in metadata["system_info"]:
                    system_frames = []
                    start_idx = sys_info["start_idx"]
                    nframes = sys_info["nframes"]
                    data_shapes = sys_info.get("data_shapes", {})
                    frame_dependent_keys = sys_info.get("frame_dependent_keys", [])

                    for i in range(start_idx, start_idx + nframes):
                        key = f"{i:{frame_idx_fmt}}".encode("ascii")
                        value = txn.get(key)
                        if value is None:
                            raise LMDBFrameError(f"Frame data not found for key: {key}")
                        frame_data = msgpack.unpackb(value, raw=False)
                        system_frames.append(frame_data)

                    # Aggregate data for one system
                    first_frame = system_frames[0]
                    is_labeled = "energies" in first_frame
                    cls = LabeledSystem if is_labeled else System

                    # Auto-register unknown data types
                    existing_dt_names = [dt.name for dt in cls.DTYPES]
                    new_dts = []
                    axis_map = {a.value: a for a in Axis}
                    for key, val in first_frame.items():
                        if key not in existing_dt_names and key in data_shapes:
                            shape_raw = data_shapes[key]
                            if shape_raw is not None:
                                shape = tuple([axis_map.get(s, s) for s in shape_raw])
                            else:
                                shape = None

                            v_arr = np.array(val)
                            new_dts.append(
                                DataType(key, type(v_arr), shape=shape, required=False)
                            )

                    if new_dts:
                        cls.register_data_type(*new_dts)

                    agg_data = {}
                    for key, val in first_frame.items():
                        if key in frame_dependent_keys:
                            agg_data[key] = np.array([d[key] for d in system_frames])
                        else:
                            agg_data[key] = val

                    systems.append(cls(data=agg_data))

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
