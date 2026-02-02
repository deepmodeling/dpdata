from __future__ import annotations

import os

import lmdb
import msgpack
import msgpack_numpy as m
import numpy as np

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

    Both single systems and multiple systems are supported via the standard
    ``dpdata`` APIs.

    Examples
    --------
    **Saving a single LabeledSystem**

    >>> import dpdata
    >>> system = dpdata.LabeledSystem("path/to/input.vasp", fmt="vasp/outcar")
    >>> system.to("lmdb", "my_single_system.lmdb")

    **Loading a single LabeledSystem**

    >>> loaded_system = dpdata.LabeledSystem("my_single_system.lmdb", fmt="lmdb")

    **Saving multiple systems to a single LMDB database**

    >>> import dpdata
    >>> system_1 = dpdata.LabeledSystem("path/to/system1/OUTCAR", fmt="vasp/outcar")
    >>> system_2 = dpdata.LabeledSystem("path/to/system2/OUTCAR", fmt="vasp/outcar")
    >>> multi_systems_obj = dpdata.MultiSystems(system_1, system_2)
    >>> multi_systems_obj.to("lmdb", "my_multi_system_db.lmdb")

    **Loading multiple systems from a single LMDB database**

    >>> import dpdata
    >>> loaded_multi_systems = dpdata.MultiSystems.from_file("my_multi_system_db.lmdb", fmt="lmdb")
    """

    def to_multi_systems(
        self, formulas, directory, map_size=1000000000, frame_idx_fmt="012d", **kwargs
    ):
        """Implement MultiSystems.to for LMDB format.

        Parameters
        ----------
        formulas : list[str]
            list of formulas
        directory : str
            directory of system
        map_size : int, optional
            Maximum size of the LMDB database in bytes. Default is 1GB.
        frame_idx_fmt : str, optional
            The format string used to encode the frame index as a key. Default is "012d".
        **kwargs : dict
            other parameters

        Yields
        ------
        tuple
            (self, formula) to be used by to_system
        """
        self._frame_idx_fmt = frame_idx_fmt
        self._global_frame_idx = 0
        self._system_info = []
        os.makedirs(directory, exist_ok=True)
        with lmdb.open(directory, map_size=map_size) as env:
            with env.begin(write=True) as txn:
                self._txn = txn
                for ff in formulas:
                    yield (self, ff)
                # Finalize metadata
                metadata = {
                    "nframes": self._global_frame_idx,
                    "system_info": self._system_info,
                    "frame_idx_fmt": self._frame_idx_fmt,
                }
                txn.put(b"__metadata__", msgpack.packb(metadata, use_bin_type=True))
                self._txn = None

    def _dump_to_txn(self, data, txn, formula, dtypes):
        from dpdata.data_type import Axis

        nframes = data["coords"].shape[0]

        # Identify symbolic shapes and frame-dependent keys
        data_shapes = {}
        frame_dependent_keys = []
        for dt in dtypes:
            if dt.name in data:
                if dt.shape is not None:
                    data_shapes[dt.name] = [
                        s.value if isinstance(s, Axis) else s for s in dt.shape
                    ]
                    if Axis.NFRAMES in dt.shape:
                        frame_dependent_keys.append(dt.name)
                else:
                    data_shapes[dt.name] = None

        # Record system info
        # natoms needs to be extracted from data
        if "atom_numbs" in data:
            natoms_list = data["atom_numbs"]
        else:
            # Fallback for systems without atom_numbs (should not happen in valid dpdata systems)
            natoms_list = []

        self._system_info.append(
            {
                "formula": formula,
                "natoms": natoms_list,
                "nframes": nframes,
                "start_idx": self._global_frame_idx,
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

            key = f"{self._global_frame_idx:{self._frame_idx_fmt}}".encode("ascii")
            value = msgpack.packb(frame_data, use_bin_type=True)
            txn.put(key, value)
            self._global_frame_idx += 1

    def to_labeled_system(self, data, file_name, **kwargs):
        """Save a single LabeledSystem to an LMDB database."""
        from dpdata.system import LabeledSystem

        if isinstance(file_name, tuple) and file_name[0] is self:
            txn, formula = self._txn, file_name[1]
            self._dump_to_txn(data, txn, formula, LabeledSystem.DTYPES)
        else:
            # Single system call: use to_multi_systems logic
            # Infer formula from data if possible, or use default
            formula = kwargs.get("formula", "unknown")
            gen = self.to_multi_systems([formula], file_name, **kwargs)
            handle = next(gen)
            self.to_labeled_system(data, handle, **kwargs)
            try:
                next(gen)
            except StopIteration:
                pass

    def to_system(self, data, file_name, **kwargs):
        """Save a single System to an LMDB database."""
        from dpdata.system import System

        if isinstance(file_name, tuple) and file_name[0] is self:
            txn, formula = self._txn, file_name[1]
            self._dump_to_txn(data, txn, formula, System.DTYPES)
        else:
            # Single system call
            formula = kwargs.get("formula", "unknown")
            gen = self.to_multi_systems([formula], file_name, **kwargs)
            handle = next(gen)
            self.to_system(data, handle, **kwargs)
            try:
                next(gen)
            except StopIteration:
                pass

    def from_multi_systems(self, file_name, map_size=1000000000, **kwargs):
        """Load multiple systems from a single LMDB database.

        Parameters
        ----------
        file_name : str
            The path to the LMDB database directory.
        map_size : int, optional
            Maximum size of the LMDB database in bytes.
        **kwargs : dict
            other parameters

        Yields
        ------
        dict
            data dictionary for each system
        """
        from dpdata.data_type import Axis, DataType
        from dpdata.system import LabeledSystem, System

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

                    yield agg_data

    def from_labeled_system(self, file_name, **kwargs):
        """Load data for a single LabeledSystem from an LMDB database."""
        if isinstance(file_name, dict):
            return file_name
        # from_multi_systems returns a generator of dicts
        gen = self.from_multi_systems(file_name, **kwargs)
        return next(gen)

    def from_system(self, file_name, **kwargs):
        """Load data for a single System from an LMDB database."""
        if isinstance(file_name, dict):
            return file_name
        # from_multi_systems returns a generator of dicts
        gen = self.from_multi_systems(file_name, **kwargs)
        return next(gen)
