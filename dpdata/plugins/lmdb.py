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

    Attributes
    ----------
    See `dpdata.format.Format` for conventions.

    How to Use
    ----------

    **1. Single System (`System` or `LabeledSystem`)**

    For single systems, the standard `dpdata` API can be used.

    - **Saving to LMDB:**
      ```python
      import dpdata

      # Load your system from any supported format
      system = dpdata.LabeledSystem("path/to/your/input.vasp", fmt="vasp/outcar")

      # Save to an LMDB database. The path will be treated as a directory.
      system.to("lmdb", "my_single_system.lmdb")
      ```

    - **Loading from LMDB:**
      ```python
      import dpdata

      # Load the system from the LMDB directory
      loaded_system = dpdata.LabeledSystem("my_single_system.lmdb", fmt="lmdb")
      ```

    **2. Multiple Systems (`MultiSystems`) in a Single Database**

    To store multiple systems in a single LMDB file for efficient sampling,
    you must call the format plugin's methods directly. This is because the
    standard `MultiSystems.to()` and `MultiSystems(fmt=...)` APIs are
    architecturally designed for a "directory of files" paradigm, where each
    system is written to a separate file. This conflicts with our goal of
    creating a single, unified database.

    - **Saving `MultiSystems` to a single LMDB file:**
      ```python
      import dpdata
      from dpdata.plugins.lmdb import LMDBFormat

      # Create your individual systems
      system_1 = dpdata.LabeledSystem("path/to/system1.log", fmt="gaussian/log")
      system_2 = dpdata.LabeledSystem("path/to/system2.out", fmt="vasp/outcar")

      multi_systems_obj = dpdata.MultiSystems(system_1, system_2)

      # Instantiate the LMDBFormat plugin and call its method directly
      lmdb_formatter = LMDBFormat()
      lmdb_formatter.to_multi_systems(
          list(multi_systems_obj.systems.values()), "my_multi_system_db.lmdb"
      )
      ```

    - **Loading `MultiSystems` from a single LMDB file:**
      ```python
      import dpdata
      from dpdata.plugins.lmdb import LMDBFormat

      # Instantiate the LMDBFormat plugin and call its method directly
      lmdb_formatter = LMDBFormat()
      loaded_multi_systems = lmdb_formatter.from_multi_systems("my_multi_system_db.lmdb")
      ```
    """

    def to_multi_systems(self, systems, file_name, **kwargs):
        """Save multiple systems to a single LMDB database."""
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
        from dpdata.system import LabeledSystem

        self.to_multi_systems([LabeledSystem(data=data)], file_name, **kwargs)

    def to_system(self, data, file_name, **kwargs):
        from dpdata.system import System

        self.to_multi_systems([System(data=data)], file_name, **kwargs)

    def from_multi_systems(self, file_name, **kwargs):
        """Load multiple systems from a single LMDB database."""
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
        # from_multi_systems returns a MultiSystems object
        multisystems_obj = self.from_multi_systems(file_name, **kwargs)
        # We need the data dictionary of the first (and only) system for from_labeled_system
        return multisystems_obj[0].data

    def from_system(self, file_name, **kwargs):
        # from_multi_systems returns a MultiSystems object
        multisystems_obj = self.from_multi_systems(file_name, **kwargs)
        # We need the data dictionary of the first (and only) system for from_system
        return multisystems_obj[0].data
