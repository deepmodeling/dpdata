from __future__ import annotations

import os
import shutil
import unittest

import lmdb
import msgpack
import msgpack_numpy as m
from comp_sys import (
    CompLabeledMultiSys,
    CompLabeledSys,
    CompSys,
    IsPBC,
    MSAllIsNoPBC,
)
from context import dpdata

from dpdata.plugins.lmdb import LMDBFormat


class TestLMDBLabeledSystem(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        self.system_1 = dpdata.LabeledSystem("poscars/OUTCAR.h2o.md", fmt="vasp/outcar")
        self.lmdb_path = "tmp_labeled.lmdb"
        self.system_1.to("lmdb", self.lmdb_path)
        self.system_2 = dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

    def tearDown(self):
        if os.path.exists(self.lmdb_path):
            shutil.rmtree(self.lmdb_path)


class TestLMDBSystem(unittest.TestCase, CompSys, IsPBC):
    def setUp(self):
        self.system_1 = dpdata.System("poscars/POSCAR.h2o.md", fmt="vasp/poscar")
        self.lmdb_path = "tmp_system.lmdb"
        self.system_1.to("lmdb", self.lmdb_path)
        self.system_2 = dpdata.System(self.lmdb_path, fmt="lmdb")
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

    def tearDown(self):
        if os.path.exists(self.lmdb_path):
            shutil.rmtree(self.lmdb_path)


class TestLMDBMultiSystems(unittest.TestCase, CompLabeledMultiSys, MSAllIsNoPBC):
    def setUp(self):
        self.lmdb_path = "tmp_multi.lmdb"
        system_1 = dpdata.LabeledSystem(
            "gaussian/methane.gaussianlog", fmt="gaussian/log"
        )
        system_2 = dpdata.LabeledSystem(
            "gaussian/methane_reordered.gaussianlog", fmt="gaussian/log"
        )
        system_3 = dpdata.LabeledSystem(
            "gaussian/methane_sub.gaussianlog", fmt="gaussian/log"
        )

        self.ms_1 = dpdata.MultiSystems(system_1, system_2, system_3)

        # Manually call the format object's to_multi_systems as dpdata.MultiSystems.to
        # is not designed for single-file multi-system formats.
        LMDBFormat().to_multi_systems(list(self.ms_1.systems.values()), self.lmdb_path)

        # Manually call the format object's from_multi_systems as dpdata.MultiSystems.from_fmt
        # is not designed for single-file multi-system formats.
        self.ms_2 = LMDBFormat().from_multi_systems(self.lmdb_path)

        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

    def tearDown(self):
        if os.path.exists(self.lmdb_path):
            shutil.rmtree(self.lmdb_path)


class TestLMDBErrorHandling(unittest.TestCase):
    def setUp(self):
        self.lmdb_path_missing_metadata = "tmp_missing_metadata.lmdb"
        self.lmdb_path_missing_frame = "tmp_missing_frame.lmdb"

        # Ensure cleanup in case of previous test failures
        if os.path.exists(self.lmdb_path_missing_metadata):
            shutil.rmtree(self.lmdb_path_missing_metadata)
        if os.path.exists(self.lmdb_path_missing_frame):
            shutil.rmtree(self.lmdb_path_missing_frame)

        # For test_load_missing_frame_data, create a valid LMDB environment
        # and write metadata, but no actual frames.
        env = lmdb.open(self.lmdb_path_missing_frame, map_size=1000000000)
        with env.begin(write=True) as txn:
            metadata = {
                "nframes": 1,
                "system_info": [
                    {
                        "formula": "H2O",
                        "natoms": [1, 2],
                        "nframes": 1,
                        "start_idx": 0,
                    }
                ],
            }
            m.patch()  # Ensure numpy patching for metadata
            txn.put(b"__metadata__", msgpack.packb(metadata, use_bin_type=True))
        env.close()

    def tearDown(self):
        if os.path.exists(self.lmdb_path_missing_metadata):
            shutil.rmtree(self.lmdb_path_missing_metadata)
        if os.path.exists(self.lmdb_path_missing_frame):
            shutil.rmtree(self.lmdb_path_missing_frame)

    def test_load_missing_metadata(self):
        # Create a valid, empty LMDB environment, then test for missing metadata
        lmdb.open(
            self.lmdb_path_missing_metadata, map_size=1000000000
        ).close()  # Creates empty LMDB environment

        with self.assertRaisesRegex(
            KeyError, "LMDB database does not contain metadata."
        ):
            LMDBFormat().from_multi_systems(self.lmdb_path_missing_metadata)

    def test_load_missing_frame_data(self):
        with self.assertRaisesRegex(
            KeyError, "Frame data not found for key: b'000000000000'"
        ):
            LMDBFormat().from_multi_systems(self.lmdb_path_missing_frame)


if __name__ == "__main__":
    unittest.main()
