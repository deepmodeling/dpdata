from __future__ import annotations

import os
import shutil
import unittest

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


if __name__ == "__main__":
    unittest.main()
