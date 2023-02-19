import os
import shutil
import unittest
from itertools import permutations

import numpy as np
from comp_sys import CompLabeledSys, CompSys, IsNoPBC, MultiSystems
from context import dpdata


class TestMixedMultiSystems(unittest.TestCase, CompLabeledSys, IsNoPBC):
    def setUp(self):
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

        system_1 = dpdata.LabeledSystem(
            "gaussian/methane.gaussianlog", fmt="gaussian/log"
        )
        system_2 = dpdata.LabeledSystem(
            "gaussian/methane_reordered.gaussianlog", fmt="gaussian/log"
        )
        system_3 = dpdata.LabeledSystem(
            "gaussian/methane_sub.gaussianlog", fmt="gaussian/log"
        )
        system_4 = dpdata.LabeledSystem(
            "gaussian/noncoveraged.gaussianlog", fmt="gaussian/log"
        )

        self.systems = dpdata.MultiSystems(system_1, system_2, system_3, system_4)
        self.systems.to_deepmd_mixed('tmp.deepmd.mixed')
        mixms = dpdata.MultiSystems().load_systems_from_file('tmp.deepmd.mixed', fmt='deepmd/mixed')
        self.system_1 = self.systems["C1H3"]
        self.system_2 = mixms["C1H3"]
        self.places = 6
        self.e_places = 6
        self.f_places = 6

    def tearDown(self):
        if os.path.exists("tmp.deepmd.npy"):
            shutil.rmtree("tmp.deepmd.npy")


if __name__ == "__main__":
    unittest.main()
