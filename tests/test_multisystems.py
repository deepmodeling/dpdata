import os
import tempfile
import unittest
from itertools import permutations

import numpy as np
from comp_sys import CompLabeledSys, IsNoPBC, MultiSystems
from context import dpdata


class TestMultiSystems(unittest.TestCase, CompLabeledSys, MultiSystems, IsNoPBC):
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

        self.systems = dpdata.MultiSystems(system_1, system_3, system_4)
        self.systems.append(system_2)
        self.system_1 = self.systems["C1H3"]
        self.system_2 = system_3

        self.system_names = ["C1H4", "C1H3"]
        self.system_sizes = {"C1H4": 2, "C1H3": 1}
        self.atom_names = ["C", "H"]

    def test_len(self):
        self.assertEqual(len(self.systems), 2)

    def test_get_nframes(self):
        self.assertEqual(self.systems.get_nframes(), 3)

    def test_str(self):
        self.assertEqual(
            str(self.systems), "MultiSystems (2 systems containing 3 frames)"
        )


class TestMultiSystemsAdd(unittest.TestCase, CompLabeledSys, MultiSystems, IsNoPBC):
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

        self.systems = dpdata.MultiSystems(system_1)
        self.systems += system_2
        self.systems += system_3
        self.systems += system_4
        for s in self.systems:
            if s.formula == "C1H3":
                self.system_1 = s
        self.system_2 = system_3

        self.system_names = ["C1H4", "C1H3"]
        self.system_sizes = {"C1H4": 2, "C1H3": 1}
        self.atom_names = ["C", "H"]


class TestMultiSystemsSorted(unittest.TestCase, MultiSystems):
    def setUp(self):
        # CH4 and O2
        system_1 = dpdata.LabeledSystem(
            "gaussian/methane.gaussianlog", fmt="gaussian/log"
        )
        system_2 = dpdata.LabeledSystem(
            "gaussian/oxygen.gaussianlog", fmt="gaussian/log"
        )
        self.systems = dpdata.MultiSystems(system_1, system_2)

        self.system_names = ["C1H4O0", "C0H0O2"]
        self.system_sizes = {"C1H4O0": 1, "C0H0O2": 1}
        self.atom_names = ["C", "H", "O"]


class TestMultiDeepmdDumpRaw(unittest.TestCase, CompLabeledSys, IsNoPBC):
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

        systems = dpdata.MultiSystems(system_1, system_2, system_3, system_4)
        path = "tmp.deepmd.multi"
        systems.to_deepmd_raw(path)
        self.system_1 = dpdata.LabeledSystem(
            os.path.join(path, "C1H3"), fmt="deepmd/raw", type_map=["C", "H"]
        )
        self.system_2 = system_3


class TestMultiDeepmdDumpComp(unittest.TestCase, CompLabeledSys, IsNoPBC):
    def setUp(self):
        self.places = 6
        self.e_places = 4
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

        systems = dpdata.MultiSystems(system_1, system_2, system_3, system_4)
        path = "tmp.deepmd.npy.multi"
        systems.to_deepmd_npy(path)
        self.system_1 = dpdata.LabeledSystem(
            os.path.join(path, "C1H3"), fmt="deepmd/npy", type_map=["C", "H"]
        )
        self.system_2 = system_3


class TestTypeMap(unittest.TestCase):
    def setUp(self):
        self.system_1 = dpdata.LabeledSystem(
            "gaussian/methane.gaussianlog", fmt="gaussian/log"
        )
        self.system_2 = dpdata.LabeledSystem(
            "gaussian/methane_reordered.gaussianlog", fmt="gaussian/log"
        )
        self.system_3 = dpdata.LabeledSystem(
            "gaussian/methane_sub.gaussianlog", fmt="gaussian/log"
        )
        self.system_4 = dpdata.LabeledSystem(
            "gaussian/noncoveraged.gaussianlog", fmt="gaussian/log"
        )

    def test_type_map(self):
        for type_map in permutations(["C", "H", "O", "N"], 4):
            systems = dpdata.MultiSystems(
                self.system_1,
                self.system_2,
                self.system_3,
                self.system_4,
                type_map=type_map,
            )
            self.assertEqual(type_map, systems.atom_names)


class TestMultiSystemsTo(unittest.TestCase, MultiSystems):
    def setUp(self):
        # CH4 and O2
        system_1 = dpdata.LabeledSystem(
            "gaussian/methane.gaussianlog", fmt="gaussian/log"
        )
        system_2 = dpdata.LabeledSystem(
            "gaussian/oxygen.gaussianlog", fmt="gaussian/log"
        )
        systems1 = dpdata.MultiSystems(system_1, system_2)
        systems1.to_deepmd_npy("tmp.multi.to")
        self.systems = dpdata.MultiSystems().from_deepmd_npy("tmp.multi.to")

        self.system_names = ["C1H4O0", "C0H0O2"]
        self.system_sizes = {"C1H4O0": 1, "C0H0O2": 1}
        self.atom_names = ["C", "H", "O"]


class TestLongFilename(unittest.TestCase):
    def test_long_filename1(self):
        system = dpdata.System(
            data={
                "atom_names": [f"TYPE{ii}" for ii in range(200)],
                "atom_numbs": [1] + [0 for _ in range(199)],
                "atom_types": np.arange(1),
                "coords": np.zeros((1, 1, 3)),
                "orig": np.zeros(3),
                "cells": np.zeros((1, 3, 3)),
            }
        )
        ms = dpdata.MultiSystems(system)
        with tempfile.TemporaryDirectory() as tmpdir:
            ms.to_deepmd_npy(tmpdir)

    def test_long_filename2(self):
        system = dpdata.System(
            data={
                "atom_names": [f"TYPE{ii}" for ii in range(200)],
                "atom_numbs": [1 for _ in range(200)],
                "atom_types": np.arange(200),
                "coords": np.zeros((1, 200, 3)),
                "orig": np.zeros(3),
                "cells": np.zeros((1, 3, 3)),
            }
        )
        ms = dpdata.MultiSystems(system)
        with tempfile.TemporaryDirectory() as tmpdir:
            ms.to_deepmd_npy(tmpdir)


if __name__ == "__main__":
    unittest.main()
