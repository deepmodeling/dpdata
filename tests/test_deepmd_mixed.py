from __future__ import annotations

import os
import shutil
import unittest
from glob import glob

import numpy as np
from comp_sys import (
    CompLabeledMultiSys,
    CompLabeledSys,
    IsNoPBC,
    MSAllIsNoPBC,
    MultiSystems,
)
from context import dpdata


class TestMixedMultiSystemsDumpLoad(
    unittest.TestCase, CompLabeledMultiSys, MultiSystems, MSAllIsNoPBC
):
    def setUp(self):
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

        # C1H4
        system_1 = dpdata.LabeledSystem(
            "gaussian/methane.gaussianlog", fmt="gaussian/log"
        )

        # C1H3
        system_2 = dpdata.LabeledSystem(
            "gaussian/methane_sub.gaussianlog", fmt="gaussian/log"
        )

        tmp_data = system_1.data.copy()
        tmp_data["atom_numbs"] = [1, 1, 1, 2]
        tmp_data["atom_names"] = ["C", "H", "A", "B"]
        tmp_data["atom_types"] = np.array([0, 1, 2, 3, 3])
        # C1H1A1B2
        system_1_modified_type_1 = dpdata.LabeledSystem(data=tmp_data)

        tmp_data = system_1.data.copy()
        tmp_data["atom_numbs"] = [1, 1, 2, 1]
        tmp_data["atom_names"] = ["C", "H", "A", "B"]
        tmp_data["atom_types"] = np.array([0, 1, 2, 2, 3])
        # C1H1A2B1
        system_1_modified_type_2 = dpdata.LabeledSystem(data=tmp_data)

        tmp_data = system_1.data.copy()
        tmp_data["atom_numbs"] = [1, 1, 1, 2]
        tmp_data["atom_names"] = ["C", "H", "A", "D"]
        tmp_data["atom_types"] = np.array([0, 1, 2, 3, 3])
        # C1H1A1C2
        system_1_modified_type_3 = dpdata.LabeledSystem(data=tmp_data)

        self.ms = dpdata.MultiSystems(
            system_1,
            system_2,
            system_1_modified_type_1,
            system_1_modified_type_2,
            system_1_modified_type_3,
        )
        self.ms.to_deepmd_npy_mixed("tmp.deepmd.mixed")
        self.place_holder_ms = dpdata.MultiSystems()
        self.place_holder_ms.from_deepmd_npy("tmp.deepmd.mixed", fmt="deepmd/npy")
        self.systems = dpdata.MultiSystems()
        self.systems.from_deepmd_npy_mixed("tmp.deepmd.mixed", fmt="deepmd/npy/mixed")
        self.ms_1 = self.ms
        self.ms_2 = self.systems
        mixed_sets = glob("tmp.deepmd.mixed/*/set.*")
        self.assertEqual(len(mixed_sets), 2)
        for i in mixed_sets:
            self.assertEqual(
                os.path.exists(os.path.join(i, "real_atom_types.npy")), True
            )

        self.system_names = [
            "C1H4A0B0D0",
            "C1H3A0B0D0",
            "C1H1A1B2D0",
            "C1H1A2B1D0",
            "C1H1A1B0D2",
        ]
        self.system_sizes = {
            "C1H4A0B0D0": 1,
            "C1H3A0B0D0": 1,
            "C1H1A1B2D0": 1,
            "C1H1A2B1D0": 1,
            "C1H1A1B0D2": 1,
        }
        self.atom_names = ["C", "H", "A", "B", "D"]

    def tearDown(self):
        if os.path.exists("tmp.deepmd.mixed"):
            shutil.rmtree("tmp.deepmd.mixed")

    def test_len(self):
        self.assertEqual(len(self.ms), 5)
        self.assertEqual(len(self.place_holder_ms), 2)
        self.assertEqual(len(self.systems), 5)

    def test_get_nframes(self):
        self.assertEqual(self.ms.get_nframes(), 5)
        self.assertEqual(self.place_holder_ms.get_nframes(), 5)
        self.assertEqual(self.systems.get_nframes(), 5)

    def test_str(self):
        self.assertEqual(str(self.ms), "MultiSystems (5 systems containing 5 frames)")
        self.assertEqual(
            str(self.place_holder_ms), "MultiSystems (2 systems containing 5 frames)"
        )
        self.assertEqual(
            str(self.systems), "MultiSystems (5 systems containing 5 frames)"
        )


class TestMixedMultiSystemsDumpLoadTypeMap(
    unittest.TestCase, CompLabeledMultiSys, MultiSystems, MSAllIsNoPBC
):
    def setUp(self):
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

        # C1H4
        system_1 = dpdata.LabeledSystem(
            "gaussian/methane.gaussianlog", fmt="gaussian/log"
        )

        # C1H3
        system_2 = dpdata.LabeledSystem(
            "gaussian/methane_sub.gaussianlog", fmt="gaussian/log"
        )

        tmp_data = system_1.data.copy()
        tmp_data["atom_numbs"] = [1, 1, 1, 2]
        tmp_data["atom_names"] = ["C", "H", "A", "B"]
        tmp_data["atom_types"] = np.array([0, 1, 2, 3, 3])
        # C1H1A1B2
        system_1_modified_type_1 = dpdata.LabeledSystem(data=tmp_data)

        tmp_data = system_1.data.copy()
        tmp_data["atom_numbs"] = [1, 1, 2, 1]
        tmp_data["atom_names"] = ["C", "H", "A", "B"]
        tmp_data["atom_types"] = np.array([0, 1, 2, 2, 3])
        # C1H1A2B1
        system_1_modified_type_2 = dpdata.LabeledSystem(data=tmp_data)

        tmp_data = system_1.data.copy()
        tmp_data["atom_numbs"] = [1, 1, 1, 2]
        tmp_data["atom_names"] = ["C", "H", "A", "D"]
        tmp_data["atom_types"] = np.array([0, 1, 2, 3, 3])
        # C1H1A1C2
        system_1_modified_type_3 = dpdata.LabeledSystem(data=tmp_data)

        self.ms = dpdata.MultiSystems(
            system_1,
            system_2,
            system_1_modified_type_1,
            system_1_modified_type_2,
            system_1_modified_type_3,
        )

        self.ms.to_deepmd_npy_mixed("tmp.deepmd.mixed")
        self.place_holder_ms = dpdata.MultiSystems()
        self.place_holder_ms.from_deepmd_npy("tmp.deepmd.mixed", fmt="deepmd/npy")

        new_type_map = ["H", "C", "D", "A", "B"]
        self.systems = dpdata.MultiSystems()
        self.systems.from_deepmd_npy_mixed(
            "tmp.deepmd.mixed", fmt="deepmd/npy/mixed", type_map=new_type_map
        )
        for kk in [ii.formula for ii in self.ms]:
            # apply type_map to each system
            self.ms[kk].apply_type_map(new_type_map)
            # revise keys in dict according because the type_map is updated.
            tmp_ss = self.ms.systems.pop(kk)
            self.ms.systems[tmp_ss.formula] = tmp_ss

        self.ms_1 = self.ms
        self.ms_2 = self.systems
        mixed_sets = glob("tmp.deepmd.mixed/*/set.*")
        self.assertEqual(len(mixed_sets), 2)
        for i in mixed_sets:
            self.assertEqual(
                os.path.exists(os.path.join(i, "real_atom_types.npy")), True
            )

        self.system_names = [
            "H4C1D0A0B0",
            "H3C1D0A0B0",
            "H1C1D0A1B2",
            "H1C1D0A2B1",
            "H1C1D2A1B0",
        ]
        self.system_sizes = {
            "H4C1D0A0B0": 1,
            "H3C1D0A0B0": 1,
            "H1C1D0A1B2": 1,
            "H1C1D0A2B1": 1,
            "H1C1D2A1B0": 1,
        }
        self.atom_names = ["H", "C", "D", "A", "B"]

    def tearDown(self):
        if os.path.exists("tmp.deepmd.mixed"):
            shutil.rmtree("tmp.deepmd.mixed")

    def test_len(self):
        self.assertEqual(len(self.ms), 5)
        self.assertEqual(len(self.place_holder_ms), 2)
        self.assertEqual(len(self.systems), 5)

    def test_get_nframes(self):
        self.assertEqual(self.ms.get_nframes(), 5)
        self.assertEqual(self.place_holder_ms.get_nframes(), 5)
        self.assertEqual(self.systems.get_nframes(), 5)

    def test_str(self):
        self.assertEqual(str(self.ms), "MultiSystems (5 systems containing 5 frames)")
        self.assertEqual(
            str(self.place_holder_ms), "MultiSystems (2 systems containing 5 frames)"
        )
        self.assertEqual(
            str(self.systems), "MultiSystems (5 systems containing 5 frames)"
        )


class TestMixedMultiSystemsDumpLoadSetSize(
    unittest.TestCase, CompLabeledMultiSys, MultiSystems, MSAllIsNoPBC
):
    def setUp(self):
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

        # C1H4
        system_1 = dpdata.LabeledSystem(
            "gaussian/methane.gaussianlog", fmt="gaussian/log"
        )

        # C1H3
        system_2 = dpdata.LabeledSystem(
            "gaussian/methane_sub.gaussianlog", fmt="gaussian/log"
        )

        tmp_data = system_1.data.copy()
        tmp_data["atom_numbs"] = [1, 1, 1, 2]
        tmp_data["atom_names"] = ["C", "H", "A", "B"]
        tmp_data["atom_types"] = np.array([0, 1, 2, 3, 3])
        # C1H1A1B2
        system_1_modified_type_1 = dpdata.LabeledSystem(data=tmp_data)

        tmp_data = system_1.data.copy()
        tmp_data["atom_numbs"] = [1, 1, 2, 1]
        tmp_data["atom_names"] = ["C", "H", "A", "B"]
        tmp_data["atom_types"] = np.array([0, 1, 2, 2, 3])
        # C1H1A2B1
        system_1_modified_type_2 = dpdata.LabeledSystem(data=tmp_data)

        tmp_data = system_1.data.copy()
        tmp_data["atom_numbs"] = [1, 1, 1, 2]
        tmp_data["atom_names"] = ["C", "H", "A", "D"]
        tmp_data["atom_types"] = np.array([0, 1, 2, 3, 3])
        # C1H1A1C2
        system_1_modified_type_3 = dpdata.LabeledSystem(data=tmp_data)

        self.ms = dpdata.MultiSystems(
            system_1,
            system_2,
            system_1_modified_type_1,
            system_1_modified_type_2,
            system_1_modified_type_3,
        )
        self.ms.to_deepmd_npy_mixed("tmp.deepmd.mixed", set_size=1)
        self.place_holder_ms = dpdata.MultiSystems()
        self.place_holder_ms.from_deepmd_npy("tmp.deepmd.mixed", fmt="deepmd/npy")
        self.systems = dpdata.MultiSystems()
        self.systems.from_deepmd_npy_mixed("tmp.deepmd.mixed", fmt="deepmd/npy/mixed")
        self.ms_1 = self.ms
        self.ms_2 = self.systems
        mixed_sets = glob("tmp.deepmd.mixed/*/set.*")
        self.assertEqual(len(mixed_sets), 5)
        for i in mixed_sets:
            self.assertEqual(
                os.path.exists(os.path.join(i, "real_atom_types.npy")), True
            )

        self.system_names = [
            "C1H4A0B0D0",
            "C1H3A0B0D0",
            "C1H1A1B2D0",
            "C1H1A2B1D0",
            "C1H1A1B0D2",
        ]
        self.system_sizes = {
            "C1H4A0B0D0": 1,
            "C1H3A0B0D0": 1,
            "C1H1A1B2D0": 1,
            "C1H1A2B1D0": 1,
            "C1H1A1B0D2": 1,
        }
        self.atom_names = ["C", "H", "A", "B", "D"]

    def tearDown(self):
        if os.path.exists("tmp.deepmd.mixed"):
            shutil.rmtree("tmp.deepmd.mixed")

    def test_len(self):
        self.assertEqual(len(self.ms), 5)
        self.assertEqual(len(self.place_holder_ms), 2)
        self.assertEqual(len(self.systems), 5)

    def test_get_nframes(self):
        self.assertEqual(self.ms.get_nframes(), 5)
        self.assertEqual(self.place_holder_ms.get_nframes(), 5)
        self.assertEqual(self.systems.get_nframes(), 5)

    def test_str(self):
        self.assertEqual(str(self.ms), "MultiSystems (5 systems containing 5 frames)")
        self.assertEqual(
            str(self.place_holder_ms), "MultiSystems (2 systems containing 5 frames)"
        )
        self.assertEqual(
            str(self.systems), "MultiSystems (5 systems containing 5 frames)"
        )


class TestMixedMultiSystemsTypeChange(
    unittest.TestCase, CompLabeledMultiSys, MultiSystems, MSAllIsNoPBC
):
    def setUp(self):
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

        # C1H4
        system_1 = dpdata.LabeledSystem(
            "gaussian/methane.gaussianlog", fmt="gaussian/log"
        )

        # C1H3
        system_2 = dpdata.LabeledSystem(
            "gaussian/methane_sub.gaussianlog", fmt="gaussian/log"
        )

        tmp_data = system_1.data.copy()
        tmp_data["atom_numbs"] = [1, 1, 1, 2]
        tmp_data["atom_names"] = ["C", "H", "A", "B"]
        tmp_data["atom_types"] = np.array([0, 1, 2, 3, 3])
        # C1H1A1B2
        system_1_modified_type_1 = dpdata.LabeledSystem(data=tmp_data)

        tmp_data = system_1.data.copy()
        tmp_data["atom_numbs"] = [1, 1, 2, 1]
        tmp_data["atom_names"] = ["C", "H", "A", "B"]
        tmp_data["atom_types"] = np.array([0, 1, 2, 2, 3])
        # C1H1A2B1
        system_1_modified_type_2 = dpdata.LabeledSystem(data=tmp_data)

        tmp_data = system_1.data.copy()
        tmp_data["atom_numbs"] = [1, 1, 1, 2]
        tmp_data["atom_names"] = ["C", "H", "A", "D"]
        tmp_data["atom_types"] = np.array([0, 1, 2, 3, 3])
        # C1H1A1C2
        system_1_modified_type_3 = dpdata.LabeledSystem(data=tmp_data)

        self.ms = dpdata.MultiSystems(
            system_1,
            system_2,
            system_1_modified_type_1,
            system_1_modified_type_2,
            system_1_modified_type_3,
            type_map=["TOKEN"],
        )
        self.ms.to_deepmd_npy_mixed("tmp.deepmd.mixed")
        self.place_holder_ms = dpdata.MultiSystems()
        self.place_holder_ms.from_deepmd_npy("tmp.deepmd.mixed", fmt="deepmd/npy")
        self.systems = dpdata.MultiSystems(type_map=["TOKEN"])
        self.systems.from_deepmd_npy_mixed("tmp.deepmd.mixed", fmt="deepmd/npy/mixed")
        self.ms_1 = self.ms
        self.ms_2 = self.systems
        mixed_sets = glob("tmp.deepmd.mixed/*/set.*")
        self.assertEqual(len(mixed_sets), 2)
        for i in mixed_sets:
            self.assertEqual(
                os.path.exists(os.path.join(i, "real_atom_types.npy")), True
            )

        self.system_names = [
            "TOKEN0C1H4A0B0D0",
            "TOKEN0C1H3A0B0D0",
            "TOKEN0C1H1A1B2D0",
            "TOKEN0C1H1A2B1D0",
            "TOKEN0C1H1A1B0D2",
        ]
        self.system_sizes = {
            "TOKEN0C1H4A0B0D0": 1,
            "TOKEN0C1H3A0B0D0": 1,
            "TOKEN0C1H1A1B2D0": 1,
            "TOKEN0C1H1A2B1D0": 1,
            "TOKEN0C1H1A1B0D2": 1,
        }
        self.atom_names = ["C", "H", "A", "B", "D"]

    def tearDown(self):
        if os.path.exists("tmp.deepmd.mixed"):
            shutil.rmtree("tmp.deepmd.mixed")

    def test_len(self):
        self.assertEqual(len(self.ms), 5)
        self.assertEqual(len(self.place_holder_ms), 2)
        self.assertEqual(len(self.systems), 5)

    def test_get_nframes(self):
        self.assertEqual(self.ms.get_nframes(), 5)
        self.assertEqual(self.place_holder_ms.get_nframes(), 5)
        self.assertEqual(self.systems.get_nframes(), 5)

    def test_str(self):
        self.assertEqual(str(self.ms), "MultiSystems (5 systems containing 5 frames)")
        self.assertEqual(
            str(self.place_holder_ms), "MultiSystems (2 systems containing 5 frames)"
        )
        self.assertEqual(
            str(self.systems), "MultiSystems (5 systems containing 5 frames)"
        )


class TestMixedSingleSystemsDump(unittest.TestCase, CompLabeledSys, IsNoPBC):
    def setUp(self):
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

        # C1H4
        self.system_1 = dpdata.LabeledSystem(
            "gaussian/methane.gaussianlog", fmt="gaussian/log"
        )
        self.system_2 = dpdata.LabeledSystem(
            "gaussian/methane.gaussianlog", fmt="gaussian/log"
        )
        # test dump
        self.system_1.to("deepmd/npy/mixed", "tmp.deepmd.mixed.single")

    def tearDown(self):
        if os.path.exists("tmp.deepmd.mixed.single"):
            shutil.rmtree("tmp.deepmd.mixed.single")


if __name__ == "__main__":
    unittest.main()
