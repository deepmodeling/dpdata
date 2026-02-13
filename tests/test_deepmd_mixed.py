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

from dpdata.data_type import (
    Axis,
    DataType,
)


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


class TestMixedSystemWithFparamAparam(
    unittest.TestCase, CompLabeledMultiSys, MultiSystems, MSAllIsNoPBC
):
    def setUp(self):
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

        new_datatypes = [
            DataType(
                "fparam",
                np.ndarray,
                shape=(Axis.NFRAMES, 2),
                required=False,
            ),
            DataType(
                "aparam",
                np.ndarray,
                shape=(Axis.NFRAMES, Axis.NATOMS, 3),
                required=False,
            ),
        ]

        for datatype in new_datatypes:
            dpdata.System.register_data_type(datatype)
            dpdata.LabeledSystem.register_data_type(datatype)

        # C1H4
        system_1 = dpdata.LabeledSystem(
            "gaussian/methane.gaussianlog", fmt="gaussian/log"
        )

        # C1H3
        system_2 = dpdata.LabeledSystem(
            "gaussian/methane_sub.gaussianlog", fmt="gaussian/log"
        )

        tmp_data_1 = system_1.data.copy()
        nframes_1 = tmp_data_1["coords"].shape[0]
        natoms_1 = tmp_data_1["atom_types"].shape[0]
        tmp_data_1["fparam"] = np.random.random([nframes_1, 2])
        tmp_data_1["aparam"] = np.random.random([nframes_1, natoms_1, 3])
        system_1_with_params = dpdata.LabeledSystem(data=tmp_data_1)

        tmp_data_2 = system_2.data.copy()
        nframes_2 = tmp_data_2["coords"].shape[0]
        natoms_2 = tmp_data_2["atom_types"].shape[0]
        tmp_data_2["fparam"] = np.random.random([nframes_2, 2])
        tmp_data_2["aparam"] = np.random.random([nframes_2, natoms_2, 3])
        system_2_with_params = dpdata.LabeledSystem(data=tmp_data_2)

        tmp_data_3 = system_1.data.copy()
        nframes_3 = tmp_data_3["coords"].shape[0]
        tmp_data_3["atom_numbs"] = [1, 1, 1, 2]
        tmp_data_3["atom_names"] = ["C", "H", "A", "B"]
        tmp_data_3["atom_types"] = np.array([0, 1, 2, 3, 3])
        natoms_3 = len(tmp_data_3["atom_types"])
        tmp_data_3["fparam"] = np.random.random([nframes_3, 2])
        tmp_data_3["aparam"] = np.random.random([nframes_3, natoms_3, 3])
        # C1H1A1B2 with params
        system_3_with_params = dpdata.LabeledSystem(data=tmp_data_3)

        self.ms = dpdata.MultiSystems(
            system_1_with_params, system_2_with_params, system_3_with_params
        )

        self.ms.to_deepmd_npy_mixed("tmp.deepmd.fparam.aparam")
        self.place_holder_ms = dpdata.MultiSystems()
        self.place_holder_ms.from_deepmd_npy(
            "tmp.deepmd.fparam.aparam", fmt="deepmd/npy"
        )
        self.systems = dpdata.MultiSystems()
        self.systems.from_deepmd_npy_mixed(
            "tmp.deepmd.fparam.aparam", fmt="deepmd/npy/mixed"
        )

        self.ms_1 = self.ms
        self.ms_2 = self.systems

        mixed_sets = glob("tmp.deepmd.fparam.aparam/*/set.*")
        for i in mixed_sets:
            self.assertEqual(
                os.path.exists(os.path.join(i, "real_atom_types.npy")), True
            )

        self.system_names = ["C1H4A0B0", "C1H3A0B0", "C1H1A1B2"]
        self.system_sizes = {"C1H4A0B0": 1, "C1H3A0B0": 1, "C1H1A1B2": 1}
        self.atom_names = ["C", "H", "A", "B"]

    def tearDown(self):
        if os.path.exists("tmp.deepmd.fparam.aparam"):
            shutil.rmtree("tmp.deepmd.fparam.aparam")

    def test_len(self):
        self.assertEqual(len(self.ms), 3)
        self.assertEqual(len(self.systems), 3)

    def test_get_nframes(self):
        self.assertEqual(self.ms.get_nframes(), 3)
        self.assertEqual(self.systems.get_nframes(), 3)

    def test_str(self):
        self.assertEqual(str(self.ms), "MultiSystems (3 systems containing 3 frames)")
        self.assertEqual(
            str(self.systems), "MultiSystems (3 systems containing 3 frames)"
        )

    def test_fparam_exists(self):
        for formula in self.system_names:
            if formula in self.ms.systems:
                self.assertTrue("fparam" in self.ms[formula].data)
            if formula in self.systems.systems:
                self.assertTrue("fparam" in self.systems[formula].data)

        for formula in self.system_names:
            if formula in self.ms.systems and formula in self.systems.systems:
                np.testing.assert_almost_equal(
                    self.ms[formula].data["fparam"],
                    self.systems[formula].data["fparam"],
                    decimal=self.places,
                )

    def test_aparam_exists(self):
        for formula in self.system_names:
            if formula in self.ms.systems:
                self.assertTrue("aparam" in self.ms[formula].data)
            if formula in self.systems.systems:
                self.assertTrue("aparam" in self.systems[formula].data)

        for formula in self.system_names:
            if formula in self.ms.systems and formula in self.systems.systems:
                np.testing.assert_almost_equal(
                    self.ms[formula].data["aparam"],
                    self.systems[formula].data["aparam"],
                    decimal=self.places,
                )


class TestMixedMultiSystemsPadding(
    unittest.TestCase, CompLabeledMultiSys, MultiSystems, MSAllIsNoPBC
):
    """Test round-trip with atom_numb_pad.

    C1H4 (5 atoms) and C1H3 (4 atoms) are both padded to 8 atoms,
    so only 1 subfolder should be created.
    """

    def setUp(self):
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

        # C1H4 (5 atoms)
        system_1 = dpdata.LabeledSystem(
            "gaussian/methane.gaussianlog", fmt="gaussian/log"
        )
        # C1H3 (4 atoms)
        system_2 = dpdata.LabeledSystem(
            "gaussian/methane_sub.gaussianlog", fmt="gaussian/log"
        )

        self.ms = dpdata.MultiSystems(system_1, system_2)
        self.ms.to_deepmd_npy_mixed("tmp.deepmd.mixed.pad", atom_numb_pad=8)
        self.systems = dpdata.MultiSystems()
        self.systems.from_deepmd_npy_mixed(
            "tmp.deepmd.mixed.pad", fmt="deepmd/npy/mixed"
        )
        self.ms_1 = self.ms
        self.ms_2 = self.systems

        self.system_names = ["C1H4", "C1H3"]
        self.system_sizes = {"C1H4": 1, "C1H3": 1}
        self.atom_names = ["C", "H"]

    def tearDown(self):
        if os.path.exists("tmp.deepmd.mixed.pad"):
            shutil.rmtree("tmp.deepmd.mixed.pad")

    def test_single_subfolder(self):
        """Both 4-atom and 5-atom systems padded to 8 -> 1 subfolder."""
        subdirs = [
            d
            for d in os.listdir("tmp.deepmd.mixed.pad")
            if os.path.isdir(os.path.join("tmp.deepmd.mixed.pad", d))
        ]
        self.assertEqual(len(subdirs), 1)
        self.assertEqual(subdirs[0], "8")

    def test_real_atom_types_on_disk(self):
        """Verify real_atom_types.npy contains -1 for virtual atoms."""
        mixed_sets = glob("tmp.deepmd.mixed.pad/*/set.*")
        for s in mixed_sets:
            rat = np.load(os.path.join(s, "real_atom_types.npy"))
            # padded to 8, so last columns should be -1
            self.assertTrue(np.any(rat == -1))
            # first columns should be >= 0
            self.assertTrue(np.all(rat[:, 0] >= 0))

    def test_loaded_natoms(self):
        """Loaded systems should have original (unpadded) atom counts."""
        for formula, sys in self.systems.systems.items():
            if "H4" in formula:
                self.assertEqual(sys.get_natoms(), 5)
            elif "H3" in formula:
                self.assertEqual(sys.get_natoms(), 4)

    def test_len(self):
        self.assertEqual(len(self.ms), 2)
        self.assertEqual(len(self.systems), 2)

    def test_get_nframes(self):
        self.assertEqual(self.ms.get_nframes(), 2)
        self.assertEqual(self.systems.get_nframes(), 2)


class TestMixedMultiSystemsPaddingMultipleGroups(
    unittest.TestCase, CompLabeledMultiSys, MultiSystems, MSAllIsNoPBC
):
    """Test padding with systems that span multiple padded groups.

    With atom_numb_pad=4: C1H3 (4 atoms) -> 4, C1H4 (5 atoms) -> 8.
    Two subfolders should be created.
    """

    def setUp(self):
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

        # C1H4 (5 atoms)
        system_1 = dpdata.LabeledSystem(
            "gaussian/methane.gaussianlog", fmt="gaussian/log"
        )
        # C1H3 (4 atoms)
        system_2 = dpdata.LabeledSystem(
            "gaussian/methane_sub.gaussianlog", fmt="gaussian/log"
        )

        self.ms = dpdata.MultiSystems(system_1, system_2)
        self.ms.to_deepmd_npy_mixed("tmp.deepmd.mixed.pad2", atom_numb_pad=4)
        self.systems = dpdata.MultiSystems()
        self.systems.from_deepmd_npy_mixed(
            "tmp.deepmd.mixed.pad2", fmt="deepmd/npy/mixed"
        )
        self.ms_1 = self.ms
        self.ms_2 = self.systems

        self.system_names = ["C1H4", "C1H3"]
        self.system_sizes = {"C1H4": 1, "C1H3": 1}
        self.atom_names = ["C", "H"]

    def tearDown(self):
        if os.path.exists("tmp.deepmd.mixed.pad2"):
            shutil.rmtree("tmp.deepmd.mixed.pad2")

    def test_two_subfolders(self):
        """4-atom -> 4, 5-atom -> 8 => 2 subfolders."""
        subdirs = sorted(
            d
            for d in os.listdir("tmp.deepmd.mixed.pad2")
            if os.path.isdir(os.path.join("tmp.deepmd.mixed.pad2", d))
        )
        self.assertEqual(len(subdirs), 2)
        self.assertIn("4", subdirs)
        self.assertIn("8", subdirs)

    def test_len(self):
        self.assertEqual(len(self.ms), 2)
        self.assertEqual(len(self.systems), 2)

    def test_get_nframes(self):
        self.assertEqual(self.ms.get_nframes(), 2)
        self.assertEqual(self.systems.get_nframes(), 2)


class TestMixedMultiSystemsPaddingTypeMap(
    unittest.TestCase, CompLabeledMultiSys, MSAllIsNoPBC
):
    """Test padding + custom type_map on reload.

    This verifies the index_map bug fix for -1 values in real_atom_types.
    """

    def setUp(self):
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

        # C1H4 (5 atoms)
        system_1 = dpdata.LabeledSystem(
            "gaussian/methane.gaussianlog", fmt="gaussian/log"
        )
        # C1H3 (4 atoms)
        system_2 = dpdata.LabeledSystem(
            "gaussian/methane_sub.gaussianlog", fmt="gaussian/log"
        )

        self.ms = dpdata.MultiSystems(system_1, system_2)
        self.ms.to_deepmd_npy_mixed("tmp.deepmd.mixed.pad.tm", atom_numb_pad=8)

        new_type_map = ["H", "C"]
        self.systems = dpdata.MultiSystems()
        self.systems.from_deepmd_npy_mixed(
            "tmp.deepmd.mixed.pad.tm",
            fmt="deepmd/npy/mixed",
            type_map=new_type_map,
        )

        # Apply same type_map to original for comparison
        for kk in [ii.formula for ii in self.ms]:
            self.ms[kk].apply_type_map(new_type_map)
            tmp_ss = self.ms.systems.pop(kk)
            self.ms.systems[tmp_ss.formula] = tmp_ss

        self.ms_1 = self.ms
        self.ms_2 = self.systems

    def tearDown(self):
        if os.path.exists("tmp.deepmd.mixed.pad.tm"):
            shutil.rmtree("tmp.deepmd.mixed.pad.tm")

    def test_len(self):
        self.assertEqual(len(self.ms), 2)
        self.assertEqual(len(self.systems), 2)

    def test_get_nframes(self):
        self.assertEqual(self.ms.get_nframes(), 2)
        self.assertEqual(self.systems.get_nframes(), 2)


class TestMixedMultiSystemsPaddingAparam(
    unittest.TestCase, CompLabeledMultiSys, MultiSystems, MSAllIsNoPBC
):
    """Test padding with custom per-atom data (aparam)."""

    def setUp(self):
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

        new_datatypes = [
            DataType(
                "fparam",
                np.ndarray,
                shape=(Axis.NFRAMES, 2),
                required=False,
            ),
            DataType(
                "aparam",
                np.ndarray,
                shape=(Axis.NFRAMES, Axis.NATOMS, 3),
                required=False,
            ),
        ]
        for datatype in new_datatypes:
            dpdata.System.register_data_type(datatype)
            dpdata.LabeledSystem.register_data_type(datatype)

        # C1H4 (5 atoms)
        system_1 = dpdata.LabeledSystem(
            "gaussian/methane.gaussianlog", fmt="gaussian/log"
        )
        # C1H3 (4 atoms)
        system_2 = dpdata.LabeledSystem(
            "gaussian/methane_sub.gaussianlog", fmt="gaussian/log"
        )

        tmp_data_1 = system_1.data.copy()
        nframes_1 = tmp_data_1["coords"].shape[0]
        natoms_1 = tmp_data_1["atom_types"].shape[0]
        tmp_data_1["fparam"] = np.random.random([nframes_1, 2])
        tmp_data_1["aparam"] = np.random.random([nframes_1, natoms_1, 3])
        system_1_with_params = dpdata.LabeledSystem(data=tmp_data_1)

        tmp_data_2 = system_2.data.copy()
        nframes_2 = tmp_data_2["coords"].shape[0]
        natoms_2 = tmp_data_2["atom_types"].shape[0]
        tmp_data_2["fparam"] = np.random.random([nframes_2, 2])
        tmp_data_2["aparam"] = np.random.random([nframes_2, natoms_2, 3])
        system_2_with_params = dpdata.LabeledSystem(data=tmp_data_2)

        self.ms = dpdata.MultiSystems(system_1_with_params, system_2_with_params)
        self.ms.to_deepmd_npy_mixed("tmp.deepmd.mixed.pad.ap", atom_numb_pad=8)
        self.systems = dpdata.MultiSystems()
        self.systems.from_deepmd_npy_mixed(
            "tmp.deepmd.mixed.pad.ap", fmt="deepmd/npy/mixed"
        )
        self.ms_1 = self.ms
        self.ms_2 = self.systems

        self.system_names = ["C1H4", "C1H3"]
        self.system_sizes = {"C1H4": 1, "C1H3": 1}
        self.atom_names = ["C", "H"]

    def tearDown(self):
        if os.path.exists("tmp.deepmd.mixed.pad.ap"):
            shutil.rmtree("tmp.deepmd.mixed.pad.ap")

    def test_single_subfolder(self):
        subdirs = [
            d
            for d in os.listdir("tmp.deepmd.mixed.pad.ap")
            if os.path.isdir(os.path.join("tmp.deepmd.mixed.pad.ap", d))
        ]
        self.assertEqual(len(subdirs), 1)

    def test_fparam_preserved(self):
        for formula in self.system_names:
            if formula in self.ms.systems and formula in self.systems.systems:
                np.testing.assert_almost_equal(
                    self.ms[formula].data["fparam"],
                    self.systems[formula].data["fparam"],
                    decimal=self.places,
                )

    def test_aparam_preserved(self):
        """Per-atom aparam should be correctly padded and unpadded."""
        for formula in self.system_names:
            if formula in self.ms.systems and formula in self.systems.systems:
                np.testing.assert_almost_equal(
                    self.ms[formula].data["aparam"],
                    self.systems[formula].data["aparam"],
                    decimal=self.places,
                )

    def test_len(self):
        self.assertEqual(len(self.ms), 2)
        self.assertEqual(len(self.systems), 2)

    def test_get_nframes(self):
        self.assertEqual(self.ms.get_nframes(), 2)
        self.assertEqual(self.systems.get_nframes(), 2)
