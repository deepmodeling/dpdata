from __future__ import annotations

import os
import unittest

import h5py  # noqa: TID253
import numpy as np
from comp_sys import (
    CompLabeledMultiSys,
    CompLabeledSys,
    CompSys,
    IsNoPBC,
    IsPBC,
    MSAllIsNoPBC,
    MultiSystems,
)
from context import dpdata


class TestDeepmdLoadDumpHDF5(unittest.TestCase, CompLabeledSys, IsPBC):
    def setUp(self):
        self.system_1 = dpdata.LabeledSystem("poscars/OUTCAR.h2o.md", fmt="vasp/outcar")
        self.system_1.to_deepmd_hdf5("tmp.deepmd.hdf5", prec=np.float64, set_size=2)

        self.system_2 = dpdata.LabeledSystem(
            "tmp.deepmd.hdf5", fmt="deepmd/hdf5", type_map=["O", "H"]
        )
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

    def tearDown(self):
        if os.path.exists("tmp.deepmd.hdf5"):
            os.remove("tmp.deepmd.hdf5")


class TestDeepmdHDF5NoLabels(unittest.TestCase, CompSys, IsPBC):
    def setUp(self):
        self.system_1 = dpdata.System("poscars/POSCAR.h2o.md", fmt="vasp/poscar")
        self.system_1.to_deepmd_hdf5("tmp.deepmd.hdf5", prec=np.float64, set_size=2)
        self.system_2 = dpdata.System(
            "tmp.deepmd.hdf5", fmt="deepmd/hdf5", type_map=["O", "H"]
        )
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

    def tearDown(self):
        if os.path.exists("tmp.deepmd.hdf5"):
            os.remove("tmp.deepmd.hdf5")


class TestHDF5Multi(unittest.TestCase, CompLabeledSys, MultiSystems, IsNoPBC):
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
        systems = dpdata.MultiSystems(system_1, system_2, system_3)
        systems.to_deepmd_hdf5("tmp.deepmd.hdf5")

        self.systems = dpdata.MultiSystems().from_deepmd_hdf5("tmp.deepmd.hdf5")
        self.system_names = ["C1H4", "C1H3"]
        self.system_sizes = {"C1H4": 2, "C1H3": 1}
        self.atom_names = ["C", "H"]
        self.system_1 = self.systems["C1H3"]
        self.system_2 = system_3

    def tearDown(self):
        if os.path.exists("tmp.deepmd.hdf5"):
            os.remove("tmp.deepmd.hdf5")


class TestHDF5MixedMulti(
    unittest.TestCase, CompLabeledMultiSys, MultiSystems, MSAllIsNoPBC
):
    def setUp(self):
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

        system_1 = dpdata.LabeledSystem(
            "gaussian/methane.gaussianlog", fmt="gaussian/log"
        )
        system_2 = dpdata.LabeledSystem(
            "gaussian/methane_sub.gaussianlog", fmt="gaussian/log"
        )

        tmp_data = system_1.data.copy()
        tmp_data["atom_numbs"] = [1, 1, 1, 2]
        tmp_data["atom_names"] = ["C", "H", "A", "B"]
        tmp_data["atom_types"] = np.array([0, 1, 2, 3, 3])
        system_3 = dpdata.LabeledSystem(data=tmp_data)

        self.ms = dpdata.MultiSystems(system_1, system_2, system_3)
        self.ms.to_deepmd_hdf5_mixed("tmp.deepmd.mixed.hdf5")
        self.systems = dpdata.MultiSystems().from_deepmd_hdf5_mixed(
            "tmp.deepmd.mixed.hdf5"
        )
        self.ms_1 = self.ms
        self.ms_2 = self.systems

        self.system_names = ["C1H4A0B0", "C1H3A0B0", "C1H1A1B2"]
        self.system_sizes = {"C1H4A0B0": 1, "C1H3A0B0": 1, "C1H1A1B2": 1}
        self.atom_names = ["C", "H", "A", "B"]

    def tearDown(self):
        if os.path.exists("tmp.deepmd.mixed.hdf5"):
            os.remove("tmp.deepmd.mixed.hdf5")

    def test_hdf5_group_layout(self):
        with h5py.File("tmp.deepmd.mixed.hdf5", "r") as f:
            self.assertEqual(set(f.keys()), {"4", "5"})
            for group in f.values():
                self.assertIn("type_map.raw", group)
                self.assertIn("set.000/real_atom_types.npy", group)


class TestHDF5MixedPadding(
    unittest.TestCase, CompLabeledMultiSys, MultiSystems, MSAllIsNoPBC
):
    def setUp(self):
        self.places = 6
        self.e_places = 6
        self.f_places = 6
        self.v_places = 6

        system_1 = dpdata.LabeledSystem(
            "gaussian/methane.gaussianlog", fmt="gaussian/log"
        )
        system_2 = dpdata.LabeledSystem(
            "gaussian/methane_sub.gaussianlog", fmt="gaussian/log"
        )

        self.ms = dpdata.MultiSystems(system_1, system_2)
        self.ms.to_deepmd_hdf5_mixed("tmp.deepmd.mixed.pad.hdf5", atom_numb_pad=8)
        self.systems = dpdata.MultiSystems().from_deepmd_hdf5_mixed(
            "tmp.deepmd.mixed.pad.hdf5"
        )
        self.ms_1 = self.ms
        self.ms_2 = self.systems

        self.system_names = ["C1H4", "C1H3"]
        self.system_sizes = {"C1H4": 1, "C1H3": 1}
        self.atom_names = ["C", "H"]

    def tearDown(self):
        if os.path.exists("tmp.deepmd.mixed.pad.hdf5"):
            os.remove("tmp.deepmd.mixed.pad.hdf5")

    def test_single_padded_group(self):
        with h5py.File("tmp.deepmd.mixed.pad.hdf5", "r") as f:
            self.assertEqual(list(f.keys()), ["8"])
            real_atom_types = f["8/set.000/real_atom_types.npy"][:]
            self.assertEqual(real_atom_types.shape[1], 8)
            self.assertTrue(np.any(real_atom_types == -1))


class TestHDF5MixedIOVariants(unittest.TestCase):
    def tearDown(self):
        for file_name in (
            "tmp.deepmd.mixed.single.hdf5",
            "tmp.deepmd.mixed.group.hdf5",
            "tmp.deepmd.mixed.object.hdf5",
            "tmp.deepmd.mixed.unlabeled.hdf5",
            "tmp.deepmd.mixed.typemap.hdf5",
            "tmp.deepmd.regular.hdf5",
        ):
            if os.path.exists(file_name):
                os.remove(file_name)

    def test_single_system_string_round_trip(self):
        system = dpdata.LabeledSystem(
            "gaussian/methane.gaussianlog", fmt="gaussian/log"
        )
        system.to("deepmd/hdf5/mixed", "tmp.deepmd.mixed.single.hdf5")

        systems = dpdata.MultiSystems().from_deepmd_hdf5_mixed(
            "tmp.deepmd.mixed.single.hdf5"
        )

        self.assertEqual(len(systems), 1)
        self.assertIn("C1H4", systems.systems)
        np.testing.assert_allclose(
            systems["C1H4"].data["coords"], system.data["coords"]
        )

    def test_hash_group_round_trip(self):
        system_1 = dpdata.LabeledSystem(
            "gaussian/methane.gaussianlog", fmt="gaussian/log"
        )
        system_2 = dpdata.LabeledSystem(
            "gaussian/methane_sub.gaussianlog", fmt="gaussian/log"
        )

        dpdata.MultiSystems(system_1, system_2).to_deepmd_hdf5_mixed(
            "tmp.deepmd.mixed.group.hdf5#mixed"
        )
        systems = dpdata.MultiSystems().from_deepmd_hdf5_mixed(
            "tmp.deepmd.mixed.group.hdf5#mixed"
        )

        self.assertEqual(set(systems.systems), {"C1H4", "C1H3"})
        with h5py.File("tmp.deepmd.mixed.group.hdf5", "r") as f:
            self.assertEqual(set(f["mixed"].keys()), {"4", "5"})

    def test_hdf5_object_round_trip(self):
        system_1 = dpdata.LabeledSystem(
            "gaussian/methane.gaussianlog", fmt="gaussian/log"
        )
        system_2 = dpdata.LabeledSystem(
            "gaussian/methane_sub.gaussianlog", fmt="gaussian/log"
        )

        with h5py.File("tmp.deepmd.mixed.object.hdf5", "w") as f:
            f.create_group("5")
            dpdata.MultiSystems(system_1, system_2).to_deepmd_hdf5_mixed(f)

        with h5py.File("tmp.deepmd.mixed.object.hdf5", "r") as f:
            systems = dpdata.MultiSystems().from_deepmd_hdf5_mixed(f)

        self.assertEqual(set(systems.systems), {"C1H4", "C1H3"})
        np.testing.assert_allclose(
            systems["C1H4"].data["forces"], system_1.data["forces"]
        )

    def test_unlabeled_round_trip(self):
        system = dpdata.System("poscars/POSCAR.h2o.md", fmt="vasp/poscar")
        system.to("deepmd/hdf5/mixed", "tmp.deepmd.mixed.unlabeled.hdf5")

        systems = dpdata.MultiSystems().load_systems_from_file(
            "tmp.deepmd.mixed.unlabeled.hdf5",
            fmt="deepmd/hdf5/mixed",
            labeled=False,
        )

        self.assertEqual(len(systems), 1)
        self.assertNotIn("energies", list(systems.systems.values())[0].data)
        np.testing.assert_allclose(
            list(systems.systems.values())[0].data["coords"], system.data["coords"]
        )

    def test_type_map_round_trip(self):
        system = dpdata.LabeledSystem(
            "gaussian/methane.gaussianlog", fmt="gaussian/log"
        )
        dpdata.MultiSystems(system).to_deepmd_hdf5_mixed(
            "tmp.deepmd.mixed.typemap.hdf5"
        )

        systems = dpdata.MultiSystems().from_deepmd_hdf5_mixed(
            "tmp.deepmd.mixed.typemap.hdf5", type_map=["H", "C"]
        )
        system_ref = system.copy()
        system_ref.apply_type_map(["H", "C"])

        self.assertEqual(set(systems.systems), {system_ref.formula})
        np.testing.assert_allclose(
            systems[system_ref.formula].data["forces"], system_ref.data["forces"]
        )

    def test_unsupported_inputs(self):
        fmt = dpdata.plugins.deepmd.DeePMDHDF5MixedFormat()

        with self.assertRaises(TypeError):
            fmt.from_system_mix(object())
        with self.assertRaises(TypeError):
            fmt.to_system({}, object())
        with self.assertRaises(TypeError):
            list(fmt.from_multi_systems(object()))
        with self.assertRaises(TypeError):
            list(fmt.to_multi_systems(["1"], object()))

    def test_regular_hdf5_groups_are_not_mixed(self):
        system_1 = dpdata.LabeledSystem(
            "gaussian/methane.gaussianlog", fmt="gaussian/log"
        )
        system_2 = dpdata.LabeledSystem(
            "gaussian/methane_sub.gaussianlog", fmt="gaussian/log"
        )
        dpdata.MultiSystems(system_1, system_2).to_deepmd_hdf5(
            "tmp.deepmd.regular.hdf5"
        )

        systems = dpdata.MultiSystems().from_deepmd_hdf5_mixed(
            "tmp.deepmd.regular.hdf5"
        )

        self.assertEqual(len(systems), 0)
