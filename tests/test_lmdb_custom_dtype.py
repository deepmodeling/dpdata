from __future__ import annotations

import os
import shutil
import unittest

import numpy as np

import dpdata
from dpdata.data_type import Axis, DataType
from dpdata.lmdb.format import LMDBFormat


class TestLMDBCustomDType(unittest.TestCase):
    def setUp(self):
        self.original_system_dtypes = dpdata.System.DTYPES
        self.original_labeled_system_dtypes = dpdata.LabeledSystem.DTYPES

        # Register custom data types as optional
        self.dt_frame = DataType(
            "frame_data", np.ndarray, shape=(Axis.NFRAMES, 2), required=False
        )
        self.dt_static = DataType("static_data", np.ndarray, shape=(2,), required=False)

        dpdata.System.register_data_type(self.dt_frame, self.dt_static)
        dpdata.LabeledSystem.register_data_type(self.dt_frame, self.dt_static)

        self.lmdb_path = "tmp_custom_dtype.lmdb"

        # Create a system with custom data
        # Assuming running from tests/ directory
        try:
            self.system = dpdata.LabeledSystem(
                "poscars/OUTCAR.h2o.md", fmt="vasp/outcar"
            )
        except FileNotFoundError:
            self.system = dpdata.LabeledSystem(
                "tests/poscars/OUTCAR.h2o.md", fmt="vasp/outcar"
            )

        nframes = self.system.get_nframes()
        self.system.data["frame_data"] = np.random.rand(nframes, 2)
        self.system.data["static_data"] = np.array([1.0, 2.0])
        self.system.check_data()

    def tearDown(self):
        dpdata.System.DTYPES = self.original_system_dtypes
        dpdata.LabeledSystem.DTYPES = self.original_labeled_system_dtypes
        if os.path.exists(self.lmdb_path):
            shutil.rmtree(self.lmdb_path)

    def test_custom_dtype_preservation(self):
        self.system.to("lmdb", self.lmdb_path)
        system_loaded = dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")

        np.testing.assert_allclose(
            system_loaded.data["frame_data"], self.system.data["frame_data"]
        )
        np.testing.assert_allclose(
            system_loaded.data["static_data"], self.system.data["static_data"]
        )

    def test_multi_systems_custom_dtype(self):
        ms = dpdata.MultiSystems(self.system)
        LMDBFormat().to_multi_systems(list(ms.systems.values()), self.lmdb_path)
        ms_loaded = LMDBFormat().from_multi_systems(self.lmdb_path)

        system_loaded = list(ms_loaded.systems.values())[0]
        np.testing.assert_allclose(
            system_loaded.data["frame_data"], self.system.data["frame_data"]
        )
        np.testing.assert_allclose(
            system_loaded.data["static_data"], self.system.data["static_data"]
        )

    def test_custom_dtype_auto_registration(self):
        # Save with custom data types registered
        self.system.to("lmdb", self.lmdb_path)

        # Simulate a clean session by unregistering the custom types
        dpdata.System.DTYPES = self.original_system_dtypes
        dpdata.LabeledSystem.DTYPES = self.original_labeled_system_dtypes

        # Verify they are currently missing
        self.assertNotIn("frame_data", [dt.name for dt in dpdata.LabeledSystem.DTYPES])

        # Load from LMDB - should trigger auto-registration
        system_loaded = dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")

        # Verify data is loaded and types are registered
        self.assertIn("frame_data", [dt.name for dt in dpdata.LabeledSystem.DTYPES])
        self.assertIn("static_data", [dt.name for dt in dpdata.LabeledSystem.DTYPES])

        np.testing.assert_allclose(
            system_loaded.data["frame_data"], self.system.data["frame_data"]
        )
        np.testing.assert_allclose(
            system_loaded.data["static_data"], self.system.data["static_data"]
        )


class TestLMDBFparamAparam(unittest.TestCase):
    def setUp(self):
        self.original_system_dtypes = dpdata.System.DTYPES
        self.original_labeled_system_dtypes = dpdata.LabeledSystem.DTYPES

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

        self.lmdb_path = "tmp_fparam_aparam.lmdb"

        try:
            self.system = dpdata.LabeledSystem(
                "poscars/OUTCAR.h2o.md", fmt="vasp/outcar"
            )
        except FileNotFoundError:
            self.system = dpdata.LabeledSystem(
                "tests/poscars/OUTCAR.h2o.md", fmt="vasp/outcar"
            )

        nframes = self.system.get_nframes()
        natoms = self.system.get_natoms()
        self.system.data["fparam"] = np.random.rand(nframes, 2)
        self.system.data["aparam"] = np.random.rand(nframes, natoms, 3)
        self.system.check_data()

    def tearDown(self):
        dpdata.System.DTYPES = self.original_system_dtypes
        dpdata.LabeledSystem.DTYPES = self.original_labeled_system_dtypes
        if os.path.exists(self.lmdb_path):
            shutil.rmtree(self.lmdb_path)

    def test_fparam_aparam_preservation(self):
        self.system.to("lmdb", self.lmdb_path)
        system_loaded = dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")

        np.testing.assert_allclose(
            system_loaded.data["fparam"], self.system.data["fparam"]
        )
        np.testing.assert_allclose(
            system_loaded.data["aparam"], self.system.data["aparam"]
        )

    def test_fparam_aparam_auto_registration(self):
        # Save with fparam/aparam registered
        self.system.to("lmdb", self.lmdb_path)

        # Simulate a clean session by restoring original DTYPES
        dpdata.System.DTYPES = self.original_system_dtypes
        dpdata.LabeledSystem.DTYPES = self.original_labeled_system_dtypes

        # Load from LMDB
        system_loaded = dpdata.LabeledSystem(self.lmdb_path, fmt="lmdb")

        # Verify auto-registration and data correctness
        self.assertIn("fparam", [dt.name for dt in dpdata.LabeledSystem.DTYPES])
        self.assertIn("aparam", [dt.name for dt in dpdata.LabeledSystem.DTYPES])

        np.testing.assert_allclose(
            system_loaded.data["fparam"], self.system.data["fparam"]
        )
        np.testing.assert_allclose(
            system_loaded.data["aparam"], self.system.data["aparam"]
        )


if __name__ == "__main__":
    unittest.main()
